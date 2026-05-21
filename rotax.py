#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rotax.py

Single-file version of the rotax rotavirus genotyping pipeline.

Run from the parent folder that contains one folder per sample.
Each sample folder should contain paired-end FASTQ files matching:
    *_R1_001.fastq.gz
    *_R2_001.fastq.gz

Typical use:
    conda activate rotax
    python3 /path/to/rotax.py --blast-db /path/to/241130_rotadb

Optional host database:
    python3 /path/to/rotax.py \
        --blast-db /path/to/241130_rotadb \
        --host-db /path/to/hostdb/hostdatabase250421

If contigs500.fasta already exists in every sample folder and you only want to rerun
genotyping:
    python3 /path/to/rotax.py --skip-assembly --blast-db /path/to/241130_rotadb

Main outputs:
    Parent folder:
        rotax_all_samples_genotyping_report.csv
        rotax_pipeline_log.txt

    Each sample folder:
        contigs.fasta
        contigs500.fasta
        contigs500_COV.fasta
        rotavirus_blast_evaluated.csv
        rotavirus_genotype_hit_summary.csv
        rotax_sample_genotyping_report.csv
        rotavirus_selected_contigs.fasta
        rotavirus_selected_orfs.fasta
        host_blast_comparison.csv, if --host-db is available

Notes:
    - The BLAST database remains external.
    - The old adapter FASTA is embedded here so adapters.fa is no longer required.
    - This script intentionally omits old/commented/obsolete parts from the repository.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path


ADAPTERS_FASTA = """>PrefixPE/1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>PrefixPE/2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
"""

WORK_FASTAS = {
    "work1": "_bbnorm100rnaviralspades",
    "work2": "_bbnorm500rnaviralspades",
    "work3": "_rnaviralspades",
    "work6": "_spadesIsolateAutocutoff",
    "work7": "_clumpifyrnaviralspades",
}

GENE_ORDER = ["G", "P", "I", "R", "C", "M", "A", "N", "T", "E", "H"]

BASE_THRESHOLDS = {
    "G": 80,
    "P": 80,
    "I": 85,
    "R": 83,
    "C": 84,
    "M": 81,
    "A": 79,
    "N": 85,
    "T": 85,
    "E": 85,
    "H": 91,
}

MIN_LENGTHS = {
    "G": 500,
    "P": 1163,
    "I": 596,
    "R": 1632,
    "C": 1319,
    "M": 1253,
    "A": 740,
    "N": 500,
    "T": 500,
    "E": 500,
    "H": 500,
}


def log(message: str, log_path: Path) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"{stamp} - {message}"
    print(line, flush=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def run_command(
    command: list[str],
    log_path: Path,
    cwd: Path | None = None,
    check: bool = False,
) -> bool:
    log("Running: " + " ".join(command), log_path)
    result = subprocess.run(command, cwd=str(cwd) if cwd else None)
    if result.returncode != 0:
        log(f"Command failed with exit code {result.returncode}: {' '.join(command)}", log_path)
        if check:
            raise subprocess.CalledProcessError(result.returncode, command)
        return False
    return True


def is_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


def write_adapter_file(sample_dir: Path) -> Path:
    adapter_path = sample_dir / "_rotax_embedded_adapters.fa"
    adapter_path.write_text(ADAPTERS_FASTA, encoding="utf-8")
    return adapter_path


def find_fastq_pair(sample_dir: Path) -> tuple[Path | None, Path | None]:
    r1 = sorted(sample_dir.glob("*_R1_001.fastq.gz"))
    r2 = sorted(sample_dir.glob("*_R2_001.fastq.gz"))
    if not r1 or not r2:
        return None, None
    return r1[0], r2[0]


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    if not path.exists():
        return records

    header: str | None = None
    seq_parts: list[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def write_fasta(records: list[tuple[str, str]], path: Path, width: int = 70) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(header.rstrip() + "\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i:i + width] + "\n")


def fasta_key(header: str) -> str:
    return header.split()[0].lstrip(">")


def get_coverage(text: str) -> float:
    try:
        return float(text.split("cov_", 1)[1].split("_", 1)[0])
    except (IndexError, ValueError):
        return 0.0


def reverse_complement(seq: str) -> str:
    table = str.maketrans(
        "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
        "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
    )
    return seq.translate(table)[::-1]


def collect_and_filter_contigs(sample_dir: Path, log_path: Path) -> None:
    contigs: list[tuple[str, str]] = []

    for work_dir, suffix in WORK_FASTAS.items():
        fasta_path = sample_dir / work_dir / "contigs.fasta"
        if not fasta_path.exists():
            log(f"File not found: {fasta_path}", log_path)
            continue

        for header, seq in parse_fasta(fasta_path):
            contigs.append((header + suffix, seq))

    write_fasta(contigs, sample_dir / "contigs.fasta")

    long_records = [(h, s) for h, s in contigs if len(s) >= 500]
    write_fasta(long_records, sample_dir / "contigs500.fasta")
    log(f"Wrote {len(contigs)} total contigs and {len(long_records)} contigs >=500 bp", log_path)


def filter_by_coverage(sample_dir: Path, min_cov: float, log_path: Path) -> None:
    input_fasta = sample_dir / "contigs500.fasta"
    output_fasta = sample_dir / "contigs500_COV.fasta"

    records = parse_fasta(input_fasta)
    kept = [(h, s) for h, s in records if get_coverage(h) >= min_cov]
    write_fasta(kept, output_fasta)
    log(f"Coverage filter: kept {len(kept)} of {len(records)} contigs with coverage >= {min_cov}", log_path)


def run_assembly(sample_dir: Path, threads: int, log_path: Path) -> bool:
    r1, r2 = find_fastq_pair(sample_dir)
    if r1 is None or r2 is None:
        log(f"No *_R1_001.fastq.gz / *_R2_001.fastq.gz pair found in {sample_dir.name}; skipping assembly", log_path)
        return False

    adapters = write_adapter_file(sample_dir)

    ok = run_command(
        [
            "trimmomatic", "PE",
            r1.name, r2.name,
            "R1_paired.fastq.gz", "R1_unpaired.fastq.gz",
            "R2_paired.fastq.gz", "R2_unpaired.fastq.gz",
            f"ILLUMINACLIP:{adapters.name}:2:30:10:1",
            "SLIDINGWINDOW:3:15",
            "LEADING:3",
            "TRAILING:3",
            "MINLEN:36",
        ],
        log_path,
        cwd=sample_dir,
    )
    if not ok:
        log(f"Trimmomatic failed in {sample_dir.name}; skipping sample", log_path)
        return False

    if not is_nonempty(sample_dir / "R1_paired.fastq.gz") or not is_nonempty(sample_dir / "R2_paired.fastq.gz"):
        log(f"Trimmomatic paired outputs are empty in {sample_dir.name}; skipping sample", log_path)
        return False

    ok = run_command(
        [
            "clumpify.sh",
            "in1=R1_paired.fastq.gz",
            "in2=R2_paired.fastq.gz",
            "out1=deduped_R1.fastq.gz",
            "out2=deduped_R2.fastq.gz",
            "dedupe",
        ],
        log_path,
        cwd=sample_dir,
    )
    if not ok:
        log(f"Clumpify failed in {sample_dir.name}; continuing without work7", log_path)

    ok = run_command(
        [
            "bbnorm.sh",
            "in1=R1_paired.fastq.gz",
            "in2=R2_paired.fastq.gz",
            "out1=R1bbnorm.fastq.gz",
            "out2=R2bbnorm.fastq.gz",
            "target=100",
            "min=5",
        ],
        log_path,
        cwd=sample_dir,
    )
    if ok and is_nonempty(sample_dir / "R1bbnorm.fastq.gz") and is_nonempty(sample_dir / "R2bbnorm.fastq.gz"):
        run_command(
            [
                "rnaviralspades.py",
                "-t", str(threads),
                "-1", "R1bbnorm.fastq.gz",
                "-2", "R2bbnorm.fastq.gz",
                "-o", "work1",
            ],
            log_path,
            cwd=sample_dir,
        )
    else:
        log(f"BBnorm target=100 output missing/empty in {sample_dir.name}; skipping work1", log_path)

    ok = run_command(
        [
            "bbnorm.sh",
            "in1=R1_paired.fastq.gz",
            "in2=R2_paired.fastq.gz",
            "out1=R1bbnorm2.fastq.gz",
            "out2=R2bbnorm2.fastq.gz",
            "target=500",
            "min=5",
        ],
        log_path,
        cwd=sample_dir,
    )
    if ok and is_nonempty(sample_dir / "R1bbnorm2.fastq.gz") and is_nonempty(sample_dir / "R2bbnorm2.fastq.gz"):
        run_command(
            [
                "rnaviralspades.py",
                "-t", str(threads),
                "-1", "R1bbnorm2.fastq.gz",
                "-2", "R2bbnorm2.fastq.gz",
                "-o", "work2",
            ],
            log_path,
            cwd=sample_dir,
        )
    else:
        log(f"BBnorm target=500 output missing/empty in {sample_dir.name}; skipping work2", log_path)

    run_command(
        [
            "rnaviralspades.py",
            "-t", str(threads),
            "-1", "R1_paired.fastq.gz",
            "-2", "R2_paired.fastq.gz",
            "-o", "work3",
        ],
        log_path,
        cwd=sample_dir,
    )

    run_command(
        [
            "spades.py",
            "-t", str(threads),
            "--isolate",
            "--cov-cutoff", "auto",
            "-1", "R1_paired.fastq.gz",
            "-2", "R2_paired.fastq.gz",
            "-o", "work6",
        ],
        log_path,
        cwd=sample_dir,
    )

    if is_nonempty(sample_dir / "deduped_R1.fastq.gz") and is_nonempty(sample_dir / "deduped_R2.fastq.gz"):
        run_command(
            [
                "rnaviralspades.py",
                "-t", str(threads),
                "-1", "deduped_R1.fastq.gz",
                "-2", "deduped_R2.fastq.gz",
                "-o", "work7",
            ],
            log_path,
            cwd=sample_dir,
        )
    else:
        log(f"Deduped output missing/empty in {sample_dir.name}; skipping work7", log_path)

    collect_and_filter_contigs(sample_dir, log_path)
    return True


def blast_rotavirus(sample_dir: Path, blast_db: str, log_path: Path) -> bool:
    input_fasta = sample_dir / "contigs500_COV.fasta"
    blast_output = sample_dir / "rotavirus_blast_raw.tsv"

    if not input_fasta.exists():
        log(f"{input_fasta} not found; cannot run rotavirus BLAST", log_path)
        return False
    if not parse_fasta(input_fasta):
        log(f"{input_fasta} contains no sequences; cannot run rotavirus BLAST", log_path)
        return False

    return run_command(
        [
            "blastn",
            "-query", str(input_fasta.name),
            "-db", blast_db,
            "-out", str(blast_output.name),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand slen",
            "-max_target_seqs", "1",
            "-strand", "both",
        ],
        log_path,
        cwd=sample_dir,
    )


def evaluate_blast(sample_dir: Path, log_path: Path) -> None:
    input_file = sample_dir / "rotavirus_blast_raw.tsv"
    output_file = sample_dir / "rotavirus_blast_evaluated.csv"
    rejection_report = sample_dir / "rotavirus_rejected_typing_hits.csv"

    blast_headers = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "sstrand", "slen",
    ]
    blast_header_line = "\t".join(blast_headers)

    report_fields = [
        "Folder Name",
        "RejectedEvaluation",
        "Reason",
        "ParsedGenotype",
        "ParsedGene",
        "ParsedPercentIdentity",
        "ParsedAlignmentLength",
        "ParsedSubjectLength",
        "IsFullLength",
        "KnownGeneCriterion",
        "KnownGeneCriterionDetail",
        "KnownGeneCriterionPassed",
        "LengthCriterion",
        "LengthCriterionDetail",
        "LengthCriterionPassed",
        "IdentityCriterion",
        "IdentityCriterionDetail",
        "IdentityCriterionPassed",
        "AcceptedCriterion",
        "AcceptedCriterionDetail",
        "AcceptedCriterionPassed",
        "BlastHeader",
        "BlastLine",
    ] + blast_headers

    def write_rejection(report_writer: csv.DictWriter, row: list[str], reason: str) -> None:
        values = {h: row[i] if i < len(row) else "" for i, h in enumerate(blast_headers)}
        sseqid = values.get("sseqid", "")
        genotype = sseqid.split("|", 1)[0] if sseqid else ""
        gene = genotype[0] if genotype else ""

        try:
            pident = float(values.get("pident", ""))
            pident_text = str(pident)
            pident_ok = True
        except ValueError:
            pident = None
            pident_text = values.get("pident", "")
            pident_ok = False

        try:
            length = int(values.get("length", ""))
            length_text = str(length)
            length_ok = True
        except ValueError:
            length = None
            length_text = values.get("length", "")
            length_ok = False

        try:
            slen = int(values.get("slen", ""))
            slen_text = str(slen)
            slen_ok = True
        except ValueError:
            slen = None
            slen_text = values.get("slen", "")
            slen_ok = False

        known_gene_passed = gene in BASE_THRESHOLDS
        known_gene_detail = (
            f"Gene must be one of {', '.join(GENE_ORDER)}; parsed gene={gene or 'NONE'}."
        )

        is_full = bool(length_ok and slen_ok and length is not None and slen is not None and length >= slen)

        if known_gene_passed and length_ok and slen_ok and length is not None and slen is not None:
            min_len = MIN_LENGTHS[gene]
            length_passed = is_full or length >= min_len
            length_detail = (
                f"Alignment is accepted for length if full-length (length >= slen) "
                f"or length >= minimum for {gene} ({min_len} bp). "
                f"Observed length={length}; slen={slen}; full_length={is_full}."
            )
        else:
            length_passed = False
            if not known_gene_passed:
                length_detail = "Cannot evaluate length threshold because gene is unknown."
            else:
                length_detail = (
                    f"Cannot evaluate length threshold because length/slen could not be parsed. "
                    f"Observed length={length_text}; slen={slen_text}."
                )

        if known_gene_passed and length_ok and slen_ok and pident_ok and pident is not None:
            threshold = BASE_THRESHOLDS[gene] if is_full else BASE_THRESHOLDS[gene] + 2
            identity_passed = pident >= threshold
            identity_detail = (
                f"Percent identity must be >= {threshold:.1f}% for {gene} "
                f"({'full-length hit' if is_full else 'partial hit: base threshold + 2'}). "
                f"Observed pident={pident:.3f}%."
            )
        else:
            identity_passed = False
            if not known_gene_passed:
                identity_detail = "Cannot evaluate percent-identity threshold because gene is unknown."
            else:
                identity_detail = (
                    f"Cannot evaluate percent-identity threshold because pident/length/slen could not be parsed. "
                    f"Observed pident={pident_text}; length={length_text}; slen={slen_text}."
                )

        accepted_passed = bool(known_gene_passed and length_passed and identity_passed)
        accepted_detail = (
            "Accepted only if known-gene, length, and percent-identity criteria all pass. "
            f"known_gene={known_gene_passed}; length={length_passed}; identity={identity_passed}."
        )

        out = {
            "Folder Name": sample_dir.name,
            "RejectedEvaluation": "not accepted",
            "Reason": reason,
            "ParsedGenotype": genotype,
            "ParsedGene": gene,
            "ParsedPercentIdentity": pident_text,
            "ParsedAlignmentLength": length_text,
            "ParsedSubjectLength": slen_text,
            "IsFullLength": str(is_full),
            "KnownGeneCriterion": "gene in accepted rotavirus segment list",
            "KnownGeneCriterionDetail": known_gene_detail,
            "KnownGeneCriterionPassed": str(known_gene_passed),
            "LengthCriterion": "length >= slen OR length >= segment-specific minimum",
            "LengthCriterionDetail": length_detail,
            "LengthCriterionPassed": str(length_passed),
            "IdentityCriterion": "pident >= segment-specific threshold",
            "IdentityCriterionDetail": identity_detail,
            "IdentityCriterionPassed": str(identity_passed),
            "AcceptedCriterion": "all criteria must pass",
            "AcceptedCriterionDetail": accepted_detail,
            "AcceptedCriterionPassed": str(accepted_passed),
            "BlastHeader": blast_header_line,
            "BlastLine": "\t".join(row),
        }
        out.update(values)
        report_writer.writerow(out)

    with input_file.open("r", encoding="utf-8", errors="replace") as infile, \
            output_file.open("w", newline="", encoding="utf-8") as outfile, \
            rejection_report.open("w", newline="", encoding="utf-8") as report_out:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile)
        report_writer = csv.DictWriter(report_out, fieldnames=report_fields)
        report_writer.writeheader()

        writer.writerow(blast_headers + ["Evaluation"])

        for row in reader:
            if len(row) < 14:
                write_rejection(report_writer, row, "BLAST row has fewer than 14 columns and cannot be evaluated")
                continue

            sseqid = row[1]
            try:
                pident = float(row[2])
                length = int(row[3])
                slen = int(row[13])
            except ValueError:
                evaluation = "not accepted"
                writer.writerow(row[:14] + [evaluation])
                write_rejection(report_writer, row, "pident, alignment length, or subject length could not be parsed")
                continue

            genotype = sseqid.split("|")[0]
            if not genotype:
                evaluation = "not accepted"
                writer.writerow(row[:14] + [evaluation])
                write_rejection(report_writer, row, "No genotype prefix found before | in sseqid")
                continue

            gene = genotype[0]
            if gene not in BASE_THRESHOLDS:
                evaluation = "not accepted"
                writer.writerow(row[:14] + [evaluation])
                write_rejection(report_writer, row, "Unknown segment/gene prefix")
                continue

            is_full = length >= slen
            threshold = BASE_THRESHOLDS[gene] if is_full else BASE_THRESHOLDS[gene] + 2

            length_passed = is_full or length >= MIN_LENGTHS[gene]
            identity_passed = pident >= threshold

            if length_passed and identity_passed:
                evaluation = genotype
            else:
                # Below-threshold or too-short hits must not be carried forward as genotypes.
                # The old multi-script pipeline sometimes created values like
                # E12000084.466 for a weak E12 hit at 84.466% identity; downstream
                # digit-stripping then turned that into E12000084466.
                evaluation = "not accepted"
                failed = []
                if not length_passed:
                    failed.append("length criterion failed")
                if not identity_passed:
                    failed.append("percent-identity criterion failed")
                write_rejection(report_writer, row, "; ".join(failed))

            writer.writerow(row[:14] + [evaluation])

    log(f"Generated {output_file.name}", log_path)
    log(f"Generated {rejection_report.name}", log_path)


def summarize_blast(sample_dir: Path, log_path: Path) -> None:
    input_file = sample_dir / "rotavirus_blast_evaluated.csv"
    output_file = sample_dir / "rotavirus_genotype_hit_summary.csv"

    genotype_counts = defaultdict(lambda: {
        "count": 0,
        "full": 0,
        "partial": 0,
        "high_cov_full": 0,
        "low_cov_full": 0,
        "high_cov_partial": 0,
        "low_cov_partial": 0,
    })

    with input_file.open("r", encoding="utf-8", errors="replace") as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            evaluation = row.get("Evaluation", "")
            if evaluation == "not accepted" or not evaluation:
                continue

            try:
                length = int(row["length"])
                slen = int(row["slen"])
            except ValueError:
                continue

            cov = get_coverage(row.get("qseqid", ""))
            is_full = length >= slen

            genotype_counts[evaluation]["count"] += 1
            if is_full:
                genotype_counts[evaluation]["full"] += 1
                if cov >= 10:
                    genotype_counts[evaluation]["high_cov_full"] += 1
                else:
                    genotype_counts[evaluation]["low_cov_full"] += 1
            else:
                genotype_counts[evaluation]["partial"] += 1
                if cov >= 10:
                    genotype_counts[evaluation]["high_cov_partial"] += 1
                else:
                    genotype_counts[evaluation]["low_cov_partial"] += 1

    with output_file.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow([
            "Genotype",
            "Total Count",
            "Full Count",
            "Partial Count",
            "High COV Count Full",
            "Low COV Count Full",
            "High COV Count Partial",
            "Low COV Count Partial",
        ])
        for genotype, counts in sorted(genotype_counts.items()):
            writer.writerow([
                genotype,
                counts["count"],
                counts["full"],
                counts["partial"],
                counts["high_cov_full"],
                counts["low_cov_full"],
                counts["high_cov_partial"],
                counts["low_cov_partial"],
            ])

    log(f"Generated {output_file.name}", log_path)


def extract_gene_and_number(genotype: str) -> tuple[str | None, str]:
    """Return the segment letter and genotype number from a genotype token.

    Handles clean tokens such as E12, GX, G?, P11, and P[11]. It also protects
    against legacy diagnostic strings such as E12000084.466 by cutting at the
    old 0000 low-identity marker before parsing.
    """
    genotype = (genotype or "").strip()
    genotype = genotype.split("0000", 1)[0]

    bracket_match = re.match(r"^([A-Za-z])\[(\d+)\]$", genotype)
    if bracket_match:
        return bracket_match.group(1), bracket_match.group(2)

    match = re.match(r"^([A-Za-z])(\d+|X|\?)$", genotype)
    if not match:
        return None, ""
    return match.group(1), match.group(2)


def format_gene_genotype(gene: str, number: str) -> str:
    if gene == "P" and number not in {"X", "?"}:
        return f"P[{number}]"
    return f"{gene}{number}"


def genotype_summary(sample_dir: Path, log_path: Path) -> None:
    input_file = sample_dir / "rotavirus_genotype_hit_summary.csv"
    output_file = sample_dir / "rotavirus_genotype_constellation_intermediate.csv"

    genotype_data: dict[str, list[dict[str, int | str]]] = {gene: [] for gene in GENE_ORDER}

    with input_file.open("r", newline="", encoding="utf-8", errors="replace") as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            genotype = row.get("Genotype", "")
            gene, number = extract_gene_and_number(genotype)
            if gene not in genotype_data or not number:
                continue

            try:
                item = {
                    "genotype": number,
                    "full_count": int(row.get("Full Count", "0")),
                    "partial_count": int(row.get("Partial Count", "0")),
                    "high_cov_full": int(row.get("High COV Count Full", "0")),
                    "low_cov_full": int(row.get("Low COV Count Full", "0")),
                    "high_cov_partial": int(row.get("High COV Count Partial", "0")),
                    "low_cov_partial": int(row.get("Low COV Count Partial", "0")),
                }
            except ValueError:
                continue

            genotype_data[gene].append(item)

    final_genotype: list[str] = []
    extra_info: list[str] = []

    for gene in GENE_ORDER:
        items = genotype_data[gene]
        if not items:
            final_genotype.append(f"{gene}X")
            continue

        top = items[:]
        for key in ["high_cov_full", "low_cov_full", "high_cov_partial", "low_cov_partial"]:
            max_value = max(int(item[key]) for item in top)
            top = [item for item in top if int(item[key]) == max_value]
            if len(top) == 1:
                break

        if len(top) == 1:
            selected_number = str(top[0]["genotype"])
            final_genotype.append(format_gene_genotype(gene, selected_number))
            selected_numbers = {selected_number}
        else:
            final_genotype.append(f"{gene}?")
            selected_numbers = {str(item["genotype"]) for item in top}
            extra_info.extend(f"{gene}{item['genotype']}" for item in top)

        extra_info.extend(
            f"{gene}{item['genotype']}"
            for item in items
            if str(item["genotype"]) not in selected_numbers
        )

    def sort_key(genotype_str: str) -> tuple[int, int, str]:
        gene = genotype_str[0] if genotype_str else ""
        number_text = "".join(filter(str.isdigit, genotype_str))
        try:
            number = int(number_text)
        except ValueError:
            number = 0
        gene_rank = GENE_ORDER.index(gene) if gene in GENE_ORDER else len(GENE_ORDER)
        return gene_rank, number, genotype_str

    unique_extra_info = sorted(set(extra_info), key=sort_key)
    genotype_string = "-".join(final_genotype)
    extra_info_string = "; ".join(unique_extra_info) if unique_extra_info else "None"

    if all(gt.endswith("X") for gt in final_genotype):
        genotype_string = ""
        extra_info_string = "None"

    with output_file.open("w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["Folder Name", "Genotype", "Extra Information"])
        writer.writerow([sample_dir.name, genotype_string, extra_info_string])

    log(f"Generated {output_file.name}", log_path)


def add_vaccine_match_columns(sample_dir: Path, log_path: Path) -> None:
    """
    Add vaccine-reference percent identity columns to the per-sample genotype report.

    Rotarix is identified by the strain/header marker RIX4414, not by the generic
    word "vaccine". RotaTeq is identified by the marker RotaTeq in the BLAST
    subject header. The values are mean percent identity among accepted BLAST
    hits with the expected genotype/segment calls:

      * Rotarix VP7: accepted G1 hits where sseqid contains RIX4414
      * Rotarix VP4: accepted P8 hits where sseqid contains RIX4414
      * RotaTeq VP6: accepted I2 hits where sseqid contains RotaTeq

    Matching is case-insensitive.
    """
    genotyping_file = sample_dir / "rotavirus_genotype_constellation_intermediate.csv"
    blast_file = sample_dir / "rotavirus_blast_evaluated.csv"
    output_file = sample_dir / "rotax_sample_genotyping_report.csv"

    with genotyping_file.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle))
        genotyping_fields = list(rows[0].keys()) if rows else ["Folder Name", "Genotype", "Extra Information"]

    with blast_file.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        blast_rows = list(csv.DictReader(handle))

    def average_pident(evaluation_value: str, header_marker: str) -> str:
        values: list[float] = []
        marker = header_marker.lower()
        for row in blast_rows:
            if row.get("Evaluation") != evaluation_value:
                continue
            if marker not in row.get("sseqid", "").lower():
                continue
            try:
                values.append(float(row.get("pident", "")))
            except ValueError:
                pass
        if not values:
            return ""
        return str(sum(values) / len(values))

    rotarix_vp7 = average_pident("G1", "RIX4414")
    rotarix_vp4 = average_pident("P8", "RIX4414")
    rotateq_vp6 = average_pident("I2", "RotaTeq")

    extra_fields = ["Rotarix VP7", "Rotarix VP4", "RotaTeq VP6"]
    fields = genotyping_fields + [field for field in extra_fields if field not in genotyping_fields]
    with output_file.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            row["Rotarix VP7"] = rotarix_vp7
            row["Rotarix VP4"] = rotarix_vp4
            row["RotaTeq VP6"] = rotateq_vp6
            writer.writerow(row)

    log(f"Generated {output_file.name}", log_path)


def selected_contigs(sample_dir: Path, log_path: Path) -> None:
    blast_csv = sample_dir / "rotavirus_blast_evaluated.csv"
    input_fasta = sample_dir / "contigs500_COV.fasta"
    output_fasta = sample_dir / "rotavirus_selected_contigs.fasta"

    full_candidates: dict[str, dict[str, float | int | str]] = {}
    partial_candidates: dict[str, dict[str, float | int | str]] = {}

    with blast_csv.open("r", encoding="utf-8", errors="replace") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            genotype = row.get("Evaluation", "")
            if genotype == "not accepted" or not genotype:
                continue

            try:
                length = int(row["length"])
                slen = int(row["slen"])
            except ValueError:
                continue

            qseqid = row["qseqid"]
            cov = get_coverage(qseqid)

            if length >= slen:
                current = full_candidates.get(genotype)
                if current is None:
                    full_candidates[genotype] = {"qseqid": qseqid, "coverage": cov}
                    continue

                current_cov = float(current["coverage"])
                if (cov >= 10 > current_cov) or ((cov >= 10) == (current_cov >= 10) and cov > current_cov):
                    full_candidates[genotype] = {"qseqid": qseqid, "coverage": cov}
            else:
                current = partial_candidates.get(genotype)
                if current is None or length > int(current["length"]):
                    partial_candidates[genotype] = {"qseqid": qseqid, "length": length}

    selected: dict[str, dict[str, str]] = {}
    for genotype in set(full_candidates) | set(partial_candidates):
        if genotype in full_candidates:
            selected[genotype] = {"qseqid": str(full_candidates[genotype]["qseqid"]), "type": "full"}
        else:
            selected[genotype] = {"qseqid": str(partial_candidates[genotype]["qseqid"]), "type": "partial"}

    fasta_records = {fasta_key(h): (h, s) for h, s in parse_fasta(input_fasta)}

    out_records: list[tuple[str, str]] = []
    for genotype in sorted(selected):
        qseqid = selected[genotype]["qseqid"]
        candidate_type = selected[genotype]["type"]
        if qseqid not in fasta_records:
            log(f"Warning: qseqid {qseqid} not found in {input_fasta.name}", log_path)
            continue

        original_header, sequence = fasta_records[qseqid]
        if candidate_type == "partial":
            new_header = f">{genotype}_PARTIAL_{original_header[1:]}"
        else:
            new_header = f">{genotype}_{original_header[1:]}"
        out_records.append((new_header, sequence))

    write_fasta(out_records, output_fasta)
    log(f"Generated {output_fasta.name} with {len(out_records)} sequences", log_path)


def selected_orfs(sample_dir: Path, log_path: Path) -> None:
    blast_csv = sample_dir / "rotavirus_blast_evaluated.csv"
    input_fasta = sample_dir / "contigs500_COV.fasta"
    output_fasta = sample_dir / "rotavirus_selected_orfs.fasta"

    full_candidates: dict[str, dict[str, float | int | str]] = {}
    partial_candidates: dict[str, dict[str, float | int | str]] = {}

    with blast_csv.open("r", encoding="utf-8", errors="replace") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            genotype = row.get("Evaluation", "")
            if genotype == "not accepted" or not genotype:
                continue

            try:
                length = int(row["length"])
                slen = int(row["slen"])
                qstart = int(row["qstart"])
                qend = int(row["qend"])
            except ValueError:
                continue

            qseqid = row["qseqid"]
            cov = get_coverage(qseqid)
            info = {
                "qseqid": qseqid,
                "coverage": cov,
                "length": length,
                "qstart": qstart,
                "qend": qend,
                "sstrand": row.get("sstrand", "+"),
            }

            if length >= slen:
                current = full_candidates.get(genotype)
                if current is None:
                    full_candidates[genotype] = info
                    continue

                current_cov = float(current["coverage"])
                if (cov >= 10 > current_cov) or ((cov >= 10) == (current_cov >= 10) and cov > current_cov):
                    full_candidates[genotype] = info
            else:
                current = partial_candidates.get(genotype)
                if current is None or length > int(current["length"]):
                    partial_candidates[genotype] = info

    selected: dict[str, dict[str, float | int | str]] = {}
    for genotype in set(full_candidates) | set(partial_candidates):
        if genotype in full_candidates:
            selected[genotype] = dict(full_candidates[genotype])
            selected[genotype]["type"] = "full"
        else:
            selected[genotype] = dict(partial_candidates[genotype])
            selected[genotype]["type"] = "partial"

    fasta_records = {fasta_key(h): (h, s) for h, s in parse_fasta(input_fasta)}

    out_records: list[tuple[str, str]] = []
    for genotype in sorted(selected):
        candidate = selected[genotype]
        qseqid = str(candidate["qseqid"])
        if qseqid not in fasta_records:
            log(f"Warning: qseqid {qseqid} not found in {input_fasta.name}", log_path)
            continue

        original_header, full_seq = fasta_records[qseqid]
        qstart = int(candidate["qstart"])
        qend = int(candidate["qend"])
        start = min(qstart, qend)
        end = max(qstart, qend)

        orf_seq = full_seq[start - 1:end]
        if str(candidate["sstrand"]).strip().lower() == "minus":
            orf_seq = reverse_complement(orf_seq)

        if candidate["type"] == "partial":
            new_header = f">{genotype}_PARTIAL_ORF_{original_header[1:]}"
        else:
            new_header = f">{genotype}_ORF_{original_header[1:]}"
        out_records.append((new_header, orf_seq))

    write_fasta(out_records, output_fasta)
    log(f"Generated {output_fasta.name} with {len(out_records)} sequences", log_path)


def run_host_blast(sample_dir: Path, host_db: str | None, log_path: Path) -> None:
    if not host_db:
        log("No --host-db provided; skipping host BLAST", log_path)
        return

    query = sample_dir / "rotavirus_selected_orfs.fasta"
    if not query.exists() or not parse_fasta(query):
        log(f"{query.name} missing/empty; skipping host BLAST", log_path)
        return

    run_command(
        [
            "blastn",
            "-query", query.name,
            "-db", host_db,
            "-out", "host_blast_comparison.csv",
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand slen",
            "-max_target_seqs", "1",
            "-strand", "both",
        ],
        log_path,
        cwd=sample_dir,
    )


def safe_remove(path: Path, log_path: Path) -> None:
    try:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    except OSError as exc:
        log(f"Could not remove {path}: {exc}", log_path)


def cleanup_sample(sample_dir: Path, log_path: Path) -> None:
    for name in ["work1", "work2", "work3", "work6", "work7"]:
        safe_remove(sample_dir / name, log_path)

    for name in [
        "R1bbnorm.fastq.gz",
        "R1bbnorm2.fastq.gz",
        "R2bbnorm.fastq.gz",
        "R2bbnorm2.fastq.gz",
        "R1_paired.fastq.gz",
        "R1_unpaired.fastq.gz",
        "R2_paired.fastq.gz",
        "R2_unpaired.fastq.gz",
        "deduped_R1.fastq.gz",
        "deduped_R2.fastq.gz",
        "rotavirus_genotype_constellation_intermediate.csv",
        "rotavirus_blast_raw.tsv",
        "_rotax_embedded_adapters.fa",
    ]:
        safe_remove(sample_dir / name, log_path)


def process_sample(
    sample_dir: Path,
    args: argparse.Namespace,
    log_path: Path,
) -> bool:
    log(f"Processing sample folder: {sample_dir.name}", log_path)

    try:
        if not args.skip_assembly:
            assembled = run_assembly(sample_dir, args.threads, log_path)
            if not assembled and not (sample_dir / "contigs500.fasta").exists():
                log(f"No assembly output for {sample_dir.name}; skipping genotyping", log_path)
                return False

        if not (sample_dir / "contigs500.fasta").exists():
            log(f"{sample_dir.name}: contigs500.fasta not found; skipping genotyping", log_path)
            return False

        filter_by_coverage(sample_dir, args.min_cov, log_path)

        if not blast_rotavirus(sample_dir, args.blast_db, log_path):
            log(f"{sample_dir.name}: rotavirus BLAST failed or produced no usable input", log_path)
            return False

        if not is_nonempty(sample_dir / "rotavirus_blast_raw.tsv"):
            log(f"{sample_dir.name}: rotavirus_blast_raw.tsv is empty; writing empty genotype report", log_path)
            write_empty_genotype_report(sample_dir)
            return True

        evaluate_blast(sample_dir, log_path)
        summarize_blast(sample_dir, log_path)
        genotype_summary(sample_dir, log_path)
        add_vaccine_match_columns(sample_dir, log_path)
        selected_contigs(sample_dir, log_path)
        selected_orfs(sample_dir, log_path)
        run_host_blast(sample_dir, args.host_db, log_path)

        return True

    except Exception as exc:
        log(f"ERROR in {sample_dir.name}: {exc}", log_path)
        return False

    finally:
        if args.cleanup:
            cleanup_sample(sample_dir, log_path)


def write_empty_genotype_report(sample_dir: Path) -> None:
    output = sample_dir / "rotax_sample_genotyping_report.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Folder Name", "Genotype", "Extra Information", "Rotarix VP7", "Rotarix VP4", "RotaTeq VP6"])
        writer.writerow([sample_dir.name, "", "None", "", "", ""])


def g_p_genotype_from_constellation(genotype: str) -> str:
    """Return the G-P genotype part from a full 11-segment genotype constellation."""
    parts = genotype.split("-") if genotype else []
    return "-".join(parts[:2]) if len(parts) >= 2 else genotype


def infer_result_category(g_p_genotype: str) -> str:
    """Return a generic result category for the final report.

    RVA is used when both G and P are resolved numeric genotypes.
    Otherwise the row is left as an unresolved rotavirus sequence result.
    """
    if "-" not in g_p_genotype:
        return "unresolved_rotavirus_sequence"

    g_part, p_part = g_p_genotype.split("-", 1)
    if re.fullmatch(r"G\d+", g_part) and (
        re.fullmatch(r"P\d+", p_part) or re.fullmatch(r"P\[\d+\]", p_part)
    ):
        return "RVA"

    return "unresolved_rotavirus_sequence"


def combine_sample_reports(parent_dir: Path, log_path: Path) -> Path | None:
    output_file = parent_dir / "rotax_all_samples_genotyping_report.csv"
    report_files = sorted(parent_dir.glob("*/rotax_sample_genotyping_report.csv"))

    if not report_files:
        log("No sample rotax_sample_genotyping_report.csv files found; no final report created", log_path)
        return None

    with output_file.open("w", encoding="utf-8", newline="") as out:
        writer: csv.DictWriter | None = None

        for report in report_files:
            with report.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames is None:
                    continue

                extra_fields = ["G_P_Genotype", "ResultCategory"]
                fieldnames = list(reader.fieldnames) + [
                    field for field in extra_fields if field not in reader.fieldnames
                ]

                if writer is None:
                    writer = csv.DictWriter(out, fieldnames=fieldnames)
                    writer.writeheader()

                for row in reader:
                    g_p_genotype = g_p_genotype_from_constellation(row.get("Genotype", ""))
                    row["G_P_Genotype"] = g_p_genotype
                    row["ResultCategory"] = infer_result_category(g_p_genotype)
                    writer.writerow(row)

    log(f"Generated final all-samples genotyping report: {output_file.name}", log_path)
    return output_file


def sample_folders(parent_dir: Path) -> list[Path]:
    folders: list[Path] = []
    for path in sorted(parent_dir.iterdir()):
        if not path.is_dir():
            continue
        if path.name.startswith("."):
            continue
        folders.append(path)
    return folders


def check_required_tools(skip_assembly: bool, host_db: str | None, log_path: Path) -> None:
    required = ["blastn"]
    if not skip_assembly:
        required.extend(["trimmomatic", "clumpify.sh", "bbnorm.sh", "rnaviralspades.py", "spades.py"])

    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        log("Missing required executable(s) in PATH: " + ", ".join(missing), log_path)
        raise SystemExit(2)

    if host_db and shutil.which("blastn") is None:
        log("blastn missing in PATH; host BLAST cannot run", log_path)
        raise SystemExit(2)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Single-file rotax pipeline for rotavirus genotyping from Illumina reads or existing contigs."
    )
    parser.add_argument(
        "--blast-db",
        default=str(Path(__file__).resolve().parent / "241130_rotadb"),
        help="Path/prefix of the rotavirus BLAST database. Default: 241130_rotadb next to this script.",
    )
    parser.add_argument(
        "--host-db",
        default=None,
        help="Optional path/prefix of the host BLAST database. If omitted, host_blast_comparison.csv is skipped.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=12,
        help="Threads for SPAdes/rnaviralSPAdes. Default: 12.",
    )
    parser.add_argument(
        "--min-cov",
        type=float,
        default=3.0,
        help="Minimum SPAdes k-mer coverage parsed from contig header. Default: 3.",
    )
    parser.add_argument(
        "--skip-assembly",
        action="store_true",
        help="Skip Trimmomatic/BBmap/SPAdes and start from existing contigs500.fasta in each sample folder.",
    )
    parser.add_argument(
        "--no-cleanup",
        dest="cleanup",
        action="store_false",
        help="Keep intermediate work folders and FASTQ files.",
    )
    parser.set_defaults(cleanup=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    parent_dir = Path.cwd()
    log_path = parent_dir / "rotax_pipeline_log.txt"
    start_time = time.time()
    start_date = time.strftime("%Y-%m-%d %H:%M:%S")

    log(f"Script started in {parent_dir}", log_path)
    log(f"Start time: {start_date}", log_path)
    log(f"Rotavirus BLAST DB: {args.blast_db}", log_path)
    if args.host_db:
        log(f"Host BLAST DB: {args.host_db}", log_path)
    else:
        log("Host BLAST DB: not provided", log_path)

    check_required_tools(args.skip_assembly, args.host_db, log_path)

    folders = sample_folders(parent_dir)
    if not folders:
        log("No sample folders found one level below current directory", log_path)
        return 1

    processed = 0
    failed = 0

    for folder in folders:
        ok = process_sample(folder, args, log_path)
        if ok:
            processed += 1
        else:
            failed += 1

    combine_sample_reports(parent_dir, log_path)

    elapsed = int(time.time() - start_time)
    log(f"Script ended in {parent_dir}", log_path)
    log(f"Time spent: {elapsed} seconds", log_path)
    log(f"Folders processed successfully: {processed}", log_path)
    log(f"Folders failed/skipped: {failed}", log_path)

    return 0 if processed > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())