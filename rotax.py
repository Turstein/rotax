#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
rotax.py

Standalone single-file version of the RotaX rotavirus genotyping pipeline.
It contains the complete pipeline and does not import rotax5 or rotax6.

Version:
    RotaX script update: 2026-11-11-v34-report-columns
    New in v34: Equine is renamed G3 equine-like; G_P_Genotype is removed.
    New in v33: LabwareGenotype holds the display call; Genotype holds the
    constellation; ResultCategory is renamed Art; no coinfection comment column.
    New in v32: Koinfeksjon(...) display with a preserved genotype constellation;
    closest-reference RotaTeq percent identities for all 11 segments.
    Update note: Each sample mapping folder contains its own GenBank
    mapping-consensus FASTA with terminal Ns removed, complete and partial
    segments, simple sample_segment IDs, and an explicit
    [completeness=complete|partial] descriptor on every definition line.
    A companion GenBank feature CSV uses coordinates shifted for the
    terminal-N trimming. Complete requires N-free 5-prime and 3-prime
    terminal regions before trimming.

How to verify the installed script:
    python3 /mnt/tempdata/rota/rotax/rotax.py --version
    head -30 /mnt/tempdata/rota/rotax/rotax.py

Run from the parent folder that contains one folder per sample.
Each sample folder should contain paired-end FASTQ files matching:
    *_R1_001.fastq.gz
    *_R2_001.fastq.gz

Default use:
    conda activate rotax
    python3 /mnt/tempdata/rota/rotax/rotax.py

By default, the rotavirus BLAST database is auto-detected in:
    /path/to/blastdb/

The full-segment nucleotide BLAST database is auto-detected in:
    /path/to/fulldb/

You can still provide another database prefix/path:
    python3 /mnt/tempdata/rota/rotax/rotax.py --blast-db /path/to/241130_rotadb

Assembly selection:
    Default is to run all assembly modes.
    To run only one or more selected modes:
        python3 /mnt/tempdata/rota/rotax/rotax.py --assemblies rnaviralspades
        python3 /mnt/tempdata/rota/rotax/rotax.py --assemblies rnaviralspades spades-isolate
        python3 /mnt/tempdata/rota/rotax/rotax.py --assemblies bbnorm100-rnaviralspades bbnorm500-rnaviralspades

Supported assembly mode names:
    bbnorm100-rnaviralspades
    bbnorm500-rnaviralspades
    rnaviralspades
    spades-isolate
    clumpify-rnaviralspades

Skip assembly and reuse an existing contigs500.fasta in each sample folder:
    python3 /mnt/tempdata/rota/rotax/rotax.py --skip-assembly --overwrite-existing

In this mode, overwrite removes old downstream reports/mapping outputs but
preserves contigs500.fasta as the required input.

FASTA-only genotyping:
    To genotype an existing FASTA file without FASTQ assembly:
        python3 /mnt/tempdata/rota/rotax/rotax.py --fasta-input selected_contigs.fasta

    The read-mapping consensus step needs paired FASTQ reads. In --fasta-input
    mode, either add --mapping-r1/--mapping-r2 or disable mapping with --no-mapping:
        python3 /mnt/tempdata/rota/rotax/rotax.py --fasta-input selected_contigs.fasta --mapping-r1 sample_R1_001.fastq.gz --mapping-r2 sample_R2_001.fastq.gz
        python3 /mnt/tempdata/rota/rotax/rotax.py --fasta-input selected_contigs.fasta --no-mapping

Main outputs:
    Parent folder:
        rotax_all_samples_genotyping_report.csv
        rotax_full_length_low_identity_x_report.csv
        rotax_mapping_summary.csv
        rotax_mapping_consensus_features.csv
        rotax_mapping_consensus_all.fasta
        rotax_mapping_consensus_full.fasta
        rotax_pipeline_log.txt

    Each sample folder:
        rotavirus_blast_evaluated.csv
        rotavirus_genotype_hit_summary.csv
        rotax_sample_genotyping_report.csv
        rotavirus_selected_contigs.fasta
        rotavirus_selected_orfs.fasta
        rotavirus_full_length_low_identity_x_hits.csv
        mapping/rotavirus_mapping_consensus.fasta
        mapping/rotavirus_mapping_consensus_full.fasta
        mapping/rotavirus_mapping_consensus_genbank.fasta
        mapping/rotavirus_mapping_consensus_genbank_features.csv
        mapping/rotavirus_mapping_consensus_partial_or_review.fasta
        mapping/rotavirus_mapping_orfs.fasta
        mapping/rotavirus_mapping_references.fasta
        mapping/rotavirus_mapping_summary.csv
        mapping/rotavirus_mapping_consensus_features.csv
        mapping/rotavirus_mapping_depth.tsv
        mapping/rotavirus_mapping.variants.vcf.gz
        mapping/rotavirus_mapping.variants.vcf.gz.csi (or .tbi)
        mapping/rotavirus_mapping_validation.html
        mapping/protein_alignments/*.protein.aln.fasta

    Optional legacy SPAdes-contig full-segment QC outputs are generated only
    with --spades-full-segment-qc. They are not used as GenBank sequences.

Notes:
    - Host database searching has been removed.
    - The old adapter FASTA is embedded here so adapters.fa is no longer required.
    - Mapping consensus v14: non-mixed samples use one selected ORF candidate per segment; possible mixed samples keep multiple ORF candidates per segment. Each selected ORF is BLASTed to the full-segment DB and only its closest same-segment full reference is indexed. Full references are oriented to positive sense before mapping. CDS coordinates on the read-mapping consensus are determined independently from the SPAdes ORF by translating all three forward consensus frames and protein-aligning them to the best same-segment genotyping-database ORF hit. The selected protein alignment is saved. N-masked codons are allowed, and one to three missing or extra amino acids are accepted as a small in-frame difference and documented. Sequences without an annotated CDS are not marked GenBankReady.
"""


from __future__ import annotations

import argparse
import csv
import gzip
import glob
import html
import io
import math
import os
import re
import shutil
import shlex
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path


ROTAX_VERSION = "2026-11-11-v34-report-columns"
ROTAX_UPDATE_NOTE = (
    "Coinfections use Koinfeksjon(...) in LabwareGenotype; Genotype retains the "
    "full constellation; ResultCategory is renamed Art and the coinfection "
    "comment column is removed; all 11 RotaTeq "
    "segment identities require a closest reported reference match. "
    "Combined genotyping report uses G3 equine-like=y for the selected G3 "
    "Equine-like hit and omits the redundant G_P_Genotype column; "
    "BBNorm uses one interleaved input stream; BWA mismatch penalty 2 avoids "
    "compensating indels; deletion-spanning reads count toward depth; IUPAC "
    "requires minor AD >= 5 and minor AD fraction >= 20 percent, otherwise the "
    "AD-major base is called; final variant VCF/AD evidence is retained; the "
    "default full-segment database is auto-detected in fulldb beside this script; "
    "mapping validation is integrated as one combined interactive HTML with "
    "the former PDF alignment/CDS bar figures embedded per segment; each sample "
    "mapping folder contains a terminal-N-trimmed GenBank mapping-consensus FASTA "
    "with complete and partial segments, simple sample_segment IDs, explicit "
    "completeness descriptors, and a companion feature CSV with trimming-adjusted "
    "coordinates; complete requires N-free 5-prime and 3-prime terminal regions "
    "before trimming; "
    "the mapping validation HTML includes an interactive SPAdes-ORF versus "
    "mapping-CDS alignment and side-by-side genotyping-database BLAST results."
)

MAPPING_BWA_MISMATCH_PENALTY = 2
IUPAC_MIN_MINOR_AD = 5
IUPAC_MIN_MINOR_FRACTION = 0.20


ADAPTERS_FASTA = """>PrefixPE/1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>PrefixPE/2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
"""

ASSEMBLY_OPTIONS = {
    "bbnorm100-rnaviralspades": {
        "work_dir": "work1",
        "suffix": "_bbnorm100rnaviralspades",
        "description": "BBnorm target=100, then rnaviralSPAdes",
    },
    "bbnorm500-rnaviralspades": {
        "work_dir": "work2",
        "suffix": "_bbnorm500rnaviralspades",
        "description": "BBnorm target=500, then rnaviralSPAdes",
    },
    "rnaviralspades": {
        "work_dir": "work3",
        "suffix": "_rnaviralspades",
        "description": "rnaviralSPAdes on trimmed paired reads without BBnorm",
    },
    "spades-isolate": {
        "work_dir": "work6",
        "suffix": "_spadesIsolateAutocutoff",
        "description": "SPAdes with --isolate and --cov-cutoff auto",
    },
    "clumpify-rnaviralspades": {
        "work_dir": "work7",
        "suffix": "_clumpifyrnaviralspades",
        "description": "Clumpify deduplication, then rnaviralSPAdes",
    },
}

DEFAULT_ASSEMBLIES = list(ASSEMBLY_OPTIONS)

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

# Segment/genotype mapping used for GenBank full-segment evaluation.
# Headers such as >T28|RVA/... are parsed from the first token before |,
# and the first letter identifies the RVA genome segment.
GENOTYPE_PREFIX_TO_SEGMENT = {
    "G": "VP7",
    "P": "VP4",
    "I": "VP6",
    "R": "VP1",
    "C": "VP2",
    "M": "VP3",
    "A": "NSP1",
    "N": "NSP2",
    "T": "NSP3",
    "E": "NSP4",
    "H": "NSP5",
}

SEGMENT_TO_GENOTYPE_PREFIX = {v: k for k, v in GENOTYPE_PREFIX_TO_SEGMENT.items()}
SEGMENT_ORDER = ["VP7", "VP4", "VP6", "VP1", "VP2", "VP3", "NSP1", "NSP2", "NSP3", "NSP4", "NSP5"]

SAMPLE_GENOTYPING_FIELDS = [
    "Folder Name", "LabwareGenotype", "Genotype", "Extra Information",
    "Rotarix VP7", "Rotarix VP4",
] + [f"RotaTeq {segment}" for segment in SEGMENT_ORDER]
ALL_GENOTYPING_FIELDS = SAMPLE_GENOTYPING_FIELDS + ["Art", "G3 equine-like"]

EXPECTED_ORF_LENGTHS = {
    "VP1": 3264,
    "VP2": 2643,
    "VP3": 2505,
    "VP4": 2328,
    "VP6": 1191,
    "VP7": 978,
    "NSP1": 1473,
    "NSP2": 951,
    "NSP3": 945,
    "NSP4": 525,
    "NSP5": 594,
}

STOP_CODONS = {"TAA", "TAG", "TGA"}

GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Conservative quality gates for protein-guided CDS annotation. Identity is
# calculated only at fully called, non-X amino-acid positions. The comparable
# fraction prevents a mostly N-masked consensus from being annotated from only
# a few matching residues, while still allowing ordinary scattered N codons.
PROTEIN_ORF_MIN_IDENTITY = 0.60
PROTEIN_ORF_MIN_COMPARABLE_FRACTION = 0.70
# One to three amino-acid insertions/deletions are accepted as small in-frame
# differences. Larger protein gaps remain a reason to suppress CDS annotation.
PROTEIN_ORF_MAX_SMALL_INDEL_AA = 3


# Files created by RotaX inside each sample folder. These are used to detect
# whether a sample already has results and to clean old results before an
# overwrite run. Raw input FASTQ files are deliberately not included here.
ROTAX_SAMPLE_OUTPUT_FILES = [
    "contigs.fasta",
    "contigs500.fasta",
    "contigs500_COV.fasta",
    "rotavirus_blast_raw.tsv",
    "rotavirus_blast_evaluated.csv",
    "rotavirus_rejected_typing_hits.csv",
    "rotavirus_genotype_hit_summary.csv",
    "rotavirus_genotype_constellation_intermediate.csv",
    "rotax_sample_genotyping_report.csv",
    "rotavirus_selected_contigs.fasta",
    "rotavirus_selected_orfs.fasta",
    "rotavirus_full_segment_blast_raw.tsv",
    "rotavirus_genbank_full_segment_query.fasta",
    "rotavirus_genbank_full_segment_report.csv",
    "rotavirus_genbank_feature_coordinates.csv",
    "rotavirus_genbank_clean_full_segments.fasta",
    "rotavirus_genbank_partial_or_review_segments.fasta",
    "rotavirus_full_length_low_identity_x_hits.csv",
    "alignments",
    "mapping",
    "_rotax_embedded_adapters.fa",
    "R1bbnorm.fastq.gz",
    "R1bbnorm2.fastq.gz",
    "R2bbnorm.fastq.gz",
    "R2bbnorm2.fastq.gz",
    "R1_paired.fastq.gz",
    "R1_paired.repaired.fastq.gz",
    "R1_paired.before_repair.fastq.gz",
    "R1_unpaired.fastq.gz",
    "R2_paired.fastq.gz",
    "R2_paired.repaired.fastq.gz",
    "R2_paired.before_repair.fastq.gz",
    "R2_unpaired.fastq.gz",
    "R1R2_singletons.fastq.gz",
    "deduped_R1.fastq.gz",
    "deduped_R2.fastq.gz",
]


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


def find_fastq_pairs(sample_dir: Path) -> list[tuple[Path, Path]]:
    """Find exact R1/R2 FASTQ pairs without mixing unrelated files.

    Earlier RotaX versions selected the first sorted *_R1_001.fastq.gz and the
    first sorted *_R2_001.fastq.gz. That is unsafe when a folder contains more
    than one pair because BBTools/Trimmomatic will then see unsynchronised reads.
    """
    pairs: list[tuple[Path, Path]] = []
    for r1 in sorted(sample_dir.glob("*_R1_001.fastq.gz")):
        r2 = sample_dir / r1.name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
        if r2.exists():
            pairs.append((r1, r2))
    return pairs


def find_fastq_pair(sample_dir: Path) -> tuple[Path | None, Path | None]:
    pairs = find_fastq_pairs(sample_dir)
    if not pairs:
        return None, None
    return pairs[0]


def fastq_record_count(path: Path) -> int:
    """Count FASTQ records in plain or gzipped FASTQ. Return -1 on failure."""
    import gzip

    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        line_count = 0
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line_count, _ in enumerate(handle, start=1):
                pass
        if line_count % 4 != 0:
            return -1
        return line_count // 4
    except OSError:
        return -1



def normalize_fastq_read_id(header_line: str) -> str:
    """Return a mate-independent FASTQ read ID for paired-read validation.

    Examples handled:
      @INST:RUN:FLOWCELL:LANE:TILE:X:Y 1:N:0:INDEX
      @INST:RUN:FLOWCELL:LANE:TILE:X:Y 2:N:0:INDEX
      @READ/1 and @READ/2
    """
    h = header_line.strip()
    if h.startswith("@"):
        h = h[1:]
    first_token = h.split()[0] if h else ""
    if first_token.endswith("/1") or first_token.endswith("/2"):
        first_token = first_token[:-2]
    return first_token


def fastq_pair_id_check(r1_path: Path, r2_path: Path) -> tuple[bool, int, str]:
    """Stream two FASTQ files and confirm mate IDs match in the same order.

    Returns (ok, pairs_checked, message). This deliberately reads the full files;
    it is slower than checking the first few records, but prevents late BBTools
    PairStreamer crashes caused by one bad/misaligned block far into the file.
    """
    import gzip

    opener1 = gzip.open if str(r1_path).endswith(".gz") else open
    opener2 = gzip.open if str(r2_path).endswith(".gz") else open
    pairs = 0
    try:
        with opener1(r1_path, "rt", encoding="utf-8", errors="replace") as h1, \
                opener2(r2_path, "rt", encoding="utf-8", errors="replace") as h2:
            while True:
                r1_header = h1.readline()
                r2_header = h2.readline()
                if not r1_header and not r2_header:
                    return True, pairs, "all paired read IDs matched"
                if not r1_header or not r2_header:
                    return False, pairs, "one FASTQ ended before the other"

                r1_seq = h1.readline()
                r1_plus = h1.readline()
                r1_qual = h1.readline()
                r2_seq = h2.readline()
                r2_plus = h2.readline()
                r2_qual = h2.readline()
                if not (r1_seq and r1_plus and r1_qual and r2_seq and r2_plus and r2_qual):
                    return False, pairs, "incomplete FASTQ record encountered"

                pairs += 1
                id1 = normalize_fastq_read_id(r1_header)
                id2 = normalize_fastq_read_id(r2_header)
                if id1 != id2:
                    return (
                        False,
                        pairs,
                        f"first mismatching pair at record {pairs}: R1={id1}; R2={id2}",
                    )
    except OSError as exc:
        return False, pairs, f"could not read FASTQ files: {exc}"


def repair_paired_fastqs(
    sample_dir: Path,
    in1: str,
    in2: str,
    out1: str,
    out2: str,
    outs: str,
    log_path: Path,
) -> bool:
    """Synchronise paired FASTQs with BBTools repair.sh.

    This prevents BBTools PairStreamer errors such as:
      AssertionError: List size mismatch
      AssertionError: <read1> != <read2>
    """
    if shutil.which("repair.sh") is None:
        log("repair.sh not found in PATH; cannot repair/synchronise paired FASTQ files", log_path)
        return False

    ok = run_command(
        [
            "repair.sh",
            f"in1={in1}",
            f"in2={in2}",
            f"out1={out1}",
            f"out2={out2}",
            f"outs={outs}",
            "repair=t",
            "overwrite=t",
        ],
        log_path,
        cwd=sample_dir,
    )
    if not ok:
        log(f"repair.sh failed for {sample_dir.name}", log_path)
        return False

    out1_path = sample_dir / out1
    out2_path = sample_dir / out2
    if not is_nonempty(out1_path) or not is_nonempty(out2_path):
        log(
            f"repair.sh produced empty repaired paired FASTQs in {sample_dir.name}; "
            "the selected R1/R2 files are probably not matching mates.",
            log_path,
        )
        return False

    n1 = fastq_record_count(out1_path)
    n2 = fastq_record_count(out2_path)
    if n1 >= 0 and n2 >= 0:
        log(f"FASTQ repair kept {n1} synchronized paired reads in {sample_dir.name}", log_path)
        if n1 != n2:
            log(f"FASTQ repair output still has unequal read counts: R1={n1}, R2={n2}", log_path)
            return False
    else:
        log("FASTQ repair completed, but record counts could not be verified", log_path)

    return True


def run_bbnorm_interleaved(
    sample_dir: Path,
    interleaved_input: Path,
    target: int,
    out1_name: str,
    out2_name: str,
    bbnorm_threads: int,
    temp_dir: Path,
    log_path: Path,
) -> bool:
    """Normalize one interleaved stream, then restore synchronized R1/R2.

    BBTools 39.81 can raise ``PairStreamer`` list-size assertions when BBNorm
    reads two synchronized gzip files concurrently.  A single interleaved
    input avoids that dual-stream race.
    """
    normalized_interleaved = temp_dir / f"target{target}.interleaved.fastq.gz"
    ok = run_command(
        [
            "bbnorm.sh",
            f"in={interleaved_input}",
            f"out={normalized_interleaved}",
            f"target={target}",
            "min=5",
            f"t={bbnorm_threads}",
            "interleaved=t",
            "overwrite=t",
        ],
        log_path,
        cwd=sample_dir,
    )
    if not ok or not is_nonempty(normalized_interleaved):
        return False

    ok = run_command(
        [
            "reformat.sh",
            f"in={normalized_interleaved}",
            f"out1={out1_name}",
            f"out2={out2_name}",
            "interleaved=t",
            "overwrite=t",
        ],
        log_path,
        cwd=sample_dir,
    )
    out1_path = sample_dir / out1_name
    out2_path = sample_dir / out2_name
    if not ok or not is_nonempty(out1_path) or not is_nonempty(out2_path):
        safe_remove(out1_path, log_path)
        safe_remove(out2_path, log_path)
        return False

    pair_ok, pair_count, pair_message = fastq_pair_id_check(out1_path, out2_path)
    if not pair_ok:
        log(
            f"BBNorm target={target} output pair validation failed for "
            f"{out1_name} + {out2_name}; checked {pair_count} pair(s); "
            f"{pair_message}",
            log_path,
        )
        safe_remove(out1_path, log_path)
        safe_remove(out2_path, log_path)
        return False

    log(
        f"BBNorm target={target} interleaved workaround produced "
        f"{pair_count} synchronized pair(s) in {out1_name} + {out2_name}",
        log_path,
    )
    return True


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


def collect_and_filter_contigs(
    sample_dir: Path,
    log_path: Path,
    assemblies: list[str] | None = None,
) -> None:
    contigs: list[tuple[str, str]] = []
    selected = assemblies or DEFAULT_ASSEMBLIES

    for assembly_name in selected:
        info = ASSEMBLY_OPTIONS[assembly_name]
        work_dir = str(info["work_dir"])
        suffix = str(info["suffix"])
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


def run_assembly(
    sample_dir: Path,
    threads: int,
    log_path: Path,
    assemblies: list[str],
    fastq_repair: bool = True,
    bbnorm_threads: int = 1,
) -> bool:
    fastq_pairs = find_fastq_pairs(sample_dir)
    if not fastq_pairs:
        log(f"No exact *_R1_001.fastq.gz / *_R2_001.fastq.gz pair found in {sample_dir.name}; skipping assembly", log_path)
        return False
    if len(fastq_pairs) > 1:
        pair_text = "; ".join(f"{a.name} + {b.name}" for a, b in fastq_pairs)
        log(
            f"More than one exact FASTQ pair found in {sample_dir.name}; RotaX expects one pair per sample folder. "
            f"Move/merge files or use one folder per pair. Found: {pair_text}",
            log_path,
        )
        return False

    r1, r2 = fastq_pairs[0]
    log(f"Selected FASTQ pair for {sample_dir.name}: {r1.name} + {r2.name}", log_path)

    selected = set(assemblies)
    log(f"BBnorm thread setting for {sample_dir.name}: t={bbnorm_threads}", log_path)
    adapters = write_adapter_file(sample_dir)

    ok = run_command(
        [
            "trimmomatic", "PE",
            r1.name, r2.name,
            "R1_paired.fastq.gz", "R1_unpaired.fastq.gz",
            "R2_paired.fastq.gz", "R2_unpaired.fastq.gz",
            f"ILLUMINACLIP:{adapters.name}:2:30:10:1",
            "SLIDINGWINDOW:4:20",
            "MINLEN:50",
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

    assembly_r1 = "R1_paired.fastq.gz"
    assembly_r2 = "R2_paired.fastq.gz"

    if fastq_repair:
        ok = repair_paired_fastqs(
            sample_dir,
            "R1_paired.fastq.gz",
            "R2_paired.fastq.gz",
            "R1_paired.repaired.fastq.gz",
            "R2_paired.repaired.fastq.gz",
            "R1R2_singletons.fastq.gz",
            log_path,
        )
        if not ok:
            log(
                f"Paired FASTQ synchronization failed in {sample_dir.name}; skipping assembly. "
                "This usually means the R1 and R2 files are not matching mates.",
                log_path,
            )
            return False
        assembly_r1 = "R1_paired.repaired.fastq.gz"
        assembly_r2 = "R2_paired.repaired.fastq.gz"
        log(
            f"Using repair.sh-synchronised paired FASTQs for downstream assembly in {sample_dir.name}: "
            f"{assembly_r1} + {assembly_r2}",
            log_path,
        )
    else:
        log(
            f"FASTQ repair is disabled for {sample_dir.name}; validating Trimmomatic paired output before assembly",
            log_path,
        )

    pair_ok, pair_count, pair_message = fastq_pair_id_check(
        sample_dir / assembly_r1,
        sample_dir / assembly_r2,
    )
    if not pair_ok:
        log(
            f"Paired FASTQ ID validation failed in {sample_dir.name} for {assembly_r1} + {assembly_r2}; "
            f"checked {pair_count} pair(s); {pair_message}. Skipping assembly to avoid BBTools/SPAdes pair errors.",
            log_path,
        )
        return False
    log(
        f"Paired FASTQ ID validation passed in {sample_dir.name}: {pair_count} synchronized pair(s) "
        f"in {assembly_r1} + {assembly_r2}",
        log_path,
    )

    def validate_current_pair_before_step(step_name: str, r1_name: str, r2_name: str) -> bool:
        ok_pair, checked, msg = fastq_pair_id_check(sample_dir / r1_name, sample_dir / r2_name)
        if ok_pair:
            log(f"{step_name}: paired FASTQ validation passed for {r1_name} + {r2_name} ({checked} pair(s))", log_path)
            return True
        log(
            f"{step_name}: paired FASTQ validation failed for {r1_name} + {r2_name}; "
            f"checked {checked} pair(s); {msg}. Skipping this step.",
            log_path,
        )
        return False

    if "clumpify-rnaviralspades" in selected:
        if not validate_current_pair_before_step("clumpify", assembly_r1, assembly_r2):
            log(f"Skipping clumpify-rnaviralspades in {sample_dir.name} because input reads are not synchronized", log_path)
        else:
            ok = run_command(
                [
                    "clumpify.sh",
                    f"in1={assembly_r1}",
                    f"in2={assembly_r2}",
                    "out1=deduped_R1.fastq.gz",
                    "out2=deduped_R2.fastq.gz",
                    "dedupe",
                ],
                log_path,
                cwd=sample_dir,
            )
        if not ok:
            log(f"Clumpify failed in {sample_dir.name}; skipping clumpify-rnaviralspades", log_path)

    bbnorm_modes = [
        (
            "bbnorm100-rnaviralspades",
            100,
            "R1bbnorm.fastq.gz",
            "R2bbnorm.fastq.gz",
            "work1",
        ),
        (
            "bbnorm500-rnaviralspades",
            500,
            "R1bbnorm2.fastq.gz",
            "R2bbnorm2.fastq.gz",
            "work2",
        ),
    ]
    selected_bbnorm_modes = [mode for mode in bbnorm_modes if mode[0] in selected]
    if selected_bbnorm_modes:
        if not validate_current_pair_before_step("bbnorm", assembly_r1, assembly_r2):
            log(
                f"Skipping BBNorm assembly modes in {sample_dir.name} because "
                "the input reads are not synchronized",
                log_path,
            )
        else:
            with tempfile.TemporaryDirectory(
                prefix="rotax_bbnorm_interleaved_",
                dir=sample_dir,
            ) as temp_name:
                bbnorm_temp_dir = Path(temp_name)
                interleaved_input = bbnorm_temp_dir / "input.interleaved.fastq.gz"
                input_ok = run_command(
                    [
                        "reformat.sh",
                        f"in1={assembly_r1}",
                        f"in2={assembly_r2}",
                        f"out={interleaved_input}",
                        "overwrite=t",
                    ],
                    log_path,
                    cwd=sample_dir,
                )
                if input_ok:
                    log(
                        f"Prepared one interleaved BBNorm input stream for "
                        f"{sample_dir.name}",
                        log_path,
                    )

                for (
                    mode_name,
                    target,
                    out1_name,
                    out2_name,
                    work_dir,
                ) in selected_bbnorm_modes:
                    ok = input_ok and run_bbnorm_interleaved(
                        sample_dir,
                        interleaved_input,
                        target,
                        out1_name,
                        out2_name,
                        bbnorm_threads,
                        bbnorm_temp_dir,
                        log_path,
                    )
                    if ok:
                        run_command(
                            [
                                "rnaviralspades.py",
                                "-t", str(threads),
                                "-1", out1_name,
                                "-2", out2_name,
                                "-o", work_dir,
                            ],
                            log_path,
                            cwd=sample_dir,
                        )
                    else:
                        log(
                            f"BBNorm target={target} failed or output was "
                            f"missing/invalid in {sample_dir.name}; skipping "
                            f"{mode_name}. Consider changing the installed "
                            "BBTools version if single-stream BBNorm also fails.",
                            log_path,
                        )

    if "rnaviralspades" in selected:
        if validate_current_pair_before_step("rnaviralspades", assembly_r1, assembly_r2):
            run_command(
                [
                    "rnaviralspades.py",
                    "-t", str(threads),
                    "-1", assembly_r1,
                    "-2", assembly_r2,
                    "-o", "work3",
                ],
                log_path,
                cwd=sample_dir,
            )
        else:
            log(f"Skipping rnaviralspades in {sample_dir.name} because input reads are not synchronized", log_path)

    if "spades-isolate" in selected:
        if validate_current_pair_before_step("spades-isolate", assembly_r1, assembly_r2):
            run_command(
                [
                    "spades.py",
                    "-t", str(threads),
                    "--isolate",
                    "--cov-cutoff", "auto",
                    "-1", assembly_r1,
                    "-2", assembly_r2,
                    "-o", "work6",
                ],
                log_path,
                cwd=sample_dir,
            )
        else:
            log(f"Skipping spades-isolate in {sample_dir.name} because input reads are not synchronized", log_path)

    if "clumpify-rnaviralspades" in selected:
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
            log(f"Deduped output missing/empty in {sample_dir.name}; skipping clumpify-rnaviralspades", log_path)

    collect_and_filter_contigs(sample_dir, log_path, assemblies)
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


def coinfection_display(constellation: str, extra_information: str) -> tuple[str, str]:
    """Return the G/P-only Labware display genotype and combined G/P calls.

    Keep the selected genotype first; append distinct alternatives numerically.
    Two resolved genotypes in any segment constitute the existing mixed signal.
    Unresolved X/? calls alone do not establish a coinfection.
    """
    calls: dict[str, list[str]] = {gene: [] for gene in GENE_ORDER}
    primary = dict(zip(GENE_ORDER, constellation.split("-"))) if constellation else {}
    for token in constellation.split("-"):
        gene, number = extract_gene_and_number(token)
        if gene in calls and number.isdigit() and number not in calls[gene]:
            calls[gene].append(number)
    alternatives: dict[str, set[str]] = {gene: set() for gene in GENE_ORDER}
    for token in extra_information.split(";"):
        gene, number = extract_gene_and_number(token)
        if gene in alternatives and number.isdigit():
            alternatives[gene].add(number)
    for gene in GENE_ORDER:
        calls[gene].extend(number for number in sorted(alternatives[gene], key=int)
                           if number not in calls[gene])

    def joined(gene: str) -> str:
        return "".join(format_gene_genotype(gene, number) for number in calls[gene])

    g_p = "-".join(joined(gene) or primary.get(gene, f"{gene}X") for gene in ("G", "P"))
    if not constellation and not any(calls.values()):
        g_p = ""
    mixed_genes = [gene for gene in GENE_ORDER if len(calls[gene]) > 1]
    if not mixed_genes:
        return g_p, g_p
    return f"Koinfeksjon({g_p})", g_p


def rotateq_closest_identities(
    sample_dir: Path, blast_rows: list[dict[str, str]],
) -> dict[str, str]:
    """Highest accepted RotaTeq identity per segment among selected contigs.

    Compare same-query, same-segment hits by bitscore then pident, regardless
    of acceptance, before checking the winning RotaTeq hit's acceptance.
    A reported exact tie with a non-RotaTeq reference is left blank for that
    contig. The input contains the references retained by the genotyping BLAST
    (currently max_target_seqs=1), not an exhaustive database comparison.
    """
    selected = {
        (str(row["qseqid"]), str(row["segment"]))
        for row in read_selected_candidate_rows(sample_dir).values()
    }
    best: dict[tuple[str, str], tuple[tuple[float, float], list[dict[str, str]]]] = {}
    for row in blast_rows:
        segment = segment_from_reference_id(row.get("sseqid", ""))
        if not segment:
            # Some supplied vaccine headers identify only the segment (G|, R|,
            # etc.). Report identity without inventing a numeric genotype.
            prefix = row.get("sseqid", "").split("|", 1)[0]
            segment = GENOTYPE_PREFIX_TO_SEGMENT.get(prefix, "")
        key = (row.get("qseqid", ""), segment)
        if key not in selected:
            continue
        try:
            rank = (float(row["bitscore"]), float(row["pident"]))
        except (KeyError, ValueError):
            continue
        if not all(math.isfinite(value) for value in rank) or not 0 <= rank[1] <= 100:
            continue
        if key not in best or rank > best[key][0]:
            best[key] = (rank, [row])
        elif rank == best[key][0]:
            best[key][1].append(row)

    identities: dict[str, float] = {}
    for (_, segment), (rank, winners) in best.items():
        if any("rotateq" not in row.get("sseqid", "").lower() for row in winners):
            continue
        if not any(row.get("Evaluation", "") not in {"", "not accepted"} for row in winners):
            continue
        identities[segment] = max(identities.get(segment, 0.0), rank[1])
    return {f"RotaTeq {segment}": str(identities[segment]) if segment in identities else ""
            for segment in SEGMENT_ORDER}


def add_vaccine_match_columns(sample_dir: Path, log_path: Path) -> None:
    """Write display calls and vaccine identities, preserving the raw constellation.

    Rotarix retains the mean identity of accepted G1/P8 RIX4414 hits.
    RotaTeq uses the highest qualifying closest-reference identity per segment.
    """
    genotyping_file = sample_dir / "rotavirus_genotype_constellation_intermediate.csv"
    blast_file = sample_dir / "rotavirus_blast_evaluated.csv"
    output_file = sample_dir / "rotax_sample_genotyping_report.csv"
    with genotyping_file.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        rows = list(csv.DictReader(handle))
    with blast_file.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        blast_rows = list(csv.DictReader(handle))

    def average_pident(evaluation_value: str, header_marker: str) -> str:
        values: list[float] = []
        for row in blast_rows:
            if row.get("Evaluation") != evaluation_value:
                continue
            if header_marker.lower() not in row.get("sseqid", "").lower():
                continue
            try:
                values.append(float(row.get("pident", "")))
            except ValueError:
                pass
        return str(sum(values) / len(values)) if values else ""

    vaccine_values = {
        "Rotarix VP7": average_pident("G1", "RIX4414"),
        "Rotarix VP4": average_pident("P8", "RIX4414"),
        **rotateq_closest_identities(sample_dir, blast_rows),
    }
    with output_file.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=SAMPLE_GENOTYPING_FIELDS)
        writer.writeheader()
        for row in rows:
            constellation = row.get("Genotype", "")
            display, _ = coinfection_display(constellation, row.get("Extra Information", ""))
            row["LabwareGenotype"] = display
            row.update(vaccine_values)
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


# ---------------------------------------------------------------------------
# GenBank full-segment evaluation
# ---------------------------------------------------------------------------

def normalize_dna(seq: str) -> str:
    return re.sub(r"\s+", "", seq).upper().replace("U", "T").replace("-", "")


def parse_genotype_token(text: str) -> str:
    """Extract a genotype token such as G3, P8, T28 or GX from a header/id."""
    token = text.lstrip(">").split()[0]
    token = token.split("|", 1)[0]
    token = token.split("_", 1)[0]
    token = token.replace("[", "").replace("]", "")
    token = token.split("PARTIAL", 1)[0].strip("_")
    match = re.match(r"^([GPIRCMANTEHgpircmanteh])(\d+|X|\?)", token)
    if not match:
        return ""
    return match.group(1).upper() + match.group(2)


def segment_from_genotype(genotype: str) -> str:
    if not genotype:
        return ""
    return GENOTYPE_PREFIX_TO_SEGMENT.get(genotype[0].upper(), "")


def segment_from_reference_id(sseqid: str) -> str:
    return segment_from_genotype(parse_genotype_token(sseqid))


def read_selected_candidate_rows(sample_dir: Path) -> dict[str, dict[str, float | int | str]]:
    """Return the same selected candidates used for selected_contigs/selected_orfs."""
    blast_csv = sample_dir / "rotavirus_blast_evaluated.csv"
    full_candidates: dict[str, dict[str, float | int | str]] = {}
    partial_candidates: dict[str, dict[str, float | int | str]] = {}

    if not blast_csv.exists():
        return {}

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
                pident = float(row["pident"])
                bitscore = float(row.get("bitscore", "0") or 0)
            except ValueError:
                continue

            qseqid = row["qseqid"]
            cov = get_coverage(qseqid)
            info: dict[str, float | int | str] = {
                "genotype": genotype,
                "segment": segment_from_genotype(genotype),
                "qseqid": qseqid,
                "coverage": cov,
                "length": length,
                "slen": slen,
                "qstart": qstart,
                "qend": qend,
                "sstart": row.get("sstart", ""),
                "send": row.get("send", ""),
                "sseqid": row.get("sseqid", ""),
                "sstrand": row.get("sstrand", "+"),
                "pident": pident,
                "bitscore": bitscore,
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
            selected[genotype]["orf_candidate_type"] = "full"
        else:
            selected[genotype] = dict(partial_candidates[genotype])
            selected[genotype]["orf_candidate_type"] = "partial"
    return selected


def terminal_homopolymer_run(seq: str, from_left: bool = True) -> tuple[str, int]:
    if not seq:
        return "", 0
    s = normalize_dna(seq)
    if not s:
        return "", 0
    if from_left:
        base = s[0]
        i = 0
        while i < len(s) and s[i] == base:
            i += 1
        return base, i
    base = s[-1]
    i = len(s) - 1
    while i >= 0 and s[i] == base:
        i -= 1
    return base, len(s) - 1 - i


def terminal_non_acgt_trim(seq: str) -> tuple[int, int]:
    s = normalize_dna(seq)
    left = 0
    while left < len(s) and s[left] not in "ACGT":
        left += 1
    right = 0
    while right < len(s) - left and s[len(s) - 1 - right] not in "ACGT":
        right += 1
    return left, right


def terminal_trim_counts(seq: str, homopolymer_threshold: int) -> tuple[int, int, str]:
    s = normalize_dna(seq)
    left_non, right_non = terminal_non_acgt_trim(s)
    left = left_non
    right = right_non
    notes: list[str] = []

    core = s[left:len(s) - right if right else len(s)]
    if core:
        lb, lr = terminal_homopolymer_run(core, from_left=True)
        rb, rr = terminal_homopolymer_run(core, from_left=False)
        if lr >= homopolymer_threshold:
            left += lr
            notes.append(f"left_terminal_homopolymer_{lb}_{lr}")
        if rr >= homopolymer_threshold and len(core) - (lr if lr >= homopolymer_threshold else 0) > rr:
            right += rr
            notes.append(f"right_terminal_homopolymer_{rb}_{rr}")

    if left_non:
        notes.append(f"left_non_acgt_trim_{left_non}")
    if right_non:
        notes.append(f"right_non_acgt_trim_{right_non}")
    return left, right, ";".join(notes)


def count_ambiguous_bases(seq: str) -> int:
    s = normalize_dna(seq)
    return sum(1 for ch in s if ch not in "ACGT")


# Conserved RVA terminal motifs in DNA/IUPAC notation.
# W = A/T and R = A/G. These match the full-segment DB builder logic.
TERMINAL_5_MOTIF = "GGCWWTWWWW"
TERMINAL_3_MOTIF_DEFAULT = "WTRTGRCC"
TERMINAL_3_MOTIF_NSP1_SPECIAL = "CTGTGAACC"

IUPAC_BASES = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def iupac_set(base: str) -> set[str]:
    return IUPAC_BASES.get(base.upper().replace("U", "T"), {base.upper().replace("U", "T")})


def motif_mismatch_count(observed: str, expected: str) -> int:
    """Compare observed DNA sequence to an IUPAC/degenerated motif.

    A position passes when the observed base's possible set intersects the
    motif position's allowed set. Missing observed positions count as
    mismatches. Ambiguous observed bases can pass if compatible with the motif.
    """
    observed = normalize_dna(observed)
    expected = expected.upper().replace("U", "T")
    mismatches = 0
    for i, expected_base in enumerate(expected):
        if i >= len(observed):
            mismatches += 1
            continue
        if not (iupac_set(observed[i]) & iupac_set(expected_base)):
            mismatches += 1
    return mismatches


def evaluate_one_terminal_motif(seq: str, expected: str, from_left: bool, allowed_mismatches: int) -> dict[str, str | int | bool]:
    s = normalize_dna(seq)
    expected = expected.upper().replace("U", "T")
    if from_left:
        observed = s[:len(expected)]
    else:
        observed = s[-len(expected):] if len(s) >= len(expected) else s
    mismatches = motif_mismatch_count(observed, expected)
    passed = mismatches <= allowed_mismatches
    return {
        "expected": expected,
        "observed": observed,
        "mismatches": mismatches,
        "status": "PASS" if passed else "FAIL",
        "pass": passed,
    }


def evaluate_terminal_motifs(seq: str, segment: str, allowed_mismatches: int) -> dict[str, str | int | bool]:
    """Evaluate conserved RVA terminal motifs on the cleaned, positive-sense sequence.

    NSP1/segment 5 is accepted if either the normal 8-mer 3-prime motif
    WTRTGRCC passes or the special 9-mer CTGTGAACC passes. This mirrors the
    full-segment reference DB filtering update after observing many valid NSP1
    GenBank records ending in AATGTGACC, where the last 8 bases match WTRTGRCC.
    """
    s = normalize_dna(seq)
    eval5 = evaluate_one_terminal_motif(s, TERMINAL_5_MOTIF, True, allowed_mismatches)
    default3 = evaluate_one_terminal_motif(s, TERMINAL_3_MOTIF_DEFAULT, False, allowed_mismatches)
    special3 = evaluate_one_terminal_motif(s, TERMINAL_3_MOTIF_NSP1_SPECIAL, False, allowed_mismatches)

    if segment == "NSP1":
        if bool(default3["pass"]):
            chosen3 = default3
            terminal3_mode = "default_8mer"
        elif bool(special3["pass"]):
            chosen3 = special3
            terminal3_mode = "nsp1_special_9mer"
        else:
            # Report the closer failing motif as the main observed/expected pair.
            chosen3 = default3 if int(default3["mismatches"]) <= int(special3["mismatches"]) else special3
            terminal3_mode = "fail_both_default_8mer_and_nsp1_special_9mer"
        pass3 = bool(default3["pass"] or special3["pass"])
    else:
        chosen3 = default3
        terminal3_mode = "default_8mer" if bool(default3["pass"]) else "fail_default_8mer"
        pass3 = bool(default3["pass"])

    pass5 = bool(eval5["pass"])
    overall = pass5 and pass3
    if overall:
        classification = f"terminal_motifs_pass:{terminal3_mode}"
    else:
        failed = []
        if not pass5:
            failed.append("5prime_motif_fail")
        if not pass3:
            failed.append("3prime_motif_fail")
        classification = ";".join(failed) + f":{terminal3_mode}"

    return {
        "terminal_motif_classification": classification,
        "terminal_motifs_pass": overall,
        "terminal_5_expected": eval5["expected"],
        "terminal_5_observed": eval5["observed"],
        "terminal_5_mismatches": eval5["mismatches"],
        "terminal_5_status": eval5["status"],
        "terminal_3_expected": chosen3["expected"],
        "terminal_3_observed": chosen3["observed"],
        "terminal_3_mismatches": chosen3["mismatches"],
        "terminal_3_status": "PASS" if pass3 else "FAIL",
        "terminal_3_mode": terminal3_mode,
        "terminal_3_default_expected": default3["expected"],
        "terminal_3_default_observed": default3["observed"],
        "terminal_3_default_mismatches": default3["mismatches"],
        "terminal_3_default_status": default3["status"],
        "terminal_3_nsp1_special_expected": special3["expected"] if segment == "NSP1" else "",
        "terminal_3_nsp1_special_observed": special3["observed"] if segment == "NSP1" else "",
        "terminal_3_nsp1_special_mismatches": special3["mismatches"] if segment == "NSP1" else "",
        "terminal_3_nsp1_special_status": special3["status"] if segment == "NSP1" else "",
    }


def validate_orf(seq: str) -> dict[str, str | int | bool]:
    s = normalize_dna(seq)
    starts_atg = len(s) >= 3 and s.startswith("ATG")
    ends_stop = len(s) >= 3 and s[-3:] in STOP_CODONS
    length_multiple_3 = bool(s) and len(s) % 3 == 0
    internal_stops = 0
    if len(s) >= 6:
        stop_check_end = len(s) - 3 if ends_stop else len(s)
        for i in range(3, stop_check_end, 3):
            if s[i:i + 3] in STOP_CODONS:
                internal_stops += 1
    valid = starts_atg and ends_stop and length_multiple_3 and internal_stops == 0
    if valid:
        status = "valid_orf"
    else:
        problems: list[str] = []
        if not starts_atg:
            problems.append("missing_ATG")
        if not ends_stop:
            problems.append("missing_stop")
        if not length_multiple_3:
            problems.append("length_not_multiple_of_3")
        if internal_stops:
            problems.append(f"internal_stop_count_{internal_stops}")
        status = ";".join(problems) if problems else "invalid_orf"
    return {
        "valid_orf": valid,
        "orf_validation_status": status,
        "starts_atg": starts_atg,
        "ends_stop": ends_stop,
        "length_multiple_3": length_multiple_3,
        "internal_stop_codons": internal_stops,
        "length": len(s),
    }


def oriented_sequence_and_cds(candidate: dict[str, float | int | str], seq: str) -> tuple[str, int, int, str]:
    """Orient the full contig according to the ORF BLAST hit and return CDS coords."""
    s = normalize_dna(seq)
    qstart = int(candidate["qstart"])
    qend = int(candidate["qend"])
    left = min(qstart, qend)
    right = max(qstart, qend)
    strand = str(candidate.get("sstrand", "+")).strip().lower()
    if strand == "minus":
        oriented = reverse_complement(s)
        cds_start = len(s) - right + 1
        cds_end = len(s) - left + 1
        orientation = "reverse_complemented_from_orf_blast"
    else:
        oriented = s
        cds_start = left
        cds_end = right
        orientation = "as_assembled_orf_blast_plus"
    return oriented, cds_start, cds_end, orientation


def safe_query_id(sample: str, genotype: str, qseqid: str) -> str:
    raw = f"{sample}|{genotype}|{qseqid}"
    return re.sub(r"[^A-Za-z0-9_.|:-]+", "_", raw)[:230]


def write_empty_genbank_reports(sample_dir: Path, args: argparse.Namespace, log_path: Path, reason: str = "no_records") -> None:
    report = sample_dir / "rotavirus_genbank_full_segment_report.csv"
    coords = sample_dir / "rotavirus_genbank_feature_coordinates.csv"
    clean_fasta = sample_dir / "rotavirus_genbank_clean_full_segments.fasta"
    review_fasta = sample_dir / "rotavirus_genbank_partial_or_review_segments.fasta"
    write_fasta([], clean_fasta)
    write_fasta([], review_fasta)
    report_fields = genbank_report_fields()
    with report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=report_fields)
        writer.writeheader()
    with coords.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=genbank_coordinate_fields())
        writer.writeheader()
    log(f"Generated empty GenBank evaluation files for {sample_dir.name}: {reason}", log_path)


def genbank_report_fields() -> list[str]:
    return [
        "Sample", "Segment", "Genotype", "Status", "SubmitToGenBank", "Reason", "Complete11SegmentSample",
        "OriginalContigID", "OriginalContigLength", "CleanLength", "Orientation", "Coverage",
        "BestFullReference", "BestFullReferenceGenotype", "BestFullReferenceSegment", "BestFullReferenceLength",
        "FullReferencePident", "FullReferenceBitscore", "ReferenceAlignedStart", "ReferenceAlignedEnd",
        "ReferenceCoveragePct", "ContigAlignedStart", "ContigAlignedEnd", "ContigCoveragePct",
        "MissingReference5primeBases", "MissingReference3primeBases", "ExtraContig5primeBases", "ExtraContig3primeBases",
        "AdjustedExtraContig5primeBases", "AdjustedExtraContig3primeBases", "TrimmedLeftBases", "TrimmedRightBases",
        "TerminalTrimNotes", "TerminalMotifClassification", "TerminalMotifsPass",
        "Terminal5Expected", "Terminal5Observed", "Terminal5Mismatches", "Terminal5Status",
        "Terminal3Expected", "Terminal3Observed", "Terminal3Mismatches", "Terminal3Status",
        "Terminal3Mode", "Terminal3DefaultExpected", "Terminal3DefaultObserved",
        "Terminal3DefaultMismatches", "Terminal3DefaultStatus",
        "Terminal3NSP1SpecialExpected", "Terminal3NSP1SpecialObserved",
        "Terminal3NSP1SpecialMismatches", "Terminal3NSP1SpecialStatus",
        "CDSStart", "CDSEnd", "CDSLength", "ExpectedORFLength", "ORFLengthDifference",
        "ORFLengthTolerance", "ORFLengthWithinTolerance",
        "ORFValidationStatus", "StartsWithATG", "EndsWithStopCodon", "LengthMultipleOf3", "InternalStopCodons",
        "HasAmbiguity", "AmbiguousBases", "LeftTerminalHomopolymerBase", "LeftTerminalHomopolymerLength",
        "RightTerminalHomopolymerBase", "RightTerminalHomopolymerLength", "ORFBlastSubject", "ORFBlastPident",
        "ORFCandidateType", "FullSegmentDB", "FullSegmentSameSegmentHitCount", "FullSegmentAllHitCount",
        "AlignmentInputFasta", "AlignmentFasta", "AlignmentReferenceCount",
    ]


def genbank_coordinate_fields() -> list[str]:
    return [
        "Sample", "Segment", "Genotype", "Status", "SubmitToGenBank", "SequenceID",
        "CleanLength", "Feature", "Start", "End", "Gene", "Product", "Note",
    ]


def run_full_segment_blast(sample_dir: Path, full_segment_db: str, query_fasta: Path, log_path: Path) -> Path | None:
    blast_output = sample_dir / "rotavirus_full_segment_blast_raw.tsv"
    if not full_segment_db:
        return None
    if not query_fasta.exists() or not parse_fasta(query_fasta):
        return None
    ok = run_command(
        [
            "blastn",
            "-query", str(query_fasta.name),
            "-db", full_segment_db,
            "-out", str(blast_output.name),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand slen qlen",
            "-max_target_seqs", "50",
            "-max_hsps", "1",
            "-strand", "both",
        ],
        log_path,
        cwd=sample_dir,
    )
    return blast_output if ok and blast_output.exists() else None


def parse_full_segment_blast(blast_tsv: Path) -> dict[str, list[dict[str, str | int | float]]]:
    hits: dict[str, list[dict[str, str | int | float]]] = defaultdict(list)
    if not blast_tsv.exists():
        return hits
    fields = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "sstrand", "slen", "qlen",
    ]
    with blast_tsv.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if len(parts) < len(fields):
                continue
            row = dict(zip(fields, parts))
            try:
                hit = {
                    "qseqid": row["qseqid"],
                    "sseqid": row["sseqid"],
                    "pident": float(row["pident"]),
                    "length": int(row["length"]),
                    "qstart": int(row["qstart"]),
                    "qend": int(row["qend"]),
                    "sstart": int(row["sstart"]),
                    "send": int(row["send"]),
                    "bitscore": float(row["bitscore"]),
                    "sstrand": row["sstrand"],
                    "slen": int(row["slen"]),
                    "qlen": int(row["qlen"]),
                }
            except ValueError:
                continue
            hits[str(hit["qseqid"])].append(hit)
    return hits


def full_segment_hit_metrics(hit: dict[str, str | int | float]) -> dict[str, int | float | str]:
    qstart = int(hit["qstart"])
    qend = int(hit["qend"])
    sstart = int(hit["sstart"])
    send = int(hit["send"])
    slen = int(hit["slen"])
    qlen = int(hit["qlen"])
    q_left = min(qstart, qend)
    q_right = max(qstart, qend)
    s_left = min(sstart, send)
    s_right = max(sstart, send)
    ref_span = s_right - s_left + 1
    query_span = q_right - q_left + 1
    return {
        "reference_aligned_start": s_left,
        "reference_aligned_end": s_right,
        "reference_coverage_pct": 100.0 * ref_span / slen if slen else 0.0,
        "contig_aligned_start": q_left,
        "contig_aligned_end": q_right,
        "contig_coverage_pct": 100.0 * query_span / qlen if qlen else 0.0,
        "missing_ref_5": s_left - 1,
        "missing_ref_3": slen - s_right,
        "extra_query_5": q_left - 1,
        "extra_query_3": qlen - q_right,
    }


def choose_best_full_segment_hit(
    hits: list[dict[str, str | int | float]],
    expected_segment: str,
) -> tuple[dict[str, str | int | float] | None, int]:
    same_segment = [h for h in hits if segment_from_reference_id(str(h["sseqid"])) == expected_segment]
    if not same_segment:
        return None, 0

    def score(hit: dict[str, str | int | float]) -> tuple[float, float]:
        # v10 ranking rule: choose the closest same-segment full reference by
        # BLAST evidence only. Do not let reference coverage or terminal/extra
        # end metrics override the BLAST best hit, because that can select a
        # biologically less similar reference with nicer terminal alignment.
        return (
            float(hit["bitscore"]),
            float(hit["pident"]),
        )

    best = sorted(same_segment, key=score, reverse=True)[0]
    return best, len(same_segment)


def ranked_same_segment_full_hits(
    hits: list[dict[str, str | int | float]],
    expected_segment: str,
    limit: int = 2,
) -> list[dict[str, str | int | float]]:
    """Return the best same-segment full-reference hits ranked by bitscore, then percent identity."""
    same_segment = [h for h in hits if segment_from_reference_id(str(h["sseqid"])) == expected_segment]

    def score(hit: dict[str, str | int | float]) -> tuple[float, float]:
        # v10 ranking rule: choose the closest same-segment full reference by
        # BLAST evidence only. Do not let reference coverage or terminal/extra
        # end metrics override the BLAST best hit, because that can select a
        # biologically less similar reference with nicer terminal alignment.
        return (
            float(hit["bitscore"]),
            float(hit["pident"]),
        )

    ranked = sorted(same_segment, key=score, reverse=True)
    out: list[dict[str, str | int | float]] = []
    seen: set[str] = set()
    for hit in ranked:
        sid = str(hit["sseqid"])
        if sid in seen:
            continue
        out.append(hit)
        seen.add(sid)
        if len(out) >= limit:
            break
    return out


def safe_file_component(text: str, max_len: int = 120) -> str:
    text = str(text or "").strip()
    text = re.sub(r"[^A-Za-z0-9_.:-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return (text or "sequence")[:max_len]


def parse_fasta_text(text: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(parts)))
            header = line
            parts = []
        else:
            parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    return records


FASTA_SOURCE_CACHE: dict[str, dict[str, tuple[str, str]]] = {}


def likely_fasta_sources_for_blast_db(db: str) -> list[Path]:
    """Return likely source FASTA files for a BLAST DB prefix.

    The alignment QC needs the actual reference sequences. ``blastdbcmd -entry``
    only works reliably when a database was created with parsed sequence IDs,
    and older RotaX/BLAST databases may not have been. Therefore we first try
    to recover the reference sequence from the source FASTA that normally sits
    next to the BLAST DB prefix.
    """
    if not db:
        return []

    p = Path(db)
    candidates: list[Path] = []

    # If the user accidentally supplied a FASTA file, use it directly.
    if p.is_file() and p.suffix.lower() in {".fa", ".fasta", ".fna", ".fas", ".txt"}:
        candidates.append(p)

    # Most RotaX DB prefixes have the source FASTA as <prefix>.fasta.
    for suffix in [".fasta", ".fa", ".fna", ".fas", ".txt"]:
        candidates.append(Path(str(p) + suffix))

    # Also check common source FASTA names in the DB directory.
    parent = p.parent if str(p.parent) else Path(".")
    stem = p.name
    common_names = [
        f"{stem}.fasta",
        f"{stem}.fa",
        "References_all_segments_plus_vaccine.fasta",
        "References_all_segments_plus_vaccine.fa",
        "rotavirus_full_segments.fasta",
        "rotavirus_full_segments.fa",
    ]
    for name in common_names:
        candidates.append(parent / name)

    out: list[Path] = []
    seen: set[Path] = set()
    for cand in candidates:
        try:
            resolved = cand.resolve()
        except OSError:
            resolved = cand
        if resolved in seen:
            continue
        seen.add(resolved)
        if cand.exists() and cand.is_file():
            out.append(cand)
    return out


def load_fasta_source_index(path: Path) -> dict[str, tuple[str, str]]:
    cache_key = str(path.resolve())
    if cache_key in FASTA_SOURCE_CACHE:
        return FASTA_SOURCE_CACHE[cache_key]

    index: dict[str, tuple[str, str]] = {}
    for header, seq in parse_fasta(path):
        key = fasta_key(header)
        clean_seq = normalize_dna(seq)
        if key and key not in index:
            index[key] = (header, clean_seq)

        # Some BLAST DBs return only the accession-like part of an ID.
        # Add a few harmless aliases without overwriting exact keys.
        token = key.split("|", 1)[0] if key else ""
        if token and token not in index:
            index[token] = (header, clean_seq)
        for part in key.split("|") if key else []:
            if part and part not in index:
                index[part] = (header, clean_seq)

    FASTA_SOURCE_CACHE[cache_key] = index
    return index


def fetch_sequence_from_fasta_sources(db: str, entry: str) -> tuple[str, str] | None:
    """Fetch a reference sequence from source FASTA files near a BLAST DB prefix."""
    if not db or not entry:
        return None

    entry = str(entry).strip()
    entry_aliases = [entry]
    if entry.startswith("lcl|"):
        entry_aliases.append(entry.split("lcl|", 1)[1])
    if " " in entry:
        entry_aliases.append(entry.split()[0])

    for fasta_path in likely_fasta_sources_for_blast_db(db):
        index = load_fasta_source_index(fasta_path)
        for alias in entry_aliases:
            if alias in index:
                return index[alias]
    return None


def fetch_blastdb_sequence(db: str, entry: str, log_path: Path) -> tuple[str, str] | None:
    """Fetch a sequence for alignment QC from source FASTA or BLAST DB.

    Prefer the source FASTA next to the BLAST DB prefix. This avoids noisy
    ``blastdbcmd -entry`` failures when the database was built without
    ``makeblastdb -parse_seqids``. Fall back to blastdbcmd only when FASTA
    lookup fails.
    """
    if not db or not entry:
        return None

    fetched = fetch_sequence_from_fasta_sources(db, entry)
    if fetched is not None:
        return fetched

    if shutil.which("blastdbcmd") is None:
        log(
            f"WARNING: could not fetch entry {entry} for alignment: no source FASTA found near {db} "
            "and blastdbcmd is not in PATH",
            log_path,
        )
        return None

    command = ["blastdbcmd", "-db", db, "-entry", entry, "-outfmt", "%f"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        log(
            f"WARNING: could not fetch entry {entry} from source FASTA or BLAST DB {db}: "
            f"{result.stderr.strip()}",
            log_path,
        )
        return None
    records = parse_fasta_text(result.stdout)
    if not records:
        log(f"WARNING: blastdbcmd returned no FASTA record for entry {entry} from {db}", log_path)
        return None
    return records[0]


def orient_reference_like_hit(seq: str, hit: dict[str, str | int | float]) -> str:
    """Orient a fetched reference sequence the same way it aligned to the candidate."""
    s = normalize_dna(seq)
    strand = str(hit.get("sstrand", "+")).strip().lower()
    if strand == "minus":
        return reverse_complement(s)
    return s


def write_candidate_reference_alignment(
    sample_dir: Path,
    args: argparse.Namespace,
    log_path: Path,
    qid: str,
    candidate: dict[str, str | int | float],
    segment: str,
    genotype: str,
    status: str,
    cleaned_sequence: str,
    full_hits: list[dict[str, str | int | float]],
) -> tuple[str, str, int]:
    """Write and, when MAFFT is available, align candidate + refs for manual inspection."""
    if not bool(getattr(args, "genbank_make_alignments", True)):
        return "", "", 0

    align_dir = sample_dir / "alignments"
    align_dir.mkdir(exist_ok=True)

    file_stem = safe_file_component(f"{segment}_{genotype}_{qid}", 180)
    input_fasta = align_dir / f"{file_stem}.fasta"
    aligned_fasta = align_dir / f"{file_stem}.aln.fasta"

    records: list[tuple[str, str]] = []
    records.append((
        f">CANDIDATE|{sample_dir.name}|{segment}|{genotype}|status={status}|source={safe_file_component(candidate.get('qseqid', ''), 80)}",
        normalize_dna(cleaned_sequence),
    ))

    full_db = str(getattr(args, "full_segment_db", "") or "")
    for idx, hit in enumerate(full_hits, start=1):
        entry = str(hit["sseqid"])
        fetched = fetch_blastdb_sequence(full_db, entry, log_path)
        if fetched is None:
            continue
        _header, ref_seq = fetched
        ref_seq = orient_reference_like_hit(ref_seq, hit)
        metrics = full_segment_hit_metrics(hit)
        records.append((
            f">FULLREF{idx}|{entry}|pident={float(hit['pident']):.3f}|refcov={float(metrics['reference_coverage_pct']):.2f}|sstrand={hit.get('sstrand', '')}",
            ref_seq,
        ))

    orf_db = str(getattr(args, "blast_db", "") or "")
    orf_entry = str(candidate.get("sseqid", "") or "")
    fetched_orf = fetch_blastdb_sequence(orf_db, orf_entry, log_path) if orf_entry else None
    if fetched_orf is not None:
        _orf_header, orf_seq = fetched_orf
        # The candidate full contig has already been oriented according to the
        # ORF BLAST hit. The ORF reference sequence from the BLAST DB is kept
        # in its database/positive-sense orientation.
        records.append((
            f">ORFREF_USED_FOR_GENOTYPE|{orf_entry}|orf_blast_pident={candidate.get('pident', '')}",
            normalize_dna(orf_seq),
        ))

    write_fasta(records, input_fasta)

    aligner = str(getattr(args, "genbank_alignment_cmd", "mafft") or "mafft")
    if len(records) >= 2 and aligner.lower() != "none" and shutil.which(aligner):
        command = [aligner, "--auto", str(input_fasta.name)]
        log("Running alignment: " + " ".join(command) + f" > {aligned_fasta.name}", log_path)
        with aligned_fasta.open("w", encoding="utf-8") as out_handle:
            result = subprocess.run(command, cwd=align_dir, stdout=out_handle)
        if result.returncode != 0:
            log(f"WARNING: alignment failed for {input_fasta}; keeping unaligned FASTA only", log_path)
            if aligned_fasta.exists():
                aligned_fasta.unlink()
            return str(input_fasta), "", len(records) - 1
    else:
        if len(records) < 2:
            log(f"WARNING: only candidate sequence available for {input_fasta}; no alignment created", log_path)
        elif aligner.lower() != "none":
            log(f"WARNING: alignment command {aligner} not found; wrote unaligned FASTA only: {input_fasta}", log_path)
        return str(input_fasta), "", len(records) - 1

    return str(input_fasta), str(aligned_fasta), len(records) - 1


def evaluate_genbank_full_segments(sample_dir: Path, args: argparse.Namespace, log_path: Path) -> None:
    report_path = sample_dir / "rotavirus_genbank_full_segment_report.csv"
    coords_path = sample_dir / "rotavirus_genbank_feature_coordinates.csv"
    clean_fasta = sample_dir / "rotavirus_genbank_clean_full_segments.fasta"
    review_fasta = sample_dir / "rotavirus_genbank_partial_or_review_segments.fasta"
    query_fasta = sample_dir / "rotavirus_genbank_full_segment_query.fasta"

    candidates = read_selected_candidate_rows(sample_dir)
    fasta_records = {fasta_key(h): (h, s) for h, s in parse_fasta(sample_dir / "contigs500_COV.fasta")}

    query_records: list[tuple[str, str]] = []
    query_map: dict[str, dict[str, str | int | float]] = {}
    for genotype in sorted(candidates):
        candidate = candidates[genotype]
        qseqid = str(candidate["qseqid"])
        if qseqid not in fasta_records:
            continue
        original_header, raw_seq = fasta_records[qseqid]
        oriented, cds_start, cds_end, orientation = oriented_sequence_and_cds(candidate, raw_seq)
        qid = safe_query_id(sample_dir.name, genotype, qseqid)
        query_records.append((f">{qid}", oriented))
        item = dict(candidate)
        item.update({
            "query_id": qid,
            "original_header": original_header,
            "oriented_sequence": oriented,
            "cds_start_untrimmed": cds_start,
            "cds_end_untrimmed": cds_end,
            "orientation": orientation,
            "original_length": len(normalize_dna(raw_seq)),
        })
        query_map[qid] = item

    write_fasta(query_records, query_fasta)

    full_db = str(getattr(args, "full_segment_db", "") or "")
    blast_hits: dict[str, list[dict[str, str | int | float]]] = defaultdict(list)
    if full_db and query_records:
        blast_path = run_full_segment_blast(sample_dir, full_db, query_fasta, log_path)
        if blast_path is not None:
            blast_hits = parse_full_segment_blast(blast_path)
    elif not query_records:
        write_empty_genbank_reports(sample_dir, args, log_path, "no_selected_candidates")
        return

    report_rows: list[dict[str, str | int | float | bool]] = []
    coord_rows: list[dict[str, str | int | float | bool]] = []
    clean_records: list[tuple[str, str]] = []
    review_records: list[tuple[str, str]] = []

    for qid in sorted(query_map):
        candidate = query_map[qid]
        genotype = str(candidate.get("genotype", ""))
        segment = str(candidate.get("segment", ""))
        seq = str(candidate["oriented_sequence"])
        cds_start = int(candidate["cds_start_untrimmed"])
        cds_end = int(candidate["cds_end_untrimmed"])
        expected_orf_len = EXPECTED_ORF_LENGTHS.get(segment, 0)
        lb, lr = terminal_homopolymer_run(seq, True)
        rb, rr = terminal_homopolymer_run(seq, False)
        left_trim, right_trim, trim_notes = terminal_trim_counts(seq, int(args.genbank_terminal_homopolymer))
        cleaned = seq[left_trim:len(seq) - right_trim if right_trim else len(seq)]
        clean_cds_start = cds_start - left_trim
        clean_cds_end = cds_end - left_trim
        ambiguous = count_ambiguous_bases(cleaned)
        has_ambiguity = ambiguous > 0
        motif_eval = evaluate_terminal_motifs(
            cleaned,
            segment,
            int(getattr(args, "genbank_terminal_motif_mismatches", 1)),
        )

        orf_seq = ""
        if 1 <= clean_cds_start <= clean_cds_end <= len(cleaned):
            orf_seq = cleaned[clean_cds_start - 1:clean_cds_end]
        orf_val = validate_orf(orf_seq)
        orf_len = int(orf_val["length"])
        orf_diff = orf_len - expected_orf_len if expected_orf_len else ""

        all_hits = blast_hits.get(qid, [])
        best_hit, same_segment_count = choose_best_full_segment_hit(all_hits, segment)

        reason_parts: list[str] = []
        status = "GENBANK_READY_FULL"
        submit = "YES"

        metrics = {
            "reference_aligned_start": "", "reference_aligned_end": "", "reference_coverage_pct": "",
            "contig_aligned_start": "", "contig_aligned_end": "", "contig_coverage_pct": "",
            "missing_ref_5": "", "missing_ref_3": "", "extra_query_5": "", "extra_query_3": "",
        }
        adj_extra5: int | str = ""
        adj_extra3: int | str = ""

        if not full_db:
            status = "NOT_EVALUATED_NO_FULL_SEGMENT_DB"
            submit = "NO"
            reason_parts.append("No --full-segment-db was provided and no default full-segment database was found")
        elif best_hit is None:
            status = "REJECTED_NO_FULL_SEGMENT_REFERENCE_HIT"
            submit = "NO"
            reason_parts.append("No same-segment hit in full-segment reference database")
        else:
            metrics = full_segment_hit_metrics(best_hit)
            adj_extra5 = max(0, int(metrics["extra_query_5"]) - left_trim)
            adj_extra3 = max(0, int(metrics["extra_query_3"]) - right_trim)
            if float(metrics["reference_coverage_pct"]) < float(args.genbank_ref_cov_full):
                status = "PARTIAL_SEGMENT"
                submit = "NO"
                reason_parts.append(
                    f"Reference coverage {float(metrics['reference_coverage_pct']):.2f}% < {float(args.genbank_ref_cov_full):.2f}%"
                )
            if int(metrics["missing_ref_5"]) > int(args.genbank_max_missing_end_bases):
                status = "PARTIAL_SEGMENT"
                submit = "NO"
                reason_parts.append(f"Missing too many 5prime reference bases: {metrics['missing_ref_5']}")
            if int(metrics["missing_ref_3"]) > int(args.genbank_max_missing_end_bases):
                status = "PARTIAL_SEGMENT"
                submit = "NO"
                reason_parts.append(f"Missing too many 3prime reference bases: {metrics['missing_ref_3']}")
            if isinstance(adj_extra5, int) and adj_extra5 > int(args.genbank_max_extra_end_bases):
                if status == "GENBANK_READY_FULL":
                    status = "FULL_BUT_NEEDS_REVIEW"
                    submit = "NO"
                reason_parts.append(f"Extra 5prime contig bases after terminal trimming: {adj_extra5}")
            if isinstance(adj_extra3, int) and adj_extra3 > int(args.genbank_max_extra_end_bases):
                if status == "GENBANK_READY_FULL":
                    status = "FULL_BUT_NEEDS_REVIEW"
                    submit = "NO"
                reason_parts.append(f"Extra 3prime contig bases after terminal trimming: {adj_extra3}")

        if left_trim or right_trim:
            if status == "GENBANK_READY_FULL":
                status = "FULL_BUT_NEEDS_REVIEW"
                submit = "NO"
            reason_parts.append("Terminal bases were trimmed; inspect clean sequence before submission")

        if has_ambiguity and not bool(args.genbank_allow_ambiguity):
            if status == "GENBANK_READY_FULL":
                status = "FULL_BUT_NEEDS_REVIEW"
                submit = "NO"
            reason_parts.append(f"Cleaned sequence contains {ambiguous} ambiguous/non-ACGT bases")

        if not bool(orf_val["valid_orf"]):
            if status == "GENBANK_READY_FULL":
                status = "FULL_BUT_NEEDS_REVIEW"
                submit = "NO"
            reason_parts.append(f"ORF validation failed: {orf_val['orf_validation_status']}")

        orf_length_tolerance = int(getattr(args, "genbank_orf_length_tolerance", 90))
        if expected_orf_len and bool(orf_val["valid_orf"]) and abs(int(orf_diff)) > orf_length_tolerance:
            if status == "GENBANK_READY_FULL":
                status = "FULL_BUT_NEEDS_REVIEW"
                submit = "NO"
            reason_parts.append(
                f"ORF length {orf_len} differs from expected {expected_orf_len} by {orf_diff}, "
                f"which is outside tolerance +/-{orf_length_tolerance}"
            )

        motif_mode = str(getattr(args, "genbank_terminal_motif_mode", "strict")).lower()
        if motif_mode != "off" and not bool(motif_eval["terminal_motifs_pass"]):
            motif_reason = (
                f"Terminal motif check failed: "
                f"5prime {motif_eval['terminal_5_status']} "
                f"observed={motif_eval['terminal_5_observed']} expected={motif_eval['terminal_5_expected']} "
                f"mismatches={motif_eval['terminal_5_mismatches']}; "
                f"3prime {motif_eval['terminal_3_status']} "
                f"observed={motif_eval['terminal_3_observed']} expected={motif_eval['terminal_3_expected']} "
                f"mismatches={motif_eval['terminal_3_mismatches']}"
            )
            if motif_mode == "strict":
                status = "PARTIAL_SEGMENT"
                submit = "NO"
            elif motif_mode == "warn" and status == "GENBANK_READY_FULL":
                status = "FULL_BUT_NEEDS_REVIEW"
                submit = "NO"
            reason_parts.append(motif_reason)

        if not reason_parts:
            reason_parts.append("Passed full-segment reference coverage, terminal motif, terminal artefact, ambiguity, and ORF checks")

        seq_id = f"{sample_dir.name}|{segment}|{genotype}|{status}"
        fasta_header = (
            f">{seq_id} original_contig={candidate.get('qseqid', '')} "
            f"cds={clean_cds_start}..{clean_cds_end} "
            f"full_reference={str(best_hit['sseqid']) if best_hit else ''}"
        )
        if status == "GENBANK_READY_FULL":
            clean_records.append((fasta_header, cleaned))
        else:
            review_records.append((fasta_header + " reason=" + ";".join(reason_parts).replace(" ", "_"), cleaned))

        best_ref = str(best_hit["sseqid"]) if best_hit else ""
        best_ref_genotype = parse_genotype_token(best_ref) if best_ref else ""
        best_ref_segment = segment_from_reference_id(best_ref) if best_ref else ""
        best_ref_len = int(best_hit["slen"]) if best_hit else ""
        best_pident = float(best_hit["pident"]) if best_hit else ""
        best_bitscore = float(best_hit["bitscore"]) if best_hit else ""

        top_full_hits = ranked_same_segment_full_hits(
            all_hits,
            segment,
            int(getattr(args, "genbank_alignment_full_ref_count", 2)),
        )
        alignment_input, alignment_fasta, alignment_ref_count = write_candidate_reference_alignment(
            sample_dir,
            args,
            log_path,
            qid,
            candidate,
            segment,
            genotype,
            status,
            cleaned,
            top_full_hits,
        )

        row = {
            "Sample": sample_dir.name,
            "Segment": segment,
            "Genotype": genotype,
            "Status": status,
            "SubmitToGenBank": submit,
            "Reason": "; ".join(reason_parts),
            "Complete11SegmentSample": "",  # filled below
            "OriginalContigID": candidate.get("qseqid", ""),
            "OriginalContigLength": candidate.get("original_length", ""),
            "CleanLength": len(cleaned),
            "Orientation": candidate.get("orientation", ""),
            "Coverage": candidate.get("coverage", ""),
            "BestFullReference": best_ref,
            "BestFullReferenceGenotype": best_ref_genotype,
            "BestFullReferenceSegment": best_ref_segment,
            "BestFullReferenceLength": best_ref_len,
            "FullReferencePident": best_pident,
            "FullReferenceBitscore": best_bitscore,
            "ReferenceAlignedStart": metrics["reference_aligned_start"],
            "ReferenceAlignedEnd": metrics["reference_aligned_end"],
            "ReferenceCoveragePct": metrics["reference_coverage_pct"],
            "ContigAlignedStart": metrics["contig_aligned_start"],
            "ContigAlignedEnd": metrics["contig_aligned_end"],
            "ContigCoveragePct": metrics["contig_coverage_pct"],
            "MissingReference5primeBases": metrics["missing_ref_5"],
            "MissingReference3primeBases": metrics["missing_ref_3"],
            "ExtraContig5primeBases": metrics["extra_query_5"],
            "ExtraContig3primeBases": metrics["extra_query_3"],
            "AdjustedExtraContig5primeBases": adj_extra5,
            "AdjustedExtraContig3primeBases": adj_extra3,
            "TrimmedLeftBases": left_trim,
            "TrimmedRightBases": right_trim,
            "TerminalTrimNotes": trim_notes,
            "TerminalMotifClassification": motif_eval["terminal_motif_classification"],
            "TerminalMotifsPass": motif_eval["terminal_motifs_pass"],
            "Terminal5Expected": motif_eval["terminal_5_expected"],
            "Terminal5Observed": motif_eval["terminal_5_observed"],
            "Terminal5Mismatches": motif_eval["terminal_5_mismatches"],
            "Terminal5Status": motif_eval["terminal_5_status"],
            "Terminal3Expected": motif_eval["terminal_3_expected"],
            "Terminal3Observed": motif_eval["terminal_3_observed"],
            "Terminal3Mismatches": motif_eval["terminal_3_mismatches"],
            "Terminal3Status": motif_eval["terminal_3_status"],
            "Terminal3Mode": motif_eval.get("terminal_3_mode", ""),
            "Terminal3DefaultExpected": motif_eval.get("terminal_3_default_expected", ""),
            "Terminal3DefaultObserved": motif_eval.get("terminal_3_default_observed", ""),
            "Terminal3DefaultMismatches": motif_eval.get("terminal_3_default_mismatches", ""),
            "Terminal3DefaultStatus": motif_eval.get("terminal_3_default_status", ""),
            "Terminal3NSP1SpecialExpected": motif_eval.get("terminal_3_nsp1_special_expected", ""),
            "Terminal3NSP1SpecialObserved": motif_eval.get("terminal_3_nsp1_special_observed", ""),
            "Terminal3NSP1SpecialMismatches": motif_eval.get("terminal_3_nsp1_special_mismatches", ""),
            "Terminal3NSP1SpecialStatus": motif_eval.get("terminal_3_nsp1_special_status", ""),
            "CDSStart": clean_cds_start,
            "CDSEnd": clean_cds_end,
            "CDSLength": orf_len,
            "ExpectedORFLength": expected_orf_len,
            "ORFLengthDifference": orf_diff,
            "ORFLengthTolerance": int(getattr(args, "genbank_orf_length_tolerance", 90)),
            "ORFLengthWithinTolerance": (
                bool(orf_val["valid_orf"])
                and (not expected_orf_len or abs(int(orf_diff)) <= int(getattr(args, "genbank_orf_length_tolerance", 90)))
            ),
            "ORFValidationStatus": orf_val["orf_validation_status"],
            "StartsWithATG": orf_val["starts_atg"],
            "EndsWithStopCodon": orf_val["ends_stop"],
            "LengthMultipleOf3": orf_val["length_multiple_3"],
            "InternalStopCodons": orf_val["internal_stop_codons"],
            "HasAmbiguity": has_ambiguity,
            "AmbiguousBases": ambiguous,
            "LeftTerminalHomopolymerBase": lb,
            "LeftTerminalHomopolymerLength": lr,
            "RightTerminalHomopolymerBase": rb,
            "RightTerminalHomopolymerLength": rr,
            "ORFBlastSubject": candidate.get("sseqid", ""),
            "ORFBlastPident": candidate.get("pident", ""),
            "ORFCandidateType": candidate.get("orf_candidate_type", ""),
            "FullSegmentDB": full_db,
            "FullSegmentSameSegmentHitCount": same_segment_count,
            "FullSegmentAllHitCount": len(all_hits),
            "AlignmentInputFasta": alignment_input,
            "AlignmentFasta": alignment_fasta,
            "AlignmentReferenceCount": alignment_ref_count,
        }
        report_rows.append(row)

        coord_rows.append({
            "Sample": sample_dir.name,
            "Segment": segment,
            "Genotype": genotype,
            "Status": status,
            "SubmitToGenBank": submit,
            "SequenceID": seq_id,
            "CleanLength": len(cleaned),
            "Feature": "source",
            "Start": 1,
            "End": len(cleaned),
            "Gene": "",
            "Product": "rotavirus A genome segment",
            "Note": "Coordinates are relative to the cleaned submitted sequence",
        })
        coord_rows.append({
            "Sample": sample_dir.name,
            "Segment": segment,
            "Genotype": genotype,
            "Status": status,
            "SubmitToGenBank": submit,
            "SequenceID": seq_id,
            "CleanLength": len(cleaned),
            "Feature": "CDS",
            "Start": clean_cds_start,
            "End": clean_cds_end,
            "Gene": segment,
            "Product": segment,
            "Note": f"ORF validation: {orf_val['orf_validation_status']}",
        })

    ready_segments = {str(r["Segment"]) for r in report_rows if r.get("Status") == "GENBANK_READY_FULL"}
    complete = "YES" if set(SEGMENT_ORDER).issubset(ready_segments) and len(ready_segments) == 11 else "NO"
    for row in report_rows:
        row["Complete11SegmentSample"] = complete

    with report_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=genbank_report_fields())
        writer.writeheader()
        writer.writerows(report_rows)
    with coords_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=genbank_coordinate_fields())
        writer.writeheader()
        writer.writerows(coord_rows)
    write_fasta(clean_records, clean_fasta)
    write_fasta(review_records, review_fasta)
    log(f"Generated {report_path.name} with {len(report_rows)} evaluated sequences", log_path)
    log(f"Generated {clean_fasta.name} with {len(clean_records)} GenBank-ready full sequences", log_path)
    log(f"Generated {review_fasta.name} with {len(review_records)} partial/review sequences", log_path)


# ---------------------------------------------------------------------------
# Full-segment read mapping consensus
# ---------------------------------------------------------------------------

def reference_id_without_leading_genotype(reference_id: str) -> str:
    """Return a reference identifier without a leading genotype token.

    Full-segment DB headers may begin with tokens such as G3|, P8| or T12|.
    Those genotype labels are useful for segment inference but must not be
    carried into mapping consensus headers/metadata because they can be wrong
    when a genotype is absent or mislabelled in the full-reference set.
    """
    ref = str(reference_id or "").strip()
    if not ref:
        return ""
    parts = ref.split("|", 1)
    if len(parts) == 2 and parse_genotype_token(parts[0]):
        return parts[1]
    return ref


def mapping_reference_accession(reference_id: str) -> str:
    """Return a concise accession-like label for human-facing mapping outputs."""
    ref = reference_id_without_leading_genotype(reference_id)
    for token in re.split(r"[|_\s]+", ref):
        token = token.strip()
        if re.fullmatch(r"[A-Za-z]{1,6}\d{4,10}(?:\.\d+)?", token):
            return token
    first = re.split(r"[|\s]+", ref, maxsplit=1)[0] if ref else "reference"
    return safe_file_component(first, max_len=40)


def unique_mapping_reference_name(segment: str, accession: str, used_names: set[str]) -> str:
    base = safe_file_component(f"{segment}_reference_{accession}", max_len=70)
    name = base
    counter = 2
    while name in used_names:
        name = safe_file_component(f"{base}_{counter}", max_len=80)
        counter += 1
    used_names.add(name)
    return name


def normalize_genotype_for_compare(genotype: str) -> str:
    token = str(genotype or "").strip().replace("[", "").replace("]", "")
    token = token.split("_", 1)[0]
    match = re.match(r"^([GPIRCMANTEHgpircmanteh])(\d+|X|\?)", token)
    if not match:
        return token.upper()
    return match.group(1).upper() + match.group(2)


def read_final_orf_genotype_by_segment(sample_dir: Path) -> dict[str, str]:
    """Read the final ORF-based genotype call from rotax_sample_genotyping_report.csv."""
    path = sample_dir / "rotax_sample_genotyping_report.csv"
    out: dict[str, str] = {}
    if not path.exists():
        return out
    with path.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            genotype = row.get("Genotype Constellation", row.get("Genotype", ""))
            for part in genotype.split("-"):
                token = normalize_genotype_for_compare(part)
                if not token:
                    continue
                segment = segment_from_genotype(token)
                if segment:
                    out[segment] = token
            break
    return out


def read_sample_extra_information(sample_dir: Path) -> str:
    """Return the per-sample Extra Information field from the ORF genotype report."""
    path = sample_dir / "rotax_sample_genotyping_report.csv"
    if not path.exists():
        return ""
    with path.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            return str(row.get("Extra Information", "") or "").strip()
    return ""


def sample_has_possible_mixed_infection(sample_dir: Path) -> bool:
    """True when Extra Information suggests additional accepted genotype signals.

    This is deliberately conservative: anything other than blank/None/NA is
    treated as possible mixed infection or otherwise mixed genotype signal.
    Such samples are still mapped, but are not automatically exported as
    GenBank-ready full mapping consensus sequences unless --allow-mixed-genbank
    is used.
    """
    extra = read_sample_extra_information(sample_dir).strip()
    return bool(extra and extra.lower() not in {"none", "na", "n/a", "null", "."})


def selected_mapping_candidates_by_segment(sample_dir: Path) -> dict[str, dict[str, str | int | float]]:
    """Select one SPAdes/ORF candidate per segment for mapping-reference choice.

    If the final ORF-based constellation has a resolved genotype for a segment,
    prefer the candidate that produced that ORF genotype. If the segment is X or
    ? or if the preferred candidate is missing, fall back to the strongest same-
    segment candidate by full/partial status, coverage, alignment length,
    identity and bitscore.
    """
    candidates = read_selected_candidate_rows(sample_dir)
    final_by_segment = read_final_orf_genotype_by_segment(sample_dir)
    grouped: dict[str, list[dict[str, str | int | float]]] = defaultdict(list)
    for genotype, candidate in candidates.items():
        item = dict(candidate)
        item["genotype"] = normalize_genotype_for_compare(str(item.get("genotype", genotype)))
        segment = str(item.get("segment", ""))
        if segment:
            grouped[segment].append(item)

    selected: dict[str, dict[str, str | int | float]] = {}
    for segment, items in grouped.items():
        preferred = final_by_segment.get(segment, "")

        def score(item: dict[str, str | int | float]) -> tuple[int, int, float, int, float, float]:
            genotype = normalize_genotype_for_compare(str(item.get("genotype", "")))
            preferred_match = int(bool(preferred and preferred not in {"GX", "PX", "IX", "RX", "CX", "MX", "AX", "NX", "TX", "EX", "HX"} and not preferred.endswith("?") and genotype == preferred))
            full_match = int(str(item.get("orf_candidate_type", "")) == "full")
            try:
                cov = float(item.get("coverage", 0) or 0)
            except (TypeError, ValueError):
                cov = 0.0
            try:
                length = int(item.get("length", 0) or 0)
            except (TypeError, ValueError):
                length = 0
            try:
                pident = float(item.get("pident", 0) or 0)
            except (TypeError, ValueError):
                pident = 0.0
            try:
                bitscore = float(item.get("bitscore", 0) or 0)
            except (TypeError, ValueError):
                bitscore = 0.0
            return (preferred_match, full_match, cov, length, pident, bitscore)

        selected[segment] = sorted(items, key=score, reverse=True)[0]
    return selected


def build_mapping_query_records_and_map(
    sample_dir: Path,
    possible_mixed: bool = False,
) -> tuple[list[tuple[str, str]], dict[str, dict[str, str | int | float]]]:
    """Build selected SPAdes ORF queries for mapping-reference choice.

    Non-mixed samples use one selected ORF/SPAdes candidate per segment. The
    selected candidate is the one that supports the final genotype call for that
    segment, with a score-based fallback if the final segment is unresolved.

    Possible mixed samples keep all accepted selected ORF candidates. This means
    multiple mapping consensuses for one segment are only produced when multiple
    accepted ORF candidates existed before mapping, not because the same ORF was
    matched to several alternative full-segment references.
    """
    if possible_mixed:
        raw_candidates = read_selected_candidate_rows(sample_dir)
        candidate_mode = "possible_mixed_all_accepted_orf_candidates"
    else:
        raw_candidates = selected_mapping_candidates_by_segment(sample_dir)
        candidate_mode = "non_mixed_one_orf_candidate_per_segment"

    fasta_records = {fasta_key(h): (h, s) for h, s in parse_fasta(sample_dir / "contigs500_COV.fasta")}
    query_records: list[tuple[str, str]] = []
    query_map: dict[str, dict[str, str | int | float]] = {}

    def candidate_sort_key(item: tuple[str, dict[str, str | int | float]]) -> tuple[int, str, str]:
        _key, candidate = item
        segment = str(candidate.get("segment", ""))
        rank = SEGMENT_ORDER.index(segment) if segment in SEGMENT_ORDER else len(SEGMENT_ORDER)
        genotype = normalize_genotype_for_compare(str(candidate.get("genotype", "")))
        return (rank, genotype, str(candidate.get("qseqid", "")))

    for _key, candidate in sorted(raw_candidates.items(), key=candidate_sort_key):
        qseqid = str(candidate.get("qseqid", ""))
        if qseqid not in fasta_records:
            continue
        original_header, raw_seq = fasta_records[qseqid]
        oriented, cds_start, cds_end, orientation = oriented_sequence_and_cds(candidate, raw_seq)
        norm_gt = normalize_genotype_for_compare(str(candidate.get("genotype", _key)))
        segment = str(candidate.get("segment", "")) or segment_from_genotype(norm_gt)
        qid = safe_query_id(sample_dir.name, norm_gt or segment, qseqid)
        # In rare cases a sanitized identifier can collide; append a small suffix.
        base_qid = qid
        suffix = 2
        while qid in query_map:
            qid = f"{base_qid}_{suffix}"
            suffix += 1
        selected_orf = ""
        if 1 <= cds_start <= cds_end <= len(oriented):
            selected_orf = normalize_dna(oriented[cds_start - 1:cds_end])
        if not selected_orf:
            continue
        query_records.append((f">{qid}", selected_orf))
        item = dict(candidate)
        item.update({
            "query_id": qid,
            "candidate_key": qid,
            "candidate_selection_mode": candidate_mode,
            "original_header": original_header,
            "oriented_sequence": oriented,
            "mapping_reference_query_sequence": selected_orf,
            "mapping_reference_query_type": "selected_spades_orf",
            "cds_start_untrimmed": cds_start,
            "cds_end_untrimmed": cds_end,
            "orientation": orientation,
            "original_length": len(normalize_dna(raw_seq)),
            "segment": segment,
            "genotype": norm_gt,
        })
        query_map[qid] = item
    return query_records, query_map

def find_mapping_fastq_pair(sample_dir: Path, args: argparse.Namespace | None = None) -> tuple[Path | None, Path | None, str]:
    """Prefer user-supplied reads, then trimmed reads, then original FASTQs.

    FASTA-input mode creates a synthetic sample folder from contigs only, so no
    reads are available unless the user explicitly supplies them with
    --mapping-r1 and --mapping-r2.
    """
    if args is not None:
        override_r1 = str(getattr(args, "mapping_r1", "") or "").strip()
        override_r2 = str(getattr(args, "mapping_r2", "") or "").strip()
        if override_r1 or override_r2:
            if not (override_r1 and override_r2):
                return None, None, "incomplete_user_supplied_fastq_pair"
            r1 = Path(override_r1).expanduser()
            r2 = Path(override_r2).expanduser()
            if is_nonempty(r1) and is_nonempty(r2):
                return r1.resolve(), r2.resolve(), "user_supplied_fastq"
            return None, None, "user_supplied_fastq_missing_or_empty"

    preferred = [
        (sample_dir / "R1_paired.fastq.gz", sample_dir / "R2_paired.fastq.gz", "trimmomatic_paired"),
        # repair.sh-resynced trimmed reads: a valid, properly-paired mapping input
        # that survives cleanup when R1_paired.fastq.gz is removed. Recognizing it
        # keeps --skip-assembly reruns working on trimmed data.
        (sample_dir / "R1_paired.repaired.fastq.gz", sample_dir / "R2_paired.repaired.fastq.gz", "trimmomatic_paired_repaired"),
        (sample_dir / "deduped_R1.fastq.gz", sample_dir / "deduped_R2.fastq.gz", "clumpify_deduped"),
    ]
    for r1, r2, label in preferred:
        if is_nonempty(r1) and is_nonempty(r2):
            return r1, r2, label
    # No trimmed/deduped reads present. Mapping must NOT silently fall back to
    # the raw, untrimmed FASTQs: that produces a plausible-looking consensus from
    # unfiltered reads and only the ReadSource column reveals it. Refuse by
    # default, and only use raw reads when the user explicitly opts in with
    # --allow-raw-fastq-mapping.
    raw1, raw2 = find_fastq_pair(sample_dir)
    raw_available = (
        raw1 is not None and raw2 is not None
        and is_nonempty(raw1) and is_nonempty(raw2)
    )
    allow_raw = bool(getattr(args, "allow_raw_fastq_mapping", False)) if args is not None else False
    if raw_available and allow_raw:
        return raw1, raw2, "raw_fastq_optin"
    if raw_available:
        # Raw reads exist but trimmed reads were expected and are missing.
        return None, None, "trimmed_reads_missing_raw_fallback_disabled"
    return None, None, "no_fastq_pair"


def run_command_stdout(command: list[str], output_path: Path, log_path: Path, cwd: Path | None = None) -> bool:
    log("Running: " + " ".join(command) + f" > {output_path.name}", log_path)
    with output_path.open("w", encoding="utf-8") as out_handle:
        result = subprocess.run(command, cwd=str(cwd) if cwd else None, stdout=out_handle)
    if result.returncode != 0:
        log(f"Command failed with exit code {result.returncode}: {' '.join(command)}", log_path)
        return False
    return True


def run_shell_command(command: str, log_path: Path, cwd: Path | None = None) -> bool:
    log("Running shell: " + command, log_path)
    result = subprocess.run(["bash", "-lc", command], cwd=str(cwd) if cwd else None)
    if result.returncode != 0:
        log(f"Command failed with exit code {result.returncode}: {command}", log_path)
        return False
    return True


def write_empty_mapping_outputs(sample_dir: Path, reason: str, log_path: Path) -> None:
    """Write only the final, human-facing mapping output files."""
    mapping_dir = sample_dir / "mapping"
    mapping_dir.mkdir(exist_ok=True)
    write_fasta([], mapping_dir / "rotavirus_mapping_consensus.fasta")
    write_fasta([], mapping_dir / "rotavirus_mapping_consensus_full.fasta")
    write_fasta([], mapping_dir / "rotavirus_mapping_consensus_genbank.fasta")
    write_fasta([], mapping_dir / "rotavirus_mapping_consensus_partial_or_review.fasta")
    write_fasta([], mapping_dir / "rotavirus_mapping_orfs.fasta")
    write_fasta([], mapping_dir / "rotavirus_mapping_references.fasta")
    for path, fields in [
        (mapping_dir / "rotavirus_mapping_summary.csv", mapping_summary_fields()),
        (mapping_dir / "rotavirus_mapping_consensus_features.csv", mapping_consensus_feature_fields()),
        (
            mapping_dir / "rotavirus_mapping_consensus_genbank_features.csv",
            mapping_genbank_feature_fields(),
        ),
    ]:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
    with (mapping_dir / "rotavirus_mapping_depth.tsv").open("w", encoding="utf-8") as handle:
        handle.write("Sample\tSegment\tORFCalledGenotype\tMappingReferenceAccession\tMappingReferenceID\tConsensusID\tReferencePosition\tDepth\tMeetsConsensusDepth\tConsensusBase\n")
    log(f"Generated empty mapping outputs for {sample_dir.name}: {reason}", log_path)


def mapping_summary_fields() -> list[str]:
    return [
        "Sample", "PossibleMixedInfection", "ExtraInformation", "Segment", "ORFCalledGenotype",
        "ConsensusID", "ConsensusStatus", "FullConsensusCalled", "GenBankReady",
        "GenBankReadyReason", "StatusReason", "MappingReferenceAccession",
        "MappingReferenceID", "MappingReferenceLength", "MappingReferencePident",
        "MappingReferenceBitscore", "ConsensusLength", "CoveredBases", "CoveragePct", "MeanDepth", "MedianDepth",
        "MinDepthForConsensus", "FirstCoveredPosition", "LastCoveredPosition",
        "MissingReference5primeBases", "MissingReference3primeBases", "ConsensusNBases",
        "MappingTerminal5Status", "MappingTerminal5Observed", "MappingTerminal5Mismatches",
        "MappingTerminal3Status", "MappingTerminal3Observed", "MappingTerminal3Mismatches",
        "MappingTerminalMotifsPass", "MappingTerminalMotifClassification",
        "CDSComplete", "CDSCompletenessReason", "ProteinAlignmentFile", "ORFSource",
        "ORFStart", "ORFEnd", "ORFStrand", "MappingORFLength", "ExpectedORFLength",
        "ORFLengthDifference", "SpadesORFLength", "SpadesVsMappingAlignmentLength",
        "SpadesVsMappingMatches", "SpadesVsMappingMismatches", "SpadesVsMappingIndelColumns",
        "SpadesVsMappingIndelEvents", "SpadesVsMappingDifferenceColumns",
        "SpadesVsMappingIdentityPctIncludingIndels", "ORFComparisonStatus", "ORFComparisonReason",
        "OriginalSpadesContigID", "CandidateKey", "ReadSource",
    ]


def mapping_consensus_feature_fields() -> list[str]:
    return [
        "Sample", "Segment", "ORFCalledGenotype", "ConsensusID", "ConsensusStatus",
        "FullConsensusCalled", "GenBankReady", "MappingReferenceAccession",
        "MappingReferenceID", "Feature", "Start", "End", "Strand", "Gene", "Product", "Note",
    ]


def mapping_genbank_feature_fields() -> list[str]:
    """Fields for features numbered against the terminal-N-trimmed FASTA."""
    return [
        "Sample", "Segment", "ORFCalledGenotype", "ConsensusID",
        "GenBankSequenceID", "ConsensusStatus", "Completeness",
        "FullConsensusCalled", "GenBankReady", "MappingReferenceAccession",
        "MappingReferenceID", "OriginalConsensusLength",
        "Terminal5NsTrimmed", "Terminal3NsTrimmed", "GenBankSequenceLength",
        "Feature", "OriginalStart", "OriginalEnd", "Start", "End", "Strand",
        "Gene", "Product", "Note",
    ]


def genbank_terminal_regions_resolved(
    seq: str,
    segment: str,
    orf_start: object = "",
    orf_end: object = "",
) -> tuple[bool, str]:
    """Require N-free terminal UTRs and conserved terminal-motif windows."""
    s = normalize_dna(seq)
    terminal_5_length = len(TERMINAL_5_MOTIF)
    terminal_3_length = (
        len(TERMINAL_3_MOTIF_NSP1_SPECIAL)
        if segment.upper() == "NSP1"
        else len(TERMINAL_3_MOTIF_DEFAULT)
    )

    try:
        cds_start = int(str(orf_start))
        cds_end = int(str(orf_end))
    except (TypeError, ValueError):
        cds_start = 0
        cds_end = 0

    if 1 <= cds_start <= cds_end <= len(s):
        terminal_5_end = max(cds_start - 1, terminal_5_length)
        terminal_3_start = min(cds_end, max(0, len(s) - terminal_3_length))
    else:
        terminal_5_end = terminal_5_length
        terminal_3_start = max(0, len(s) - terminal_3_length)

    terminal_5 = s[:terminal_5_end]
    terminal_3 = s[terminal_3_start:]
    reasons: list[str] = []
    if len(s) < terminal_5_length or "N" in terminal_5:
        reasons.append("5-prime terminal region is incomplete or contains N")
    if len(s) < terminal_3_length or "N" in terminal_3:
        reasons.append("3-prime terminal region is incomplete or contains N")

    if reasons:
        return False, "; ".join(reasons)
    return True, "5-prime and 3-prime terminal regions are N-free"


def trim_terminal_ns(seq: str) -> tuple[str, int, int]:
    """Remove only terminal N calls and return sequence, left trim, right trim."""
    normalized = normalize_dna(seq)
    left_trim = len(normalized) - len(normalized.lstrip("N"))
    right_trim = len(normalized) - len(normalized.rstrip("N"))
    if left_trim + right_trim > len(normalized):
        right_trim = max(0, len(normalized) - left_trim)
    end = len(normalized) - right_trim if right_trim else len(normalized)
    return normalized[left_trim:end], left_trim, right_trim


def genbank_mapping_export(
    sample: str,
    consensus_records: list[tuple[str, str]],
    summary_rows: list[dict[str, str]],
) -> tuple[list[tuple[str, str]], dict[str, dict[str, object]]]:
    """Build the trimmed GenBank FASTA and coordinate-remapping metadata."""
    summary_by_id = {
        row.get("ConsensusID", ""): row
        for row in summary_rows
        if row.get("ConsensusID", "")
    }
    id_counts: defaultdict[str, int] = defaultdict(int)
    records: list[tuple[str, str]] = []
    trim_metadata: dict[str, dict[str, object]] = {}

    for header, seq in consensus_records:
        key = fasta_key(header)
        key_parts = key.split("|")
        row = summary_by_id.get(key, {})
        record_sample = str(row.get("Sample", "") or sample)
        segment = str(
            row.get("Segment", "")
            or (key_parts[1] if len(key_parts) > 1 else "segment")
        )
        terminals_resolved, _terminal_reason = genbank_terminal_regions_resolved(
            seq,
            segment,
            row.get("ORFStart", ""),
            row.get("ORFEnd", ""),
        )
        completeness = (
            "complete"
            if (
                str(row.get("FullConsensusCalled", "")).upper() == "YES"
                and terminals_resolved
            )
            else "partial"
        )

        base_id = f"{record_sample}_{segment}"
        id_counts[base_id] += 1
        occurrence = id_counts[base_id]
        record_id = base_id if occurrence == 1 else f"{base_id}_{occurrence}"
        trimmed_seq, left_trim, right_trim = trim_terminal_ns(seq)
        records.append(
            (f">{record_id} [completeness={completeness}]", trimmed_seq)
        )
        trim_metadata[key] = {
            "GenBankSequenceID": record_id,
            "Completeness": completeness,
            "OriginalConsensusLength": len(normalize_dna(seq)),
            "Terminal5NsTrimmed": left_trim,
            "Terminal3NsTrimmed": right_trim,
            "GenBankSequenceLength": len(trimmed_seq),
        }

    return records, trim_metadata


def genbank_mapping_records(
    sample: str,
    consensus_records: list[tuple[str, str]],
    summary_rows: list[dict[str, str]],
) -> list[tuple[str, str]]:
    """Return terminal-N-trimmed per-sample GenBank FASTA records."""
    records, _trim_metadata = genbank_mapping_export(
        sample,
        consensus_records,
        summary_rows,
    )
    return records


def genbank_mapping_feature_rows(
    feature_rows: list[dict[str, object]],
    trim_metadata: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    """Shift feature coordinates to match the terminal-N-trimmed FASTA."""
    adjusted_rows: list[dict[str, object]] = []
    for row in feature_rows:
        consensus_id = str(row.get("ConsensusID", "") or "")
        trim = trim_metadata.get(consensus_id)
        if trim is None:
            continue
        try:
            original_start = int(str(row.get("Start", "") or ""))
            original_end = int(str(row.get("End", "") or ""))
            left_trim = int(trim["Terminal5NsTrimmed"])
            trimmed_length = int(trim["GenBankSequenceLength"])
        except (TypeError, ValueError):
            continue

        feature = str(row.get("Feature", "") or "")
        if feature == "source":
            adjusted_start = 1
            adjusted_end = trimmed_length
        else:
            adjusted_start = original_start - left_trim
            adjusted_end = original_end - left_trim
            if adjusted_end < 1 or adjusted_start > trimmed_length:
                continue
            adjusted_start = max(1, adjusted_start)
            adjusted_end = min(trimmed_length, adjusted_end)

        adjusted_rows.append({
            "Sample": row.get("Sample", ""),
            "Segment": row.get("Segment", ""),
            "ORFCalledGenotype": row.get("ORFCalledGenotype", ""),
            "ConsensusID": consensus_id,
            "GenBankSequenceID": trim["GenBankSequenceID"],
            "ConsensusStatus": row.get("ConsensusStatus", ""),
            "Completeness": trim["Completeness"],
            "FullConsensusCalled": row.get("FullConsensusCalled", ""),
            "GenBankReady": row.get("GenBankReady", ""),
            "MappingReferenceAccession": row.get("MappingReferenceAccession", ""),
            "MappingReferenceID": row.get("MappingReferenceID", ""),
            "OriginalConsensusLength": trim["OriginalConsensusLength"],
            "Terminal5NsTrimmed": trim["Terminal5NsTrimmed"],
            "Terminal3NsTrimmed": trim["Terminal3NsTrimmed"],
            "GenBankSequenceLength": trimmed_length,
            "Feature": feature,
            "OriginalStart": original_start,
            "OriginalEnd": original_end,
            "Start": adjusted_start,
            "End": adjusted_end,
            "Strand": row.get("Strand", ""),
            "Gene": row.get("Gene", ""),
            "Product": row.get("Product", ""),
            "Note": row.get("Note", ""),
        })

    return adjusted_rows


def write_mapping_genbank_exports(
    mapping_dir: Path,
    sample: str,
    consensus_records: list[tuple[str, str]],
    summary_rows: list[dict[str, object]],
    feature_rows: list[dict[str, object]],
) -> tuple[int, int]:
    """Write matching trimmed GenBank FASTA and adjusted-coordinate features."""
    genbank_records, trim_metadata = genbank_mapping_export(
        sample,
        consensus_records,
        summary_rows,
    )
    adjusted_features = genbank_mapping_feature_rows(
        feature_rows,
        trim_metadata,
    )
    write_fasta(
        genbank_records,
        mapping_dir / "rotavirus_mapping_consensus_genbank.fasta",
    )
    with (
        mapping_dir / "rotavirus_mapping_consensus_genbank_features.csv"
    ).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=mapping_genbank_feature_fields(),
        )
        writer.writeheader()
        writer.writerows(adjusted_features)
    return len(genbank_records), len(adjusted_features)


def parse_depth_file(depth_path: Path, reference_lengths: dict[str, int], min_depth: int) -> dict[str, dict[str, int | float | str]]:
    depths: dict[str, list[int]] = {ref: [0] * length for ref, length in reference_lengths.items()}
    if depth_path.exists():
        with depth_path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                ref, pos_text, depth_text = parts[:3]
                if ref not in depths:
                    continue
                try:
                    pos = int(pos_text)
                    depth = int(depth_text)
                except ValueError:
                    continue
                if 1 <= pos <= len(depths[ref]):
                    depths[ref][pos - 1] = depth

    stats: dict[str, dict[str, int | float | str]] = {}
    for ref, values in depths.items():
        length = len(values)
        covered_positions = [idx + 1 for idx, depth in enumerate(values) if depth >= min_depth]
        covered = len(covered_positions)
        first = covered_positions[0] if covered_positions else ""
        last = covered_positions[-1] if covered_positions else ""
        missing5 = int(first) - 1 if isinstance(first, int) else length
        missing3 = length - int(last) if isinstance(last, int) else length
        sorted_values = sorted(values)
        if values:
            mean_depth = sum(values) / len(values)
            mid = len(sorted_values) // 2
            if len(sorted_values) % 2:
                median_depth = float(sorted_values[mid])
            else:
                median_depth = (sorted_values[mid - 1] + sorted_values[mid]) / 2.0
        else:
            mean_depth = 0.0
            median_depth = 0.0
        stats[ref] = {
            "reference_length": length,
            "covered_bases": covered,
            "coverage_pct": 100.0 * covered / length if length else 0.0,
            "mean_depth": mean_depth,
            "median_depth": median_depth,
            "first_covered": first,
            "last_covered": last,
            "missing5": missing5,
            "missing3": missing3,
        }
    return stats


def read_depth_vectors(depth_path: Path, reference_lengths: dict[str, int]) -> dict[str, list[int]]:
    """Return per-reference raw read-depth vectors in mapping coordinates."""
    depths: dict[str, list[int]] = {
        ref_name: [0] * length for ref_name, length in reference_lengths.items()
    }
    if not depth_path.exists():
        return depths
    with depth_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            ref_name, pos_text, depth_text = parts[:3]
            if ref_name not in depths:
                continue
            try:
                pos = int(pos_text)
                depth = int(depth_text)
            except ValueError:
                continue
            if 1 <= pos <= len(depths[ref_name]):
                depths[ref_name][pos - 1] = depth
    return depths


def write_low_coverage_bed(depth_path: Path, reference_lengths: dict[str, int], min_depth: int, bed_path: Path) -> None:
    depths: dict[str, list[int]] = {ref: [0] * length for ref, length in reference_lengths.items()}
    if depth_path.exists():
        with depth_path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                ref, pos_text, depth_text = parts[:3]
                if ref not in depths:
                    continue
                try:
                    pos = int(pos_text)
                    depth = int(depth_text)
                except ValueError:
                    continue
                if 1 <= pos <= len(depths[ref]):
                    depths[ref][pos - 1] = depth

    with bed_path.open("w", encoding="utf-8") as out:
        for ref, values in depths.items():
            start: int | None = None
            for idx, depth in enumerate(values):
                if depth < min_depth and start is None:
                    start = idx
                elif depth >= min_depth and start is not None:
                    out.write(f"{ref}\t{start}\t{idx}\n")
                    start = None
            if start is not None:
                out.write(f"{ref}\t{start}\t{len(values)}\n")


def consensus_record_status(
    depth_stats: dict[str, int | float | str],
    args: argparse.Namespace,
) -> tuple[str, str, list[str]]:
    reasons: list[str] = []
    cov = float(depth_stats.get("coverage_pct", 0.0) or 0.0)
    missing5 = int(depth_stats.get("missing5", 0) or 0)
    missing3 = int(depth_stats.get("missing3", 0) or 0)
    if cov < float(getattr(args, "mapping_ref_cov_full", 99.0)):
        reasons.append(f"Reference coverage {cov:.2f}% < {float(getattr(args, 'mapping_ref_cov_full', 99.0)):.2f}%")
    if missing5 > int(getattr(args, "mapping_max_missing_end_bases", 5)):
        reasons.append(f"Missing too many 5prime reference bases: {missing5}")
    if missing3 > int(getattr(args, "mapping_max_missing_end_bases", 5)):
        reasons.append(f"Missing too many 3prime reference bases: {missing3}")
    if reasons:
        return "PARTIAL_MAPPING_CONSENSUS", "partial", reasons
    return "FULL_MAPPING_CONSENSUS", "full", ["Mapping consensus covers full reference within thresholds"]


def best_segment_blast_hits_for_mapping(
    sample_dir: Path,
    args: argparse.Namespace,
    query_records: list[tuple[str, str]],
    query_map: dict[str, dict[str, str | int | float]],
    log_path: Path,
) -> dict[str, list[dict[str, str | int | float]]]:
    """Return the closest same-segment full-reference hit for each selected ORF.

    This BLAST is deliberately separate from the optional whole-contig
    full-segment QC. Only the extracted selected SPAdes ORF sequence is used to
    choose the mapping reference.
    """
    query_fasta = sample_dir / "rotavirus_mapping_orf_reference_query.fasta"
    blast_path = sample_dir / "rotavirus_mapping_orf_reference_blast.tsv"
    write_fasta(query_records, query_fasta)
    full_db = str(getattr(args, "full_segment_db", "") or "")
    blast_hits: dict[str, list[dict[str, str | int | float]]] = defaultdict(list)

    if full_db and query_records:
        ok = run_command(
            [
                "blastn",
                "-query", query_fasta.name,
                "-db", full_db,
                "-out", blast_path.name,
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand slen qlen",
                "-max_target_seqs", "50",
                "-max_hsps", "1",
                "-strand", "both",
            ],
            log_path,
            cwd=sample_dir,
        )
        if ok and blast_path.exists():
            blast_hits = parse_full_segment_blast(blast_path)

    requested_count = int(getattr(args, "mapping_full_ref_count", 1) or 1)
    if requested_count != 1:
        log(
            f"{sample_dir.name}: --mapping-full-ref-count={requested_count} was supplied, but mapping uses exactly "
            "one closest full-segment reference per selected ORF candidate.",
            log_path,
        )

    selected_hits: dict[str, list[dict[str, str | int | float]]] = {}
    for qid, candidate in query_map.items():
        segment = str(candidate.get("segment", ""))
        selected_hits[qid] = ranked_same_segment_full_hits(blast_hits.get(qid, []), segment, limit=1)
    return selected_hits


def extract_spades_orf_for_candidate(candidate: dict[str, str | int | float]) -> str:
    seq = normalize_dna(str(candidate.get("oriented_sequence", "")))
    try:
        start = int(candidate.get("cds_start_untrimmed", 0) or 0)
        end = int(candidate.get("cds_end_untrimmed", 0) or 0)
    except (TypeError, ValueError):
        return ""
    if 1 <= start <= end <= len(seq):
        return seq[start - 1:end]
    return ""


def build_mapping_orf_alignment_subject(consensus_records: list[tuple[str, str]], output_path: Path) -> dict[str, tuple[str, str]]:
    by_segment: dict[str, tuple[str, str]] = {}
    for header, seq in consensus_records:
        key = fasta_key(header)
        parts = key.split("|")
        segment = ""
        for part in parts:
            if part in SEGMENT_ORDER:
                segment = part
                break
        if segment:
            by_segment[segment] = (header, normalize_dna(seq))
    write_fasta(consensus_records, output_path)
    return by_segment


def _segment_from_subject_title(stitle: str) -> str:
    """Pull the segment tag from a genotyping-DB defline.

    Genotyping-DB deflines look like::

        G1|KC579514:28..1008|USA|1989 orf_crop_source=... source=VP7

    The segment is the value of the LAST ``source=`` field (an earlier
    ``orf_crop_source=`` also matches a naive search).
    """
    matches = re.findall(r"source=(\w+)", stitle or "")
    return matches[-1] if matches else ""


def run_genotype_orf_blast(
    mapping_dir: Path,
    consensus_fasta: Path,
    blast_db: str,
    segment_by_consensus_id: dict[str, str],
    log_path: Path,
) -> dict[str, dict[str, str | int | float]]:
    """Return the best same-segment genotyping-DB ORF hit per consensus.

    Mapping consensuses are required to be positive-sense sequences, so this
    search is deliberately restricted to the plus strand. Subjects in the
    genotyping database are validated, coding-cropped ORFs. Hits are filtered to
    the expected segment and ranked by BLAST bitscore first and percent identity
    second. The selected subject is later fetched, translated, and aligned to all
    three forward translations of the consensus to determine the CDS coordinates.
    """
    if not blast_db or not consensus_fasta.exists() or not parse_fasta(consensus_fasta):
        return {}

    with tempfile.TemporaryDirectory(prefix="rotax_genotype_orf_", dir=mapping_dir) as tmp_name:
        tmp_dir = Path(tmp_name)
        blast_tsv = tmp_dir / "genotype_orf_coordinates.tsv"
        ok = run_command(
            [
                "blastn",
                "-query", str(consensus_fasta),
                "-db", blast_db,
                "-out", str(blast_tsv),
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qlen slen stitle",
                "-max_target_seqs", "50",
                "-max_hsps", "1",
                "-strand", "plus",
            ],
            log_path,
        )
        if not ok or not blast_tsv.exists():
            return {}

        fields = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "sstrand", "qlen", "slen", "stitle",
        ]
        best: dict[str, dict[str, str | int | float]] = {}
        with blast_tsv.open("r", encoding="utf-8", errors="replace") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for parts in reader:
                if len(parts) < len(fields):
                    continue
                if len(parts) > len(fields):
                    parts = parts[:len(fields) - 1] + ["\t".join(parts[len(fields) - 1:])]
                row = dict(zip(fields, parts))
                query_key = row["qseqid"]
                want_segment = segment_by_consensus_id.get(query_key, "")
                subject_segment = (
                    _segment_from_subject_title(row.get("stitle", ""))
                    or segment_from_reference_id(row.get("sseqid", ""))
                )
                if want_segment and subject_segment != want_segment:
                    continue
                try:
                    hit: dict[str, str | int | float] = {
                        "qseqid": row["qseqid"],
                        "sseqid": row["sseqid"],
                        "pident": float(row["pident"]),
                        "length": int(row["length"]),
                        "qstart": int(row["qstart"]),
                        "qend": int(row["qend"]),
                        "sstart": int(row["sstart"]),
                        "send": int(row["send"]),
                        "bitscore": float(row["bitscore"]),
                        "sstrand": row["sstrand"],
                        "qlen": int(row["qlen"]),
                        "slen": int(row["slen"]),
                        "subject_segment": subject_segment,
                    }
                except ValueError:
                    continue

                current = best.get(query_key)
                hit_rank = (float(hit["bitscore"]), float(hit["pident"]))
                current_rank = (
                    (float(current["bitscore"]), float(current["pident"]))
                    if current is not None
                    else None
                )
                if current_rank is None or hit_rank > current_rank:
                    best[query_key] = hit

        log(
            f"{mapping_dir.parent.name}: selected best same-segment genotyping-DB ORF hits "
            f"for {len(best)} of {len(segment_by_consensus_id)} mapping consensus sequence(s)",
            log_path,
        )
        return best

def translate_dna(seq: str) -> str:
    """Translate a nucleotide sequence in frame 0.

    Codons containing N or another ambiguous symbol become X. Stop codons are
    represented by ``*``. The caller selects the reading frame by slicing the
    nucleotide sequence before calling this function.
    """
    s = normalize_dna(seq)
    return "".join(GENETIC_CODE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def align_two_proteins(
    reference_protein: str,
    frame_protein: str,
    mapping_dir: Path,
    aligner: str,
    log_path: Path,
) -> tuple[str, str] | None:
    """Pairwise-align a reference protein and one consensus-frame translation."""
    aligner = (aligner or "").strip()
    if not reference_protein or not frame_protein:
        return None
    if not aligner or aligner.lower() == "none":
        log("Protein-guided CDS annotation is disabled because the alignment command is 'none'", log_path)
        return None
    if shutil.which(aligner) is None:
        log(f"Protein-guided CDS annotation requires {aligner}, but it was not found in PATH", log_path)
        return None

    with tempfile.TemporaryDirectory(prefix="rotax_protein_orf_", dir=mapping_dir) as tmp_name:
        tmp_dir = Path(tmp_name)
        input_fasta = tmp_dir / "protein_pair.fasta"
        aligned_fasta = tmp_dir / "protein_pair.aln.fasta"
        write_fasta([
            (">REFERENCE_ORF", reference_protein),
            (">CONSENSUS_FRAME", frame_protein),
        ], input_fasta)
        command = [
            aligner,
            "--quiet",
            "--anysymbol",
            "--globalpair",
            "--maxiterate", "1000",
            str(input_fasta),
        ]
        with aligned_fasta.open("w", encoding="utf-8") as out_handle:
            result = subprocess.run(
                command,
                stdout=out_handle,
                stderr=subprocess.PIPE,
                text=True,
            )
        if result.returncode != 0 or not is_nonempty(aligned_fasta):
            stderr = (result.stderr or "").strip()
            log(
                f"Protein alignment failed with exit code {result.returncode}"
                + (f": {stderr}" if stderr else ""),
                log_path,
            )
            return None

        records = parse_fasta(aligned_fasta)
        aligned = {fasta_key(header): seq for header, seq in records}
        ref_aln = aligned.get("REFERENCE_ORF", "")
        frame_aln = aligned.get("CONSENSUS_FRAME", "")
        if not ref_aln or not frame_aln or len(ref_aln) != len(frame_aln):
            return None
        return ref_aln, frame_aln


def protein_alignment_placement(ref_aln: str, frame_aln: str) -> dict[str, int | float] | None:
    """Score an aligned frame and locate the reference-protein span within it.

    ``missing_amino_acids`` counts reference residues aligned to a gap in the
    consensus translation. ``extra_amino_acids`` counts consensus residues
    aligned to a gap in the reference. These are amino-acid/codon differences,
    not nucleotide-level frameshifts.
    """
    if not ref_aln or len(ref_aln) != len(frame_aln):
        return None
    ref_columns = [i for i, residue in enumerate(ref_aln) if residue != "-"]
    if not ref_columns:
        return None

    first_col = ref_columns[0]
    last_col = ref_columns[-1]
    frame_start_residue = sum(1 for residue in frame_aln[:first_col] if residue != "-")
    frame_end_residue_exclusive = sum(
        1 for residue in frame_aln[:last_col + 1] if residue != "-"
    )
    span_ref = ref_aln[first_col:last_col + 1]
    span_frame = frame_aln[first_col:last_col + 1]

    missing_amino_acids = sum(
        1 for ref_residue, frame_residue in zip(span_ref, span_frame)
        if ref_residue != "-" and frame_residue == "-"
    )
    extra_amino_acids = sum(
        1 for ref_residue, frame_residue in zip(span_ref, span_frame)
        if ref_residue == "-" and frame_residue != "-"
    )
    internal_gap_columns = missing_amino_acids + extra_amino_acids
    comparable = sum(
        1 for ref_residue, frame_residue in zip(span_ref, span_frame)
        if ref_residue not in {"-", "X", "*"}
        and frame_residue not in {"-", "X", "*"}
    )
    matches = sum(
        1 for ref_residue, frame_residue in zip(span_ref, span_frame)
        if ref_residue == frame_residue and ref_residue not in {"-", "X", "*"}
    )
    reference_residues = sum(1 for residue in span_ref if residue != "-")
    frame_residues_in_span = sum(1 for residue in span_frame if residue != "-")
    identity = matches / comparable if comparable else 0.0
    comparable_fraction = comparable / reference_residues if reference_residues else 0.0

    return {
        "frame_start_residue": frame_start_residue,
        "frame_end_residue_exclusive": frame_end_residue_exclusive,
        "internal_gap_columns": internal_gap_columns,
        "missing_amino_acids": missing_amino_acids,
        "extra_amino_acids": extra_amino_acids,
        "matches": matches,
        "comparable": comparable,
        "reference_residues": reference_residues,
        "frame_residues_in_span": frame_residues_in_span,
        "identity": identity,
        "comparable_fraction": comparable_fraction,
    }

def annotate_consensus_cds_by_protein(
    consensus_seq: str,
    genotype_hit: dict[str, str | int | float] | None,
    blast_db: str,
    segment: str,
    consensus_id: str,
    mapping_dir: Path,
    aligner: str,
    log_path: Path,
    depth_values: list[int] | None = None,
) -> dict[str, object]:
    """Locate and quality-check a CDS by reference-protein frame alignment.

    The consensus is assumed to be positive sense. The best same-segment ORF hit
    from the genotyping BLAST database is fetched and translated. The three
    forward reading frames of the mapping consensus are translated and aligned
    separately to the reference protein. The best alignment is always retained
    in ``mapping/protein_alignments`` when alignment succeeds, including when a
    later quality rule prevents CDS annotation.

    N-containing codons become X and are ignored for protein identity. One to
    three missing or extra amino acids are accepted as a small in-frame protein
    difference and are explicitly documented. The consensus CDS endpoint is
    derived from the number of consensus codons in the selected alignment and a
    nearby in-frame terminal stop, rather than being forced to the nucleotide
    length of the reference ORF.
    """
    result: dict[str, object] = {
        "annotate": False,
        "start": "",
        "end": "",
        "strand": "+",
        "cds_nt": "",
        "orf_source": "not_annotated",
        "cds_complete": "NO",
        "reason": "",
        "frame": "",
        "frame_identity": 0.0,
        "comparable_fraction": 0.0,
        "protein_alignment_file": "",
        "missing_amino_acids": 0,
        "extra_amino_acids": 0,
    }

    consensus = normalize_dna(consensus_seq)
    if len(consensus) < 6:
        result["reason"] = "Consensus is too short to contain a CDS; CDS not annotated"
        return result
    if genotype_hit is None:
        result["reason"] = "No same-segment genotyping-database ORF hit was found; CDS not annotated"
        return result
    if str(genotype_hit.get("sstrand", "plus")).strip().lower() == "minus":
        result["reason"] = (
            "The best genotyping-database ORF hit is on the minus strand, although mapping consensuses "
            "must be positive sense; CDS not annotated"
        )
        return result

    reference_id = str(genotype_hit.get("sseqid", ""))
    fetched = fetch_blastdb_sequence(blast_db, reference_id, log_path)
    if fetched is None:
        result["reason"] = (
            f"The nucleotide sequence of best genotyping-database ORF "
            f"{reference_id or '[unknown]'} could not be retrieved; CDS not annotated"
        )
        return result

    _reference_header, reference_orf_nt = fetched
    reference_orf_nt = normalize_dna(reference_orf_nt)
    if len(reference_orf_nt) < 6 or len(reference_orf_nt) % 3 != 0:
        result["reason"] = (
            f"Best genotyping-database ORF {reference_id} has invalid nucleotide length "
            f"{len(reference_orf_nt)}; CDS not annotated"
        )
        return result

    reference_protein = translate_dna(reference_orf_nt)
    if not reference_protein.endswith("*"):
        result["reason"] = (
            f"Best genotyping-database ORF {reference_id} does not end in a translated stop codon; "
            "CDS not annotated"
        )
        return result
    reference_protein_core = reference_protein[:-1]
    if not reference_protein_core or "*" in reference_protein_core:
        result["reason"] = (
            f"Best genotyping-database ORF {reference_id} contains an invalid internal stop; "
            "CDS not annotated"
        )
        return result

    best: dict[str, int | float] | None = None
    best_alignment: tuple[str, str] | None = None
    for frame in (0, 1, 2):
        # Convert stop codons to X for alignment. Fully called premature stops are
        # checked directly in the selected nucleotide span below.
        frame_protein = translate_dna(consensus[frame:]).replace("*", "X")
        aligned = align_two_proteins(
            reference_protein_core,
            frame_protein,
            mapping_dir,
            aligner,
            log_path,
        )
        if aligned is None:
            continue
        placement = protein_alignment_placement(aligned[0], aligned[1])
        if placement is None:
            continue
        placement["frame"] = frame
        candidate_rank = (
            int(placement["matches"]),
            int(placement["comparable"]),
            -int(placement["internal_gap_columns"]),
            float(placement["identity"]),
        )
        if best is None:
            best = placement
            best_alignment = aligned
        else:
            best_rank = (
                int(best["matches"]),
                int(best["comparable"]),
                -int(best["internal_gap_columns"]),
                float(best["identity"]),
            )
            if candidate_rank > best_rank:
                best = placement
                best_alignment = aligned

    if best is None or best_alignment is None:
        result["reason"] = (
            "None of the three forward consensus translations could be aligned to the best "
            "genotyping-database ORF; CDS not annotated"
        )
        return result

    frame = int(best["frame"])
    matches = int(best["matches"])
    comparable = int(best["comparable"])
    identity = float(best["identity"])
    comparable_fraction = float(best["comparable_fraction"])
    internal_gap_columns = int(best["internal_gap_columns"])
    missing_amino_acids = int(best["missing_amino_acids"])
    extra_amino_acids = int(best["extra_amino_acids"])
    result["frame"] = frame + 1
    result["frame_identity"] = round(identity, 6)
    result["comparable_fraction"] = round(comparable_fraction, 6)
    result["missing_amino_acids"] = missing_amino_acids
    result["extra_amino_acids"] = extra_amino_acids

    protein_alignment_dir = mapping_dir / "protein_alignments"
    protein_alignment_dir.mkdir(exist_ok=True)
    alignment_stem = safe_file_component(f"{segment}_{consensus_id}", 180)
    protein_alignment_path = protein_alignment_dir / f"{alignment_stem}.protein.aln.fasta"
    write_fasta([
        (
            f">REFERENCE_ORF|segment={segment}|id={reference_id}",
            best_alignment[0],
        ),
        (
            f">MAPPING_CONSENSUS_TRANSLATION|segment={segment}|id={consensus_id}|"
            f"positive_sense_frame={frame + 1}",
            best_alignment[1],
        ),
    ], protein_alignment_path)
    relative_alignment_path = str(protein_alignment_path.relative_to(mapping_dir))
    result["protein_alignment_file"] = relative_alignment_path
    log(
        f"{segment}: saved selected protein alignment to {protein_alignment_path}",
        log_path,
    )

    if comparable == 0:
        result["reason"] = (
            "The candidate CDS contains no fully called codons that can verify the reading frame; "
            f"CDS not annotated. Protein alignment: {relative_alignment_path}"
        )
        return result
    if comparable_fraction < PROTEIN_ORF_MIN_COMPARABLE_FRACTION:
        result["reason"] = (
            f"Only {comparable_fraction * 100:.1f}% of reference amino-acid positions are comparable "
            f"(< {PROTEIN_ORF_MIN_COMPARABLE_FRACTION * 100:.0f}%); the consensus is too N-masked "
            f"to verify a complete CDS, so it was not annotated. Protein alignment: "
            f"{relative_alignment_path}"
        )
        return result
    if identity < PROTEIN_ORF_MIN_IDENTITY:
        result["reason"] = (
            f"Best-frame protein identity to {reference_id} is {identity * 100:.1f}% "
            f"(< {PROTEIN_ORF_MIN_IDENTITY * 100:.0f}%); likely wrong frame or low-quality sequence, "
            f"so the CDS was not annotated. Protein alignment: {relative_alignment_path}"
        )
        return result
    if internal_gap_columns > PROTEIN_ORF_MAX_SMALL_INDEL_AA:
        result["reason"] = (
            f"The best protein alignment to {reference_id} contains {internal_gap_columns} amino-acid "
            f"gap column(s): {missing_amino_acids} missing and {extra_amino_acids} extra in the mapping "
            f"consensus relative to the reference. More than {PROTEIN_ORF_MAX_SMALL_INDEL_AA} total "
            f"missing/extra amino acids is not accepted; CDS not annotated. Protein alignment: "
            f"{relative_alignment_path}"
        )
        return result

    alignment_start0 = frame + 3 * int(best["frame_start_residue"])
    alignment_coding_end0 = frame + 3 * int(best["frame_end_residue_exclusive"])

    def codon_at(position0: int) -> str:
        if position0 < 0 or position0 + 3 > len(consensus):
            return ""
        return consensus[position0:position0 + 3]

    # A terminal 1-3-aa difference may place the actual ATG or stop up to nine
    # nucleotides from the alignment-derived boundary. Prefer the exact boundary,
    # then the nearest upstream/downstream in-frame codon.
    nearby_codon_offsets = [0, -3, 3, -6, 6, -9, 9]
    start0: int | None = None
    for offset in nearby_codon_offsets:
        candidate = alignment_start0 + offset
        if codon_at(candidate) == "ATG":
            start0 = candidate
            break
    if start0 is None:
        aligned_start_codon = codon_at(alignment_start0)
        if aligned_start_codon and any(base not in "ACGT" for base in aligned_start_codon):
            start0 = alignment_start0
        else:
            result["reason"] = (
                f"No in-frame ATG was found within three codons of the protein-alignment start "
                f"({alignment_start0 + 1}); CDS not annotated. Protein alignment: "
                f"{relative_alignment_path}"
            )
            return result

    stop_start0: int | None = None
    for offset in nearby_codon_offsets:
        candidate = alignment_coding_end0 + offset
        if codon_at(candidate) in STOP_CODONS:
            stop_start0 = candidate
            break
    if stop_start0 is None:
        aligned_stop_codon = codon_at(alignment_coding_end0)
        if aligned_stop_codon and any(base not in "ACGT" for base in aligned_stop_codon):
            stop_start0 = alignment_coding_end0
        else:
            result["reason"] = (
                f"No in-frame stop codon was found within three codons of the protein-alignment end "
                f"({alignment_coding_end0 + 1}); CDS not annotated. Protein alignment: "
                f"{relative_alignment_path}"
            )
            return result

    end0 = stop_start0 + 3
    if start0 < 0 or end0 > len(consensus) or end0 <= start0:
        result["reason"] = (
            "The complete protein-guided CDS span is not present in the mapping consensus "
            f"(5-prime or 3-prime truncation); CDS not annotated. Protein alignment: "
            f"{relative_alignment_path}"
        )
        return result
    if (end0 - start0) % 3 != 0:
        result["reason"] = (
            "The inferred CDS boundaries are not in the same reading frame; CDS not annotated. "
            f"Protein alignment: {relative_alignment_path}"
        )
        return result

    cds_nt = consensus[start0:end0]
    reference_core_aa_length = len(reference_protein_core)
    consensus_core_aa_length = len(cds_nt) // 3 - 1
    aa_length_difference = consensus_core_aa_length - reference_core_aa_length
    if abs(aa_length_difference) > PROTEIN_ORF_MAX_SMALL_INDEL_AA:
        direction = "extra" if aa_length_difference > 0 else "missing"
        result["reason"] = (
            f"The inferred CDS encodes {consensus_core_aa_length} amino acids before the stop, versus "
            f"{reference_core_aa_length} in {reference_id} ({abs(aa_length_difference)} {direction}); "
            f"only one to {PROTEIN_ORF_MAX_SMALL_INDEL_AA} missing or extra amino acids are accepted. "
            f"CDS not annotated. Protein alignment: {relative_alignment_path}"
        )
        return result

    ambiguous_offsets = [
        offset for offset, base in enumerate(cds_nt)
        if base not in "ACGT"
    ]
    if ambiguous_offsets:
        if depth_values is None:
            result["reason"] = (
                "Ambiguous bases occur within the inferred CDS, but no raw read-depth vector is available "
                "to distinguish covered low-depth Ns from true zero-coverage gaps; CDS not annotated. "
                f"Protein alignment: {relative_alignment_path}"
            )
            return result
        if len(depth_values) != len(consensus):
            result["reason"] = (
                f"Ambiguous bases occur within the inferred CDS, but consensus length {len(consensus)} "
                f"differs from mapping-depth coordinate length {len(depth_values)}; N positions cannot be "
                f"verified as covered, so the CDS was not annotated. Protein alignment: "
                f"{relative_alignment_path}"
            )
            return result
        zero_depth_positions = [
            start0 + offset + 1
            for offset in ambiguous_offsets
            if depth_values[start0 + offset] <= 0
        ]
        if zero_depth_positions:
            preview = ",".join(str(pos) for pos in zero_depth_positions[:10])
            if len(zero_depth_positions) > 10:
                preview += ",..."
            result["reason"] = (
                f"The inferred CDS contains {len(zero_depth_positions)} zero-read-depth gap position(s) "
                f"({preview}); this is an uncovered ORF gap rather than merely a low-depth N call, so "
                f"the CDS was not annotated. Protein alignment: {relative_alignment_path}"
            )
            return result

    start_codon = cds_nt[:3]
    stop_codon = cds_nt[-3:]
    if set(start_codon) <= set("ACGT") and start_codon != "ATG":
        result["reason"] = (
            f"The inferred CDS start is fully called as {start_codon}, not ATG; CDS not annotated. "
            f"Protein alignment: {relative_alignment_path}"
        )
        return result
    if set(stop_codon) <= set("ACGT") and stop_codon not in STOP_CODONS:
        result["reason"] = (
            f"The inferred terminal codon is fully called as {stop_codon}, not a stop codon; "
            f"CDS not annotated. Protein alignment: {relative_alignment_path}"
        )
        return result

    for nt_offset in range(3, len(cds_nt) - 3, 3):
        codon = cds_nt[nt_offset:nt_offset + 3]
        if set(codon) <= set("ACGT") and GENETIC_CODE.get(codon) == "*":
            codon_number = nt_offset // 3 + 1
            result["reason"] = (
                f"A fully called premature stop occurs at CDS codon {codon_number}; "
                f"CDS not annotated. Protein alignment: {relative_alignment_path}"
            )
            return result

    ambiguous_bases = len(ambiguous_offsets)
    ambiguous_codons = sum(
        1 for nt_offset in range(0, len(cds_nt), 3)
        if any(base not in "ACGT" for base in cds_nt[nt_offset:nt_offset + 3])
    )
    reason = (
        f"CDS annotated from positive-sense frame {frame + 1} by protein alignment to best "
        f"same-segment genotyping-database ORF {reference_id}; {matches}/{comparable} comparable "
        f"residues matched ({identity * 100:.1f}% identity; "
        f"{comparable_fraction * 100:.1f}% of reference residues comparable). "
    )

    small_difference_present = (
        internal_gap_columns > 0
        or aa_length_difference != 0
        or start0 != alignment_start0
        or stop_start0 != alignment_coding_end0
    )
    if small_difference_present:
        length_word = "extra" if aa_length_difference > 0 else "missing"
        if aa_length_difference == 0:
            length_statement = "the final protein lengths are equal"
        else:
            length_statement = (
                f"the consensus protein has {abs(aa_length_difference)} {length_word} amino acid(s) "
                "relative to the reference"
            )
        reason += (
            f"Small in-frame protein difference accepted: the alignment contains "
            f"{missing_amino_acids} amino acid(s) missing and {extra_amino_acids} extra in the "
            f"consensus relative to the reference; {length_statement}. Differences of one to "
            f"{PROTEIN_ORF_MAX_SMALL_INDEL_AA} amino acids are allowed. "
        )

    if ambiguous_bases:
        reason += (
            f"The complete {len(cds_nt)}-nt CDS coordinate span is present and was retained "
            f"despite {ambiguous_bases} ambiguous base(s) in {ambiguous_codons} codon(s); raw depth "
            "confirmed nonzero read coverage at every ambiguous CDS position, and coordinates remain "
            "protein-guided across the N-masked codons. "
        )
    else:
        reason += "No ambiguous bases occur within the annotated CDS. "
    reason += f"Protein alignment: {relative_alignment_path}."

    result.update({
        "annotate": True,
        "start": start0 + 1,
        "end": end0,
        "strand": "+",
        "cds_nt": cds_nt,
        "orf_source": "genotype_db_best_hit_protein_frame_alignment",
        "cds_complete": "YES",
        "reason": reason,
    })
    log(
        f"{segment}: CDS annotated at {start0 + 1}..{end0} on positive-sense frame {frame + 1} "
        f"using {reference_id} ({identity * 100:.1f}% protein identity; "
        f"missing_aa={missing_amino_acids}; extra_aa={extra_amino_acids})",
        log_path,
    )
    return result

def count_alignment_differences(aligned_a: str, aligned_b: str) -> dict[str, int | float]:
    matches = 0
    mismatches = 0
    indel_columns = 0
    indel_events = 0
    compared = 0
    previous_indel = False
    for a, b in zip(aligned_a, aligned_b):
        if a == "-" and b == "-":
            continue
        is_indel = (a == "-" or b == "-")
        if is_indel:
            indel_columns += 1
            if not previous_indel:
                indel_events += 1
            previous_indel = True
            continue
        previous_indel = False
        compared += 1
        if a.upper() == b.upper():
            matches += 1
        else:
            mismatches += 1
    identity_no_indels = 100.0 * matches / compared if compared else 0.0
    total_columns = matches + mismatches + indel_columns
    identity_with_indels = 100.0 * matches / total_columns if total_columns else 0.0
    return {
        "AlignmentLength": len(aligned_a),
        "Matches": matches,
        "Mismatches": mismatches,
        "IndelColumns": indel_columns,
        "IndelEvents": indel_events,
        "ComparedNonGapColumns": compared,
        "IdentityPctNoIndels": identity_no_indels,
        "IdentityPctIncludingIndels": identity_with_indels,
    }


def align_orfs_and_count(
    mapping_dir: Path,
    sample: str,
    segment: str,
    spades_orf: str,
    mapping_orf: str,
    log_path: Path,
    aligner: str,
) -> tuple[dict[str, int | float], str, str]:
    """Align SPAdes and mapping ORFs, return counts, and discard alignment files."""
    with tempfile.TemporaryDirectory(prefix=f"rotax_{safe_file_component(segment)}_orf_", dir=mapping_dir) as tmp_name:
        tmp_dir = Path(tmp_name)
        input_fasta = tmp_dir / "spades_vs_mapping_orf.fasta"
        aligned_fasta = tmp_dir / "spades_vs_mapping_orf.aln.fasta"
        write_fasta([
            (f">{sample}|{segment}|SPADES_ORF", normalize_dna(spades_orf)),
            (f">{sample}|{segment}|MAPPING_ORF", normalize_dna(mapping_orf)),
        ], input_fasta)

        if aligner.lower() != "none" and shutil.which(aligner):
            command = [aligner, "--auto", str(input_fasta)]
            log("Running ORF comparison alignment: " + " ".join(command), log_path)
            with aligned_fasta.open("w", encoding="utf-8") as out_handle:
                result = subprocess.run(command, stdout=out_handle)
            if result.returncode == 0 and aligned_fasta.exists():
                records = parse_fasta(aligned_fasta)
                if len(records) >= 2:
                    counts = count_alignment_differences(records[0][1], records[1][1])
                    return counts, "OK", "Aligned with MAFFT; temporary alignment removed"
            log(f"WARNING: ORF comparison alignment failed for {segment}; trying direct comparison if possible", log_path)

        if len(spades_orf) == len(mapping_orf):
            counts = count_alignment_differences(normalize_dna(spades_orf), normalize_dna(mapping_orf))
            return counts, "OK", "Same-length direct comparison; no external alignment used"

    empty_counts = {
        "AlignmentLength": 0,
        "Matches": 0,
        "Mismatches": 0,
        "IndelColumns": 0,
        "IndelEvents": 0,
        "ComparedNonGapColumns": 0,
        "IdentityPctNoIndels": 0.0,
        "IdentityPctIncludingIndels": 0.0,
    }
    return empty_counts, "NOT_COMPARED", "No aligner available and ORF lengths differ"


def cleanup_mapping_intermediates(mapping_dir: Path, sample_dir: Path, ref_fasta: Path, log_path: Path) -> None:
    """Remove mapping intermediates while retaining final VCF/AD evidence."""
    for path in [
        sample_dir / "rotavirus_mapping_orf_reference_query.fasta",
        sample_dir / "rotavirus_mapping_orf_reference_blast.tsv",
        sample_dir / "rotavirus_mapping_full_segment_query.fasta",
        sample_dir / "rotavirus_full_segment_blast_raw.tsv",
        mapping_dir / "rotavirus_mapping.sorted.bam",
        mapping_dir / "rotavirus_mapping.sorted.bam.bai",
        mapping_dir / "rotavirus_mapping.low_coverage.bed",
        mapping_dir / "rotavirus_mapping_consensus.raw.fasta",
        mapping_dir / "rotavirus_mapping_depth.raw.tsv",
    ]:
        safe_remove(path, log_path)
    for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]:
        safe_remove(Path(str(ref_fasta) + suffix), log_path)
    # Iterative-polish intermediates: polished mapping targets, their BWA indices,
    # .fai indices, and per-round variant calls.
    for tmp_path in mapping_dir.glob("rotavirus_mapping_references.polish*"):
        safe_remove(tmp_path, log_path)
    for tmp_path in mapping_dir.glob("rotavirus_mapping.polish*.vcf.gz*"):
        safe_remove(tmp_path, log_path)
    for tmp_path in mapping_dir.glob("rotax_*_orf_*"):
        safe_remove(tmp_path, log_path)
    for tmp_path in mapping_dir.glob("rotax_orf_coordinate_*"):
        safe_remove(tmp_path, log_path)
    log(
        f"{sample_dir.name}: retained final variant evidence "
        f"{mapping_dir / 'rotavirus_mapping.variants.vcf.gz'} and index",
        log_path,
    )


def write_descriptive_mapping_depth(
    raw_depth_path: Path,
    output_path: Path,
    sample: str,
    ref_by_name: dict[str, dict[str, object]],
    consensus_id_by_ref: dict[str, str],
    consensus_by_ref: dict[str, str],
    min_depth: int,
) -> None:
    """Replace internal BWA names with descriptive reference and consensus columns."""
    with output_path.open("w", encoding="utf-8") as out:
        out.write(
            "Sample\tSegment\tORFCalledGenotype\tMappingReferenceAccession\tMappingReferenceID\t"
            "ConsensusID\tReferencePosition\tDepth\tMeetsConsensusDepth\tConsensusBase\n"
        )
        if not raw_depth_path.exists():
            return
        with raw_depth_path.open("r", encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                parts = raw.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                ref_name, pos_text, depth_text = parts[:3]
                meta = ref_by_name.get(ref_name)
                if meta is None:
                    continue
                try:
                    pos = int(pos_text)
                    depth = int(depth_text)
                except ValueError:
                    continue
                consensus_seq = consensus_by_ref.get(ref_name, "")
                consensus_base = consensus_seq[pos - 1] if 1 <= pos <= len(consensus_seq) else ""
                out.write(
                    "\t".join([
                        sample,
                        str(meta.get("segment", "")),
                        str(meta.get("orf_called_genotype", "")),
                        str(meta.get("reference_accession", "")),
                        str(meta.get("reference_id_without_genotype", "")),
                        consensus_id_by_ref.get(ref_name, ""),
                        str(pos),
                        str(depth),
                        "YES" if depth >= min_depth else "NO",
                        consensus_base,
                    ]) + "\n"
                )


def run_mapping_consensus(sample_dir: Path, args: argparse.Namespace, log_path: Path) -> None:
    """Create read-derived full-segment consensuses from ORF-selected references.

    The SPAdes/ORF genotyping workflow is not changed. Each selected SPAdes ORF
    is BLASTed against the full-segment database, and the closest same-segment
    full reference is selected by bitscore first and percent identity second.
    References are oriented to positive sense, and reads are mapped once against
    all selected references for the sample. CDS coordinates are then determined
    from the best same-segment genotyping-DB ORF hit by protein-aligning all three
    forward translations of each positive-sense mapping consensus.
    """
    if not bool(getattr(args, "mapping_consensus", True)):
        return

    mapping_dir = sample_dir / "mapping"
    mapping_dir.mkdir(exist_ok=True)
    ref_fasta = mapping_dir / "rotavirus_mapping_references.fasta"

    full_db = str(getattr(args, "full_segment_db", "") or "")
    if not full_db:
        write_empty_mapping_outputs(sample_dir, "No full-segment database available for mapping-reference selection", log_path)
        return

    possible_mixed = sample_has_possible_mixed_infection(sample_dir)
    extra_info = read_sample_extra_information(sample_dir)
    allow_mixed_genbank = bool(getattr(args, "allow_mixed_genbank", False))
    if possible_mixed and not allow_mixed_genbank:
        log(
            f"{sample_dir.name}: possible mixed infection/genotype signal detected: {extra_info}. "
            "Mapping will be performed, but consensuses will remain review-only for GenBank by default.",
            log_path,
        )

    r1, r2, read_source = find_mapping_fastq_pair(sample_dir, args)
    if r1 is None or r2 is None:
        if getattr(args, "fasta_input", None):
            reason = (
                f"FASTA-input mode has no reads for mapping ({read_source}). "
                "Re-run with --mapping-r1 and --mapping-r2."
            )
        elif read_source == "trimmed_reads_missing_raw_fallback_disabled":
            reason = (
                "Trimmed/deduped reads (R1_paired.fastq.gz / deduped_R1.fastq.gz) "
                "were not found, so mapping was refused rather than silently using "
                "the raw, untrimmed FASTQs. This usually happens when a rerun with "
                "--skip-assembly --overwrite-existing deleted the trimmed reads "
                "without regenerating them. Re-run without --skip-assembly to "
                "regenerate trimmed reads, or pass --allow-raw-fastq-mapping to "
                "deliberately map the raw FASTQs."
            )
            log(f"{sample_dir.name}: {reason}", log_path)
        else:
            reason = f"No paired FASTQ files available for mapping ({read_source})"
        write_empty_mapping_outputs(sample_dir, reason, log_path)
        return

    query_records, query_map = build_mapping_query_records_and_map(sample_dir, possible_mixed=possible_mixed)
    if not query_records or not query_map:
        write_empty_mapping_outputs(sample_dir, "No selected ORF candidates available for mapping", log_path)
        return

    try:
        selected_hits_by_query = best_segment_blast_hits_for_mapping(sample_dir, args, query_records, query_map, log_path)

        ref_records: list[tuple[str, str]] = []
        ref_by_name: dict[str, dict[str, object]] = {}
        seen_ref_sequences: dict[tuple[str, str], str] = {}
        used_ref_names: set[str] = set()

        def qid_sort_key(qid: str) -> tuple[int, str, str]:
            candidate = query_map[qid]
            segment = str(candidate.get("segment", ""))
            rank = SEGMENT_ORDER.index(segment) if segment in SEGMENT_ORDER else len(SEGMENT_ORDER)
            return (rank, normalize_genotype_for_compare(str(candidate.get("genotype", ""))), qid)

        for qid in sorted(query_map, key=qid_sort_key):
            candidate = query_map[qid]
            segment = str(candidate.get("segment", ""))
            genotype = normalize_genotype_for_compare(str(candidate.get("genotype", "")))
            hits = selected_hits_by_query.get(qid, [])
            if not hits:
                log(f"{sample_dir.name} {segment} {genotype}: no same-segment full-reference hit found", log_path)
                continue

            hit = hits[0]
            fetched = fetch_blastdb_sequence(full_db, str(hit["sseqid"]), log_path)
            if fetched is None:
                log(f"{sample_dir.name} {segment} {genotype}: selected reference could not be fetched", log_path)
                continue
            _ref_header, ref_seq = fetched
            clean_ref = orient_reference_like_hit(ref_seq, hit)
            if str(hit.get("sstrand", "+")).strip().lower() == "minus":
                log(
                    f"{sample_dir.name} {segment} {genotype}: reverse-complemented selected full-segment "
                    "reference before mapping so the mapping consensus is positive sense",
                    log_path,
                )
            if not clean_ref:
                log(f"{sample_dir.name} {segment} {genotype}: selected reference sequence was empty", log_path)
                continue

            dedup_key = (segment, clean_ref)
            if dedup_key in seen_ref_sequences:
                existing_name = seen_ref_sequences[dedup_key]
                aliases = ref_by_name[existing_name].setdefault("candidate_aliases", [])
                if isinstance(aliases, list):
                    aliases.append(qid)
                log(
                    f"{sample_dir.name} {segment} {genotype}: selected the same reference sequence as "
                    f"{existing_name}; one shared mapping consensus will be produced",
                    log_path,
                )
                continue

            reference_id = reference_id_without_leading_genotype(str(hit["sseqid"]))
            accession = mapping_reference_accession(str(hit["sseqid"]))
            ref_name = unique_mapping_reference_name(segment, accession, used_ref_names)
            seen_ref_sequences[dedup_key] = ref_name
            ref_records.append((
                f">{ref_name} segment={segment} source_reference={reference_id}",
                clean_ref,
            ))
            ref_by_name[ref_name] = {
                "candidate": candidate,
                "candidate_key": qid,
                "candidate_aliases": [],
                "segment": segment,
                "orf_called_genotype": genotype,
                "reference_id_without_genotype": reference_id,
                "reference_accession": accession,
                "length": len(clean_ref),
                "full_reference_pident": float(hit["pident"]),
                "full_reference_bitscore": float(hit["bitscore"]),
            }

        if not ref_records:
            write_empty_mapping_outputs(sample_dir, "No usable mapping references selected", log_path)
            return

        write_fasta(ref_records, ref_fasta)
        reference_lengths = {fasta_key(h): len(normalize_dna(seq)) for h, seq in ref_records}
        log(
            f"{sample_dir.name}: mapping {len(ref_records)} descriptive full-reference sequence(s) "
            f"selected from {len(query_map)} ORF candidate(s); possible_mixed={possible_mixed}",
            log_path,
        )

        bam = mapping_dir / "rotavirus_mapping.sorted.bam"
        vcf = mapping_dir / "rotavirus_mapping.variants.vcf.gz"
        raw_depth_path = mapping_dir / "rotavirus_mapping_depth.raw.tsv"
        final_depth_path = mapping_dir / "rotavirus_mapping_depth.tsv"
        lowcov_bed = mapping_dir / "rotavirus_mapping.low_coverage.bed"
        raw_consensus = mapping_dir / "rotavirus_mapping_consensus.raw.fasta"
        consensus_fasta = mapping_dir / "rotavirus_mapping_consensus.fasta"

        min_depth = int(getattr(args, "mapping_min_depth", 5))
        polish_rounds = max(0, int(getattr(args, "mapping_polish_rounds", 1)))

        def bwa_map_sort_index(target_fasta: Path) -> bool:
            """(Re)map the sample reads to target_fasta into the shared sorted BAM."""
            ok_local = run_command(["bwa", "index", target_fasta.name], log_path, cwd=mapping_dir)
            cmd = (
                f"bwa mem -B {MAPPING_BWA_MISMATCH_PENALTY} -t {int(args.threads)} "
                f"{shlex.quote(target_fasta.name)} "
                f"{shlex.quote(str(r1))} {shlex.quote(str(r2))} "
                f"| samtools sort -@ {int(args.threads)} -o {shlex.quote(bam.name)} -"
            )
            ok_local = ok_local and run_shell_command(cmd, log_path, cwd=mapping_dir)
            ok_local = ok_local and run_command(["samtools", "index", bam.name], log_path, cwd=mapping_dir)
            return ok_local

        # Round 0: map reads against the ORF-selected GenBank reference.
        current_target = ref_fasta
        ok = bwa_map_sort_index(current_target)
        if not ok:
            write_empty_mapping_outputs(sample_dir, "Read mapping failed", log_path)
            return

        # Iterative remap-to-consensus polish. Each round rebuilds a SNP-corrected
        # mapping target that stays in the ORIGINAL reference's coordinate frame:
        # substitutions from the reads are applied, indels are excluded, and
        # low-depth positions are kept as the reference base (not masked to N).
        # Remapping the same reads to this sample-specific target lets read ends
        # that soft-clipped against the divergent GenBank reference realign, which
        # extends terminal coverage. Because the target length is unchanged, all
        # downstream coverage / missing-end metrics remain in reference coordinates.
        for polish_i in range(1, polish_rounds + 1):
            polish_vcf = mapping_dir / f"rotavirus_mapping.polish{polish_i}.vcf.gz"
            polished_ref = mapping_dir / f"rotavirus_mapping_references.polish{polish_i}.fasta"
            # Call variants and keep SNP records only (-v snps). Filtering to SNPs
            # here guarantees the polished target stays in reference coordinates:
            # indels would shift downstream offsets, and `bcftools consensus -e
            # TYPE="indel"` is known to miss mixed multiallelic SNP+indel records,
            # so we drop indels at the call-set level instead.
            call_cmd = (
                f"bcftools mpileup -f {shlex.quote(current_target.name)} -Ou {shlex.quote(bam.name)} "
                f"| bcftools call --ploidy 1 -mv -Ou "
                f"| bcftools view -v snps -Oz -o {shlex.quote(polish_vcf.name)}"
            )
            ok = run_shell_command(call_cmd, log_path, cwd=mapping_dir)
            ok = ok and run_command(["bcftools", "index", "-f", polish_vcf.name], log_path, cwd=mapping_dir)
            # SNP-only, unmasked consensus => coordinate-stable polished target.
            polish_cons_cmd = (
                f"bcftools consensus -f {shlex.quote(current_target.name)} {shlex.quote(polish_vcf.name)}"
            )
            ok = ok and run_command_stdout(["bash", "-lc", polish_cons_cmd], polished_ref, log_path, cwd=mapping_dir)
            if not ok or not is_nonempty(polished_ref):
                log(f"{sample_dir.name}: mapping polish round {polish_i} failed; keeping previous target", log_path)
                bwa_map_sort_index(current_target)
                break
            if not bwa_map_sort_index(polished_ref):
                log(f"{sample_dir.name}: remap in polish round {polish_i} failed; keeping previous target", log_path)
                bwa_map_sort_index(current_target)
                break
            current_target = polished_ref
            log(f"{sample_dir.name}: completed mapping polish round {polish_i}", log_path)

        # Final depth + variant call + masked consensus against the final target.
        ok = run_command_stdout(
            ["samtools", "depth", "-J", "-aa", bam.name],
            raw_depth_path,
            log_path,
            cwd=mapping_dir,
        )
        if not ok:
            write_empty_mapping_outputs(sample_dir, "Read mapping failed", log_path)
            return

        depth_stats = parse_depth_file(raw_depth_path, reference_lengths, min_depth)
        depth_vectors = read_depth_vectors(raw_depth_path, reference_lengths)
        write_low_coverage_bed(raw_depth_path, reference_lengths, min_depth, lowcov_bed)

        # Keep an IUPAC genotype only if the less-supported called allele has
        # at least five high-quality AD reads and at least 20% of AD across the
        # two called alleles. Otherwise +setGT -n X chooses the AD-major allele.
        failing_heterozygote = (
            'GT="het" && '
            f"(MIN(FMT/AD[GT])<{IUPAC_MIN_MINOR_AD} || "
            "MIN(FMT/AD[GT])/SUM(FMT/AD[GT])"
            f"<{IUPAC_MIN_MINOR_FRACTION:g})"
        )
        mpileup_cmd = (
            f"bcftools mpileup -f {shlex.quote(current_target.name)} -Ou {shlex.quote(bam.name)} "
            f"| bcftools call --ploidy 2 -mv -Ou "
            f"| bcftools +setGT -Ou -- -t q -n X -i {shlex.quote(failing_heterozygote)} "
            f"| bcftools +fill-tags -Ou -- -t AC,AN "
            f"| bcftools view -Oz -o {shlex.quote(vcf.name)}"
        )
        ok = run_shell_command(mpileup_cmd, log_path, cwd=mapping_dir)
        ok = ok and run_command(["bcftools", "index", "-f", vcf.name], log_path, cwd=mapping_dir)
        consensus_cmd = (
            f"bcftools consensus -f {shlex.quote(current_target.name)} --mask {shlex.quote(lowcov_bed.name)} "
            f"--mask-with N --iupac-codes {shlex.quote(vcf.name)}"
        )
        ok = ok and run_command_stdout(["bash", "-lc", consensus_cmd], raw_consensus, log_path, cwd=mapping_dir)
        if not ok or not raw_consensus.exists():
            write_empty_mapping_outputs(sample_dir, "Consensus generation failed", log_path)
            return

        raw_records = parse_fasta(raw_consensus)
        consensus_records: list[tuple[str, str]] = []
        full_consensus_records: list[tuple[str, str]] = []
        review_consensus_records: list[tuple[str, str]] = []
        review_consensus_ids: set[str] = set()
        consensus_by_ref: dict[str, str] = {}
        consensus_id_by_ref: dict[str, str] = {}
        status_by_ref: dict[str, tuple[str, str, list[str]]] = {}
        ready_by_ref: dict[str, tuple[str, str]] = {}
        used_consensus_ids: set[str] = set()

        for header, seq in raw_records:
            ref_name = fasta_key(header)
            if ref_name not in ref_by_name:
                continue
            meta = ref_by_name[ref_name]
            segment = str(meta.get("segment", ""))
            genotype = str(meta.get("orf_called_genotype", ""))
            accession = str(meta.get("reference_accession", ""))
            status, full_or_partial, reasons = consensus_record_status(depth_stats.get(ref_name, {}), args)
            genbank_ready = "YES"
            ready_reason = "Full mapping consensus and no mixed-infection restriction"
            if full_or_partial != "full":
                genbank_ready = "NO"
                ready_reason = "Mapping consensus is partial"
            elif possible_mixed and not allow_mixed_genbank:
                genbank_ready = "NO"
                ready_reason = "Possible mixed infection; conservative GenBank readiness disabled"
                status = "FULL_MAPPING_CONSENSUS_MIXED_REVIEW"
                reasons = list(reasons) + ["Possible mixed infection; full sequence retained for review"]
            elif possible_mixed and allow_mixed_genbank:
                ready_reason = "Mixed-sample override enabled with --allow-mixed-genbank; CDS annotation still required"

            base_id = f"{sample_dir.name}|{segment}|{genotype}|mapping_consensus|ref={accession}"
            consensus_id = base_id
            alt_number = 2
            while consensus_id in used_consensus_ids:
                consensus_id = f"{base_id}|alt={alt_number}"
                alt_number += 1
            used_consensus_ids.add(consensus_id)

            clean_seq = normalize_dna(seq)
            fasta_header = f">{consensus_id} status={status} orientation=positive_sense"
            consensus_records.append((fasta_header, clean_seq))
            if full_or_partial == "full":
                full_consensus_records.append((fasta_header, clean_seq))
            if full_or_partial != "full" or genbank_ready != "YES":
                review_consensus_records.append((fasta_header, clean_seq))
                review_consensus_ids.add(consensus_id)
            consensus_by_ref[ref_name] = clean_seq
            consensus_id_by_ref[ref_name] = consensus_id
            status_by_ref[ref_name] = (status, full_or_partial, reasons)
            ready_by_ref[ref_name] = (genbank_ready, ready_reason)

        write_fasta(consensus_records, consensus_fasta)
        write_fasta(full_consensus_records, mapping_dir / "rotavirus_mapping_consensus_full.fasta")

        # Keep the selected SPAdes ORF only for the diagnostic SPAdes-versus-
        # mapping-ORF comparison. It no longer determines CDS coordinates.
        spades_orf_by_ref: dict[str, str] = {}
        for ref_name, meta in ref_by_name.items():
            candidate = meta.get("candidate", {})
            if isinstance(candidate, dict):
                orf = extract_spades_orf_for_candidate(candidate)
                if orf:
                    spades_orf_by_ref[ref_name] = normalize_dna(orf)

        # BLAST each positive-sense mapping consensus against the validated ORF
        # genotyping database. The selected best same-segment hit is fetched and
        # protein-aligned to all three forward consensus frames below.
        segment_by_consensus_id = {
            consensus_id_by_ref[rn]: str(ref_by_name[rn].get("segment", ""))
            for rn in consensus_id_by_ref
            if rn in ref_by_name
        }
        genotype_orf_hits = run_genotype_orf_blast(
            mapping_dir,
            consensus_fasta,
            str(getattr(args, "blast_db", "") or ""),
            segment_by_consensus_id,
            log_path,
        )

        mapping_orf_records: list[tuple[str, str]] = []
        summary_rows: list[dict[str, str | int | float]] = []
        feature_rows: list[dict[str, str | int | float]] = []

        def ref_sort_key(ref_name: str) -> tuple[int, str, str]:
            segment = str(ref_by_name[ref_name].get("segment", ""))
            rank = SEGMENT_ORDER.index(segment) if segment in SEGMENT_ORDER else len(SEGMENT_ORDER)
            return (rank, str(ref_by_name[ref_name].get("orf_called_genotype", "")), ref_name)

        for ref_name in sorted(ref_by_name, key=ref_sort_key):
            consensus_seq = consensus_by_ref.get(ref_name, "")
            if not consensus_seq:
                continue
            meta = ref_by_name[ref_name]
            candidate = meta.get("candidate", {}) if isinstance(meta.get("candidate", {}), dict) else {}
            segment = str(meta.get("segment", ""))
            genotype = str(meta.get("orf_called_genotype", ""))
            consensus_id = consensus_id_by_ref.get(ref_name, "")
            status, full_or_partial, reasons = status_by_ref.get(
                ref_name,
                ("PARTIAL_MAPPING_CONSENSUS", "partial", ["No depth statistics found"]),
            )
            genbank_ready, ready_reason = ready_by_ref.get(
                ref_name,
                ("NO", "No GenBank-readiness decision found"),
            )
            stats = depth_stats.get(ref_name, {})
            expected_orf_len = EXPECTED_ORF_LENGTHS.get(segment, 0)
            n_bases = sum(1 for ch in consensus_seq if ch.upper() not in "ACGT")

            gt_hit = genotype_orf_hits.get(consensus_id)
            cds_annotation = annotate_consensus_cds_by_protein(
                consensus_seq,
                gt_hit,
                str(getattr(args, "blast_db", "") or ""),
                segment,
                consensus_id,
                mapping_dir,
                str(getattr(args, "mapping_orf_alignment_cmd", "mafft") or "mafft"),
                log_path,
                depth_values=depth_vectors.get(ref_name),
            )
            orf_start = cds_annotation["start"]
            orf_end = cds_annotation["end"]
            orf_strand = str(cds_annotation["strand"])
            orf_source = str(cds_annotation["orf_source"])
            mapping_orf = str(cds_annotation["cds_nt"])
            mapping_orf_len: int | str = len(mapping_orf) if mapping_orf else ""
            orf_diff: int | str = (
                int(mapping_orf_len) - expected_orf_len
                if mapping_orf and expected_orf_len
                else ""
            )
            cds_annotated = bool(cds_annotation["annotate"])
            cds_complete = str(cds_annotation["cds_complete"])
            cds_reason = str(cds_annotation["reason"])
            terminals_resolved, terminal_resolution_reason = (
                genbank_terminal_regions_resolved(
                    consensus_seq,
                    segment,
                    orf_start,
                    orf_end,
                )
            )

            # Final GenBank readiness requires both the earlier segment-level
            # checks and a successfully annotated CDS. A full consensus without
            # a CDS or with an N in either terminal region is retained for review
            # but is never marked GenBankReady.
            if not cds_annotated:
                if genbank_ready == "YES":
                    genbank_ready = "NO"
                    ready_reason = f"CDS was not annotated: {cds_reason}"
                else:
                    ready_reason = f"{ready_reason}; CDS was not annotated: {cds_reason}"
            elif not terminals_resolved:
                if genbank_ready == "YES":
                    genbank_ready = "NO"
                    ready_reason = terminal_resolution_reason
                else:
                    ready_reason = (
                        f"{ready_reason}; {terminal_resolution_reason}"
                    )
            elif genbank_ready == "YES":
                ready_reason = (
                    f"{ready_reason}; CDS annotated; {terminal_resolution_reason}"
                )

            if genbank_ready != "YES" and consensus_id not in review_consensus_ids:
                review_header = f">{consensus_id} status={status} orientation=positive_sense"
                review_consensus_records.append((review_header, consensus_seq))
                review_consensus_ids.add(consensus_id)

            if mapping_orf:
                mapping_orf_records.append((
                    f">{sample_dir.name}|{segment}|{genotype}|mapping_ORF|ref={meta.get('reference_accession', '')} "
                    f"coords={orf_start}..{orf_end} strand={orf_strand} orf_source={orf_source}",
                    mapping_orf,
                ))

            spades_orf = spades_orf_by_ref.get(ref_name, "")
            if spades_orf and mapping_orf:
                counts, comparison_status, comparison_reason = align_orfs_and_count(
                    mapping_dir,
                    sample_dir.name,
                    segment,
                    spades_orf,
                    mapping_orf,
                    log_path,
                    str(getattr(args, "mapping_orf_alignment_cmd", "mafft") or "mafft"),
                )
            else:
                counts = {
                    "AlignmentLength": 0,
                    "Matches": 0,
                    "Mismatches": 0,
                    "IndelColumns": 0,
                    "IndelEvents": 0,
                    "ComparedNonGapColumns": 0,
                    "IdentityPctNoIndels": 0.0,
                    "IdentityPctIncludingIndels": 0.0,
                }
                comparison_status = "NOT_COMPARED"
                comparison_reason = "Missing SPAdes ORF or mapping ORF sequence"

            difference_columns = int(counts["Mismatches"]) + int(counts["IndelColumns"])

            # CDS completeness and its explanatory note come directly from
            # the protein-guided annotation. N-masked codons are allowed when
            # the full reference-defined CDS span remains present; broken or
            # insufficiently supported ORFs are left unannotated.

            # Mapping references are forced to positive sense before read
            # mapping, so the consensus itself can be evaluated directly against
            # the conserved positive-sense RVA terminal motifs.
            terminal_motif_eval = evaluate_terminal_motifs(
                consensus_seq,
                segment,
                int(getattr(args, "genbank_terminal_motif_mismatches", 1)),
            )

            summary_rows.append({
                "Sample": sample_dir.name,
                "PossibleMixedInfection": "YES" if possible_mixed else "NO",
                "ExtraInformation": extra_info,
                "Segment": segment,
                "ORFCalledGenotype": genotype,
                "ConsensusID": consensus_id,
                "ConsensusStatus": status,
                "FullConsensusCalled": "YES" if full_or_partial == "full" else "NO",
                "GenBankReady": genbank_ready,
                "GenBankReadyReason": ready_reason,
                "StatusReason": "; ".join(reasons),
                "MappingReferenceAccession": meta.get("reference_accession", ""),
                "MappingReferenceID": meta.get("reference_id_without_genotype", ""),
                "MappingReferenceLength": meta.get("length", ""),
                "MappingReferencePident": meta.get("full_reference_pident", ""),
                "MappingReferenceBitscore": meta.get("full_reference_bitscore", ""),
                "ConsensusLength": len(consensus_seq),
                "CoveredBases": stats.get("covered_bases", ""),
                "CoveragePct": stats.get("coverage_pct", ""),
                "MeanDepth": stats.get("mean_depth", ""),
                "MedianDepth": stats.get("median_depth", ""),
                "MinDepthForConsensus": min_depth,
                "FirstCoveredPosition": stats.get("first_covered", ""),
                "LastCoveredPosition": stats.get("last_covered", ""),
                "MissingReference5primeBases": stats.get("missing5", ""),
                "MissingReference3primeBases": stats.get("missing3", ""),
                "ConsensusNBases": n_bases,
                "MappingTerminal5Status": terminal_motif_eval["terminal_5_status"],
                "MappingTerminal5Observed": terminal_motif_eval["terminal_5_observed"],
                "MappingTerminal5Mismatches": terminal_motif_eval["terminal_5_mismatches"],
                "MappingTerminal3Status": terminal_motif_eval["terminal_3_status"],
                "MappingTerminal3Observed": terminal_motif_eval["terminal_3_observed"],
                "MappingTerminal3Mismatches": terminal_motif_eval["terminal_3_mismatches"],
                "MappingTerminalMotifsPass": "YES" if terminal_motif_eval["terminal_motifs_pass"] else "NO",
                "MappingTerminalMotifClassification": terminal_motif_eval["terminal_motif_classification"],
                "CDSComplete": cds_complete,
                "CDSCompletenessReason": cds_reason,
                "ProteinAlignmentFile": cds_annotation.get("protein_alignment_file", ""),
                "ORFSource": orf_source,
                "ORFStart": orf_start,
                "ORFEnd": orf_end,
                "ORFStrand": orf_strand,
                "MappingORFLength": mapping_orf_len,
                "ExpectedORFLength": expected_orf_len,
                "ORFLengthDifference": orf_diff,
                "SpadesORFLength": len(spades_orf),
                "SpadesVsMappingAlignmentLength": counts["AlignmentLength"],
                "SpadesVsMappingMatches": counts["Matches"],
                "SpadesVsMappingMismatches": counts["Mismatches"],
                "SpadesVsMappingIndelColumns": counts["IndelColumns"],
                "SpadesVsMappingIndelEvents": counts["IndelEvents"],
                "SpadesVsMappingDifferenceColumns": difference_columns,
                "SpadesVsMappingIdentityPctIncludingIndels": counts["IdentityPctIncludingIndels"],
                "ORFComparisonStatus": comparison_status,
                "ORFComparisonReason": comparison_reason,
                "OriginalSpadesContigID": candidate.get("qseqid", "") if isinstance(candidate, dict) else "",
                "CandidateKey": meta.get("candidate_key", ""),
                "ReadSource": read_source,
            })

            common_feature = {
                "Sample": sample_dir.name,
                "Segment": segment,
                "ORFCalledGenotype": genotype,
                "ConsensusID": consensus_id,
                "ConsensusStatus": status,
                "FullConsensusCalled": "YES" if full_or_partial == "full" else "NO",
                "GenBankReady": genbank_ready,
                "MappingReferenceAccession": meta.get("reference_accession", ""),
                "MappingReferenceID": meta.get("reference_id_without_genotype", ""),
            }
            feature_rows.append({
                **common_feature,
                "Feature": "source",
                "Start": 1,
                "End": len(consensus_seq),
                "Strand": "+",
                "Gene": "",
                "Product": "rotavirus A genome segment",
                "Note": (
                    f"{full_or_partial} positive-sense read-mapping consensus; selected reference="
                    f"{meta.get('reference_accession', '')}; "
                    f"CDS {'annotated' if cds_annotated else 'NOT annotated'}: {cds_reason}"
                ),
            })
            if cds_annotated and orf_start and orf_end:
                feature_rows.append({
                    **common_feature,
                    "Feature": "CDS",
                    "Start": orf_start,
                    "End": orf_end,
                    "Strand": "+",
                    "Gene": segment,
                    "Product": segment,
                    "Note": cds_reason,
                })

        genbank_record_count, genbank_feature_count = write_mapping_genbank_exports(
            mapping_dir,
            sample_dir.name,
            consensus_records,
            summary_rows,
            feature_rows,
        )
        write_fasta(mapping_orf_records, mapping_dir / "rotavirus_mapping_orfs.fasta")
        write_fasta(
            review_consensus_records,
            mapping_dir / "rotavirus_mapping_consensus_partial_or_review.fasta",
        )
        with (mapping_dir / "rotavirus_mapping_summary.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=mapping_summary_fields())
            writer.writeheader()
            writer.writerows(summary_rows)
        with (mapping_dir / "rotavirus_mapping_consensus_features.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=mapping_consensus_feature_fields())
            writer.writeheader()
            writer.writerows(feature_rows)
        write_descriptive_mapping_depth(
            raw_depth_path,
            final_depth_path,
            sample_dir.name,
            ref_by_name,
            consensus_id_by_ref,
            consensus_by_ref,
            min_depth,
        )
        log(
            f"Generated simplified mapping outputs in {mapping_dir}: {len(summary_rows)} consensus summary row(s), "
            f"{len(full_consensus_records)} full consensus sequence(s), and "
            f"{genbank_record_count} terminal-N-trimmed per-sample GenBank "
            f"sequence(s) with {genbank_feature_count} adjusted feature row(s)",
            log_path,
        )
    finally:
        if bool(getattr(args, "cleanup", True)):
            cleanup_mapping_intermediates(mapping_dir, sample_dir, ref_fasta, log_path)
        else:
            log(
                f"{sample_dir.name}: --no-cleanup enabled; retaining mapping "
                f"intermediates in {mapping_dir}",
                log_path,
            )


# ---------------------------------------------------------------------------
# Integrated mapping-reference validation HTML
# ---------------------------------------------------------------------------

MAPPING_VALIDATION_IUPAC: dict[str, tuple[str, ...]] = {
    "A": ("A",), "C": ("C",), "G": ("G",), "T": ("T",),
    "R": ("A", "G"), "Y": ("C", "T"), "S": ("C", "G"), "W": ("A", "T"),
    "K": ("G", "T"), "M": ("A", "C"),
    "B": ("C", "G", "T"), "D": ("A", "G", "T"),
    "H": ("A", "C", "T"), "V": ("A", "C", "G"),
    # N is unknown/masked and must not be represented as supported ambiguity.
    "N": (),
}
MAPPING_VALIDATION_COLORS = {
    "A": "#2b9348", "C": "#277da1", "G": "#f8961e", "T": "#d62828",
}
MAPPING_VALIDATION_COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


def mapping_validation_header_fields(header: str) -> dict[str, str]:
    key = fasta_key(header)
    parts = key.split("|")
    fields = {
        "id": key,
        "sample": parts[0] if parts else "",
        "segment": parts[1].upper() if len(parts) > 1 else "",
        "genotype": parts[2] if len(parts) > 2 else "",
        "reference": "",
        "status": "",
    }
    for match in re.finditer(r"(\w+)=([^\s|]+)", header):
        name = match.group(1).lower()
        if name == "ref":
            fields["reference"] = match.group(2)
        elif name == "status":
            fields["status"] = match.group(2)
    return fields


def mapping_validation_load_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def mapping_validation_load_depth(
    path: Path,
) -> dict[str, dict[int, dict[str, object]]]:
    tracks: dict[str, dict[int, dict[str, object]]] = {}
    if not path.exists():
        return tracks
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            consensus_id = (row.get("ConsensusID", "") or "").strip()
            if not consensus_id:
                continue
            try:
                position = int(row.get("ReferencePosition", ""))
                depth = max(0, int(float(row.get("Depth", "0") or 0)))
            except (TypeError, ValueError):
                continue
            base = normalize_dna(row.get("ConsensusBase", "") or "N")[:1] or "N"
            tracks.setdefault(consensus_id, {})[position] = {
                "depth": depth,
                "base": base,
                "meets": (row.get("MeetsConsensusDepth", "") or "").upper() == "YES",
            }
    return tracks


def mapping_validation_depth_sequence(entries: dict[int, dict[str, object]]) -> str:
    if not entries:
        return ""
    return "".join(
        str(entries.get(pos, {}).get("base", "N"))
        for pos in range(1, max(entries) + 1)
    )


def mapping_validation_reference_records(
    path: Path,
) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    for header, sequence in parse_fasta(path):
        name = fasta_key(header)
        values = {
            match.group(1): match.group(2)
            for match in re.finditer(r"(\w+)=([^\s]+)", header)
        }
        segment = values.get("segment", "")
        if not segment and "_reference_" in name:
            segment = name.split("_reference_", 1)[0]
        source_reference = values.get("source_reference", "")
        records.append({
            "name": name,
            "segment": segment.upper(),
            "source_reference": source_reference,
            "accession": mapping_reference_accession(source_reference or name),
            "sequence": normalize_dna(sequence),
        })
    return records


def mapping_validation_choose_reference(
    references: list[dict[str, str]],
    segment: str,
    accession: str,
    reference_id: str,
) -> dict[str, str] | None:
    same_segment = [record for record in references if record["segment"] == segment]
    for record in same_segment:
        if reference_id and record["source_reference"] == reference_id:
            return record
    wanted = accession.split(".")[0].upper()
    for record in same_segment:
        if record["accession"].split(".")[0].upper() == wanted:
            return record
    return same_segment[0] if len(same_segment) == 1 else None


def mapping_validation_load_allele_depths(
    path: Path,
) -> dict[str, dict[int, dict[str, int]]]:
    by_reference: dict[str, dict[int, dict[str, int]]] = {}
    if not path.exists():
        return by_reference
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, position_text, _variant_id, ref, alt = fields[:5]
            alleles = [ref.upper()] + [value.upper() for value in alt.split(",")]
            if any(len(allele) != 1 or allele not in "ACGT" for allele in alleles):
                continue
            format_fields = fields[8].split(":")
            sample_fields = fields[9].split(":")
            if "AD" not in format_fields:
                continue
            ad_index = format_fields.index("AD")
            if ad_index >= len(sample_fields):
                continue
            try:
                position = int(position_text)
                depths = [max(0, int(value)) for value in sample_fields[ad_index].split(",")]
            except ValueError:
                continue
            if len(depths) < len(alleles):
                continue
            by_reference.setdefault(chrom, {})[position] = {
                allele: count for allele, count in zip(alleles, depths) if count > 0
            }
    return by_reference


def mapping_validation_longest_cds(record) -> dict[str, object] | None:
    best: dict[str, object] | None = None
    best_length = -1
    for feature in record.features:
        if feature.type != "CDS":
            continue
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        length = end - start + 1
        if length <= best_length:
            continue
        best_length = length
        product_values = feature.qualifiers.get("product", [])
        gene_values = feature.qualifiers.get("gene", [])
        best = {
            "start": start,
            "end": end,
            "strand": "+" if (feature.location.strand or 1) >= 0 else "-",
            "gene": gene_values[0] if gene_values else "",
            "product": product_values[0] if product_values else "",
        }
    return best


def mapping_validation_fetch_genbank(
    accession: str,
    email: str,
    api_key: str,
    cache: dict[str, tuple[str, dict[str, object] | None] | None],
    notes: list[str],
) -> tuple[str, dict[str, object] | None] | None:
    cache_key = accession.split(".")[0].upper()
    if cache_key in cache:
        return cache[cache_key]
    if not email:
        cache[cache_key] = None
        return None
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        notes.append("Biopython is unavailable; reference CDS annotations could not be fetched.")
        cache[cache_key] = None
        return None
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    last_error = ""
    for attempt in range(1, 4):
        try:
            handle = Entrez.efetch(
                db="nucleotide", id=accession, rettype="gb", retmode="text"
            )
            text = handle.read()
            handle.close()
            record = SeqIO.read(io.StringIO(text), "genbank")
            result = (
                normalize_dna(str(record.seq)),
                mapping_validation_longest_cds(record),
            )
            cache[cache_key] = result
            time.sleep(0.11 if api_key else 0.34)
            return result
        except Exception as exc:  # noqa: BLE001
            last_error = str(exc)
            time.sleep(attempt)
    notes.append(f"Could not fetch GenBank {accession}: {last_error}")
    cache[cache_key] = None
    return None


def mapping_validation_aligned_coord_map(sequence: str) -> list[int]:
    coordinates: list[int] = []
    position = 0
    for base in sequence:
        if base == "-":
            coordinates.append(0)
        else:
            position += 1
            coordinates.append(position)
    return coordinates


def mapping_validation_run_alignment(
    consensus: str,
    reference: str,
    depth_sequence: str,
    log_path: Path,
) -> tuple[str, str, str, bool, bool] | None:
    if shutil.which("mafft") is None:
        log("Mapping validation HTML skipped: mafft was not found in PATH", log_path)
        return None
    with tempfile.TemporaryDirectory(prefix="rotax_mapping_validation_") as temp_name:
        temp_dir = Path(temp_name)
        input_fasta = temp_dir / "input.fasta"
        records = [
            (">consensus_track", consensus),
            (">reference_track", reference),
        ]
        if depth_sequence:
            records.append((">depth_coordinate_track", depth_sequence))
        write_fasta(records, input_fasta)
        result = subprocess.run(
            ["mafft", "--auto", "--adjustdirection", str(input_fasta)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            log(
                "Mapping validation MAFFT failed: " + result.stderr.strip()[-500:],
                log_path,
            )
            return None
        parsed = parse_fasta_text(result.stdout)
    aligned: dict[str, str] = {}
    reversed_records: set[str] = set()
    for header, sequence in parsed:
        key = fasta_key(header)
        if key.startswith("_R_"):
            key = key[3:]
            reversed_records.add(key)
        # Keep MAFFT's gap columns.  The pipeline-wide normalize_dna() helper
        # intentionally removes "-" for unaligned input sequences.
        aligned[key] = re.sub(r"\s+", "", sequence).upper().replace("U", "T")
    if "consensus_track" not in aligned or "reference_track" not in aligned:
        return None
    alignment_length = max(len(value) for value in aligned.values())
    return (
        aligned["consensus_track"].ljust(alignment_length, "-"),
        aligned["reference_track"].ljust(alignment_length, "-"),
        aligned.get("depth_coordinate_track", "").ljust(alignment_length, "-"),
        "reference_track" in reversed_records,
        "depth_coordinate_track" in reversed_records,
    )


def mapping_validation_run_orf_pair_alignment(
    mapping_cds: str,
    spades_orf: str,
    log_path: Path,
) -> tuple[str, str] | None:
    """Align the mapping CDS and selected SPAdes nucleotide ORF for HTML review."""
    if not mapping_cds or not spades_orf:
        return None
    if shutil.which("mafft") is None:
        log(
            "SPAdes-versus-mapping HTML frame skipped: mafft was not found in PATH",
            log_path,
        )
        return None
    with tempfile.TemporaryDirectory(prefix="rotax_spades_mapping_html_") as temp_name:
        input_fasta = Path(temp_name) / "orf_pair.fasta"
        write_fasta(
            [
                (">mapping_cds_track", normalize_dna(mapping_cds)),
                (">spades_orf_track", normalize_dna(spades_orf)),
            ],
            input_fasta,
        )
        result = subprocess.run(
            ["mafft", "--auto", str(input_fasta)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            log(
                "SPAdes-versus-mapping validation MAFFT failed: "
                + result.stderr.strip()[-500:],
                log_path,
            )
            return None
        parsed = {
            fasta_key(header): re.sub(r"\s+", "", sequence).upper().replace("U", "T")
            for header, sequence in parse_fasta_text(result.stdout)
        }
    if "mapping_cds_track" not in parsed or "spades_orf_track" not in parsed:
        return None
    alignment_length = max(
        len(parsed["mapping_cds_track"]),
        len(parsed["spades_orf_track"]),
    )
    return (
        parsed["mapping_cds_track"].ljust(alignment_length, "-"),
        parsed["spades_orf_track"].ljust(alignment_length, "-"),
    )


def mapping_validation_selected_spades_orfs(path: Path) -> dict[str, str]:
    """Return selected SPAdes ORFs keyed by their original contig ID."""
    records: dict[str, str] = {}
    for header, sequence in parse_fasta(path):
        key = fasta_key(header)
        if "_ORF_" not in key:
            continue
        contig_id = key.split("_ORF_", 1)[1]
        records[contig_id] = normalize_dna(sequence)
    return records


def mapping_validation_best_spades_hits(path: Path) -> dict[str, dict[str, str]]:
    """Return the highest-scoring recorded SPAdes genotyping hit per contig."""
    best: dict[str, dict[str, str]] = {}
    for row in mapping_validation_load_csv(path):
        query = (row.get("qseqid", "") or "").strip()
        if not query:
            continue
        try:
            rank = (
                float(row.get("bitscore", "0") or 0),
                float(row.get("pident", "0") or 0),
            )
            current = best.get(query)
            current_rank = (
                (
                    float(current.get("bitscore", "0") or 0),
                    float(current.get("pident", "0") or 0),
                )
                if current is not None
                else None
            )
        except (TypeError, ValueError):
            continue
        if current_rank is None or rank > current_rank:
            best[query] = row
    return best


def mapping_validation_hit_value(
    hit: dict[str, object] | None,
    key: str,
) -> str:
    if not hit:
        return ""
    value = hit.get(key, "")
    return str(value) if value is not None else ""


def mapping_validation_blast_comparison_html(
    reported_genotype: str,
    spades_hit: dict[str, object] | None,
    mapping_hit: dict[str, object] | None,
) -> str:
    """Render side-by-side genotype-database hits and a concise agreement verdict."""
    spades_subject = mapping_validation_hit_value(spades_hit, "sseqid")
    mapping_subject = mapping_validation_hit_value(mapping_hit, "sseqid")
    spades_genotype = spades_subject.split("|", 1)[0] if spades_subject else ""
    mapping_genotype = mapping_subject.split("|", 1)[0] if mapping_subject else ""

    if not spades_hit or not mapping_hit:
        verdict_class = "unavailable"
        verdict = "Comparison unavailable: one or both BLAST hits are missing"
    elif spades_subject == mapping_subject:
        verdict_class = "same"
        verdict = "Same best genotyping-database sequence"
    elif spades_genotype == mapping_genotype:
        verdict_class = "same-genotype"
        verdict = "Same genotype, but a different best database sequence"
    else:
        verdict_class = "different"
        verdict = "Different genotype assignments"

    def row_html(label: str, hit: dict[str, object] | None) -> str:
        subject = mapping_validation_hit_value(hit, "sseqid") or "unavailable"
        genotype = subject.split("|", 1)[0] if subject != "unavailable" else "—"
        pident = mapping_validation_hit_value(hit, "pident") or "—"
        length = mapping_validation_hit_value(hit, "length") or "—"
        bitscore = mapping_validation_hit_value(hit, "bitscore") or "—"
        return (
            "<tr>"
            f"<th>{html.escape(label)}</th>"
            f"<td>{html.escape(subject)}</td>"
            f"<td>{html.escape(genotype)}</td>"
            f"<td>{html.escape(pident)}</td>"
            f"<td>{html.escape(length)}</td>"
            f"<td>{html.escape(bitscore)}</td>"
            "</tr>"
        )

    return (
        '<div class="genotype-origin">'
        "The reported RotaX genotype "
        f"<b>{html.escape(reported_genotype or 'unavailable')}</b> is based on the "
        "selected SPAdes sequence. The mapping CDS is BLASTed separately only for "
        "this validation comparison.</div>"
        '<table class="blast-comparison"><thead><tr>'
        "<th>Query</th><th>Best genotype-DB hit</th><th>Genotype</th>"
        "<th>Identity (%)</th><th>Aligned nt</th><th>Bitscore</th>"
        "</tr></thead><tbody>"
        + row_html("Selected SPAdes ORF", spades_hit)
        + row_html("Final mapping CDS", mapping_hit)
        + "</tbody></table>"
        f'<div class="blast-verdict {verdict_class}">{html.escape(verdict)}</div>'
    )


def mapping_validation_percentile(values: list[int], fraction: float) -> int:
    if not values:
        return 1
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, math.ceil(fraction * len(ordered)) - 1))
    return max(1, ordered[index])


def mapping_validation_logo_height(
    depth: int,
    normalizer: int,
    scale: str,
    maximum: float = 64.0,
) -> float:
    if depth <= 0:
        return 0.0
    if scale == "linear":
        fraction = min(1.0, depth / max(1, normalizer))
    elif scale == "sqrt":
        fraction = min(1.0, math.sqrt(depth / max(1, normalizer)))
    else:
        fraction = min(
            1.0, math.log1p(depth) / math.log1p(max(1, normalizer))
        )
    return max(5.0, maximum * fraction)


def mapping_validation_display_orf(
    orf: dict[str, object] | None,
    sequence_length: int,
    reversed_sequence: bool,
) -> dict[str, object] | None:
    if orf is None or not reversed_sequence:
        return orf
    start = int(orf.get("start", 0) or 0)
    end = int(orf.get("end", 0) or 0)
    transformed = dict(orf)
    transformed["start"] = sequence_length + 1 - end
    transformed["end"] = sequence_length + 1 - start
    transformed["strand"] = "-" if orf.get("strand") == "+" else "+"
    return transformed


def mapping_validation_overview_svg(item: dict[str, object]) -> str:
    """Return the former PDF-style alignment/CDS bars as an inline SVG."""
    consensus = str(item["aligned_consensus"])
    reference = str(item["aligned_reference"])
    alignment_length = max(len(consensus), len(reference), 1)
    consensus = consensus.ljust(alignment_length, "-")
    reference = reference.ljust(alignment_length, "-")
    consensus_coords = mapping_validation_aligned_coord_map(consensus)
    reference_coords = mapping_validation_aligned_coord_map(reference)
    fields = item.get("fields")
    summary = item.get("summary")
    fields = fields if isinstance(fields, dict) else {}
    summary = summary if isinstance(summary, dict) else {}

    width = 1120
    height = 270
    plot_left = 190.0
    plot_right = 1090.0
    plot_width = plot_right - plot_left
    thin_height = 7.0
    thick_height = 27.0
    lanes = (
        (
            "consensus",
            consensus,
            consensus_coords,
            item.get("consensus_orf"),
            88.0,
            "#1f6fb4",
            "My consensus",
            " ".join(
                value for value in (
                    str(fields.get("genotype", "")),
                    str(summary.get("ConsensusStatus") or fields.get("status", "")),
                ) if value
            ),
        ),
        (
            "reference",
            reference,
            reference_coords,
            item.get("reference_orf"),
            172.0,
            "#c1440e",
            "GenBank/reference",
            str(item.get("accession", "")),
        ),
    )
    chunks = [
        f'<svg class="overview-svg" viewBox="0 0 {width} {height}" '
        'role="img" aria-label="Consensus and reference alignment overview with CDS bars">',
        '<rect width="1120" height="270" fill="#fff"/>',
    ]

    def x_for_column(column: float) -> float:
        return plot_left + (column / alignment_length) * plot_width

    for (
        lane_name,
        sequence,
        coordinate_map,
        orf_value,
        lane_y,
        color,
        label,
        sublabel,
    ) in lanes:
        non_gap_columns = [
            column for column, base in enumerate(sequence) if base != "-"
        ]
        if not non_gap_columns:
            continue
        first = non_gap_columns[0]
        last = non_gap_columns[-1]
        first_x = x_for_column(first)
        last_x = x_for_column(last + 1)
        chunks.append(
            f'<rect x="{first_x:.2f}" y="{lane_y - thin_height / 2:.2f}" '
            f'width="{max(0.75, last_x - first_x):.2f}" height="{thin_height}" '
            f'fill="{color}" opacity=".86"/>'
        )

        orf = orf_value if isinstance(orf_value, dict) else None
        orf_columns: tuple[int, int] | None = None
        if orf is not None:
            try:
                start = min(int(orf.get("start", 0)), int(orf.get("end", 0)))
                end = max(int(orf.get("start", 0)), int(orf.get("end", 0)))
            except (TypeError, ValueError):
                start = end = 0
            matching = [
                column for column, position in enumerate(coordinate_map)
                if start <= position <= end and position > 0
            ]
            if matching:
                orf_columns = (matching[0], matching[-1])
                orf_x = x_for_column(matching[0])
                orf_end_x = x_for_column(matching[-1] + 1)
                chunks.append(
                    f'<rect x="{orf_x:.2f}" y="{lane_y - thick_height / 2:.2f}" '
                    f'width="{max(0.75, orf_end_x - orf_x):.2f}" height="{thick_height}" '
                    f'fill="{color}" stroke="#1e2933" stroke-width=".8" opacity=".96"/>'
                )

        # Preserve the black N/internal-gap markings used by the PDF figure.
        black_run_start: int | None = None
        for column in range(first, last + 2):
            black = (
                column <= last
                and (sequence[column] == "N" or sequence[column] == "-")
            )
            if black and black_run_start is None:
                black_run_start = column
            elif not black and black_run_start is not None:
                black_run_end = column - 1
                in_orf = (
                    orf_columns is not None
                    and black_run_start >= orf_columns[0]
                    and black_run_end <= orf_columns[1]
                )
                mark_height = thick_height if in_orf else thin_height
                mark_x = x_for_column(black_run_start)
                mark_end_x = x_for_column(column)
                chunks.append(
                    f'<rect x="{mark_x:.2f}" y="{lane_y - mark_height / 2:.2f}" '
                    f'width="{max(0.9, mark_end_x - mark_x):.2f}" '
                    f'height="{mark_height}" fill="#111"/>'
                )
                black_run_start = None

        chunks.extend([
            f'<text class="overview-lane-label" x="{plot_left - 18:.0f}" '
            f'y="{lane_y - 4:.0f}" fill="{color}">{html.escape(label)}</text>',
            f'<text class="overview-lane-detail" x="{plot_left - 18:.0f}" '
            f'y="{lane_y + 12:.0f}" fill="{color}">{html.escape(sublabel)}</text>',
        ])
        if orf is not None:
            cds_text = (
                f"CDS {start}-{end} "
                f"({orf.get('strand', '+')})"
            )
        else:
            cds_text = "CDS: unavailable"
        chunks.append(
            f'<text class="overview-cds-label" x="{plot_left - 18:.0f}" '
            f'y="{lane_y + 28:.0f}" fill="{color}">{html.escape(cds_text)}</text>'
        )

        label_columns = {first, last}
        if orf_columns is not None:
            label_columns.update(orf_columns)
        for column in sorted(label_columns):
            if not (0 <= column < len(coordinate_map)):
                continue
            position = coordinate_map[column]
            if position <= 0:
                continue
            label_x = x_for_column(column)
            text_y = lane_y - 21 if lane_name == "consensus" else lane_y + 28
            chunks.append(
                f'<text class="overview-position" x="{label_x:.2f}" '
                f'y="{text_y:.2f}">{position}</text>'
            )

    # Shared ruler: labels use reference coordinates where present.
    ruler_y = 230.0
    chunks.append(
        f'<line x1="{plot_left}" y1="{ruler_y}" x2="{plot_right}" '
        f'y2="{ruler_y}" stroke="#73808b" stroke-width="1"/>'
    )
    for tick_index in range(11):
        column = round(tick_index * (alignment_length - 1) / 10)
        position = (
            reference_coords[column]
            if reference_coords[column] > 0 else consensus_coords[column]
        )
        tick_x = x_for_column(column)
        chunks.extend([
            f'<line x1="{tick_x:.2f}" y1="{ruler_y}" x2="{tick_x:.2f}" '
            f'y2="{ruler_y + 6}" stroke="#73808b"/>',
            f'<text class="overview-ruler-label" x="{tick_x:.2f}" '
            f'y="{ruler_y + 19}">{position or ""}</text>',
        ])
    chunks.extend([
        '<text class="overview-axis-title" x="640" y="266">'
        'Position (ungapped; reference coordinates where available)</text>',
        '<text class="overview-legend" x="1090" y="24">black = N or internal gap</text>',
        "</svg>",
    ])
    return "".join(chunks)


def mapping_validation_svg(
    item: dict[str, object],
    depth_scale: str,
) -> tuple[str, dict[str, int]]:
    consensus = str(item["aligned_consensus"])
    reference = str(item["aligned_reference"])
    depth_track = str(item["aligned_depth"])
    entries = item["depth_entries"]
    allele_depths = item["allele_depths"]
    if not isinstance(entries, dict) or not isinstance(allele_depths, dict):
        raise TypeError("Invalid mapping validation depth tracks")
    alignment_length = len(consensus)
    consensus_coords = mapping_validation_aligned_coord_map(consensus)
    reference_coords = mapping_validation_aligned_coord_map(reference)
    depth_coords = mapping_validation_aligned_coord_map(depth_track)
    depth_length = max(entries) if entries else 0
    normalizer = mapping_validation_percentile(
        [
            int(entry.get("depth", 0))
            for entry in entries.values()
            if isinstance(entry, dict) and int(entry.get("depth", 0)) > 0
        ],
        0.99,
    )
    consensus_orf = item.get("consensus_orf")
    reference_orf = item.get("reference_orf")
    reference_reversed = bool(item.get("reference_reversed"))
    depth_reversed = bool(item.get("depth_reversed"))
    reference_length = sum(base != "-" for base in reference)
    column_width = 18
    logo_bottom = 126.0
    chunks = [
        f'<svg class="sequence-svg" viewBox="0 0 {max(1, alignment_length * column_width)} 190" '
        f'width="{max(1, alignment_length * column_width)}" height="190" '
        'role="img" aria-label="Aligned consensus, read-depth logo, and mapping reference">'
    ]
    stats = {
        "alignment": alignment_length, "mismatch": 0, "ambiguity": 0,
        "n": 0, "gap": 0, "low": 0, "zero": 0, "normalizer": normalizer,
    }
    for column in range(alignment_length):
        x = column * column_width
        consensus_base = consensus[column]
        reference_base = reference[column]
        consensus_position = consensus_coords[column]
        reference_display_position = reference_coords[column]
        reference_position = (
            reference_length + 1 - reference_display_position
            if reference_reversed and reference_display_position else reference_display_position
        )
        depth_position = depth_coords[column]
        if depth_position and depth_reversed:
            depth_position = depth_length + 1 - depth_position
        entry = entries.get(depth_position) if depth_position else None
        depth = int(entry.get("depth", 0)) if isinstance(entry, dict) else 0
        meets = bool(entry.get("meets")) if isinstance(entry, dict) else False
        ambiguous = (
            consensus_base in MAPPING_VALIDATION_IUPAC
            and len(MAPPING_VALIDATION_IUPAC[consensus_base]) > 1
        )
        n_call = consensus_base == "N"
        mismatch = (
            consensus_base in "ACGT" and reference_base in "ACGT"
            and consensus_base != reference_base
        )
        gap = (consensus_base == "-") != (reference_base == "-")
        low = entry is not None and not meets
        zero = entry is not None and depth == 0
        stats["ambiguity"] += int(ambiguous)
        stats["n"] += int(n_call)
        stats["mismatch"] += int(mismatch)
        stats["gap"] += int(gap)
        stats["low"] += int(low)
        stats["zero"] += int(zero)

        classes = ["base-column"]
        if ambiguous or n_call or mismatch or gap or low:
            classes.append("issue")
        if ambiguous:
            classes.append("ambiguous")
        if n_call:
            classes.append("n-call")

        in_consensus_orf = (
            isinstance(consensus_orf, dict) and consensus_position
            and int(consensus_orf.get("start", 0)) <= consensus_position
            <= int(consensus_orf.get("end", 0))
        )
        in_reference_orf = (
            isinstance(reference_orf, dict) and reference_display_position
            and int(reference_orf.get("start", 0)) <= reference_display_position
            <= int(reference_orf.get("end", 0))
        )
        backgrounds: list[str] = []
        if in_consensus_orf:
            backgrounds.append('<rect y="19" width="18" height="111" fill="#edf6fc"/>')
        if in_reference_orf:
            backgrounds.append('<rect y="130" width="18" height="38" fill="#fff4e8"/>')
        if mismatch:
            backgrounds.append('<rect width="18" height="190" fill="#ffdede" fill-opacity=".72"/>')
        elif gap:
            backgrounds.append('<rect width="18" height="190" fill="#eee5ff" fill-opacity=".72"/>')
        elif zero or n_call:
            backgrounds.append('<rect width="18" height="190" fill="#e9ecef" fill-opacity=".72"/>')
        elif ambiguous:
            backgrounds.append('<rect width="18" height="190" fill="#fff0b8" fill-opacity=".72"/>')
        elif low:
            backgrounds.append('<rect width="18" height="190" fill="#ffe8cf" fill-opacity=".72"/>')

        components = MAPPING_VALIDATION_IUPAC.get(consensus_base, ())
        counts = allele_depths.get(depth_position, {}) if depth_position else {}
        if counts and depth_reversed:
            counts = {
                MAPPING_VALIDATION_COMPLEMENT.get(base, base): count
                for base, count in counts.items()
            }
        selected_counts = {base: int(counts.get(base, 0)) for base in components}
        selected_total = sum(selected_counts.values())
        proportion_note = ""
        if len(components) > 1 and selected_total:
            proportions = {
                base: selected_counts[base] / selected_total for base in components
            }
            proportion_note = "VCF AD " + ", ".join(
                f"{base}:{selected_counts[base]}" for base in components
            )
        elif len(components) > 1:
            proportions = {base: 1 / len(components) for base in components}
            proportion_note = "equal IUPAC split; VCF AD unavailable"
        else:
            proportions = {components[0]: 1.0} if components else {}

        title = [
            f"alignment column {column + 1}",
            f"consensus {consensus_position or '-'}: {consensus_base}",
            f"GenBank/reference {reference_position or '-'}: {reference_base}",
            f"mapping-depth position {depth_position or '-'}: "
            f"{depth if entry is not None else 'n/a'}",
        ]
        if proportion_note:
            title.append(proportion_note)
        if n_call:
            title.append("N is unknown/masked, not an ambiguity call")
        chunks.append(
            f'<g class="{" ".join(classes)}" transform="translate({x},0)" '
            f'data-x="{x}" data-cons-pos="{consensus_position or ""}">'
            f"<title>{html.escape(chr(10).join(title))}</title>{''.join(backgrounds)}"
        )
        if consensus_position and consensus_position % 50 == 0:
            chunks.append(
                f'<line x1="9" y1="0" x2="9" y2="8" stroke="#6c757d"/>'
                f'<text class="coord" x="9" y="17">{consensus_position}</text>'
            )
        elif consensus_position and consensus_position % 10 == 0:
            chunks.append('<line x1="9" y1="2" x2="9" y2="7" stroke="#adb5bd"/>')
        shown_consensus = "·" if consensus_base == "-" else consensus_base
        chunks.append(
            f'<text class="sequence-letter{" ambiguous-call" if ambiguous else ""}" '
            f'x="9" y="39">{shown_consensus}</text>'
        )
        height = mapping_validation_logo_height(depth, normalizer, depth_scale)
        y_bottom = logo_bottom
        for base in sorted(components, key=lambda value: proportions.get(value, 0)):
            base_height = height * proportions.get(base, 0)
            if base_height < 2:
                continue
            chunks.append(
                f'<text class="logo-letter" x="9" y="{y_bottom:.2f}" '
                f'font-size="{base_height:.2f}" textLength="14" '
                f'lengthAdjust="spacingAndGlyphs" '
                f'fill="{MAPPING_VALIDATION_COLORS[base]}">{base}</text>'
            )
            y_bottom -= base_height
        if zero and consensus_base != "-":
            chunks.append(
                '<line x1="3" y1="125" x2="15" y2="125" '
                'stroke="#9aa0a6" stroke-width="2"/>'
            )
        shown_reference = "·" if reference_base == "-" else reference_base
        chunks.append(
            f'<text class="sequence-letter reference{" mismatch-letter" if mismatch else ""}" '
            f'x="9" y="154">{shown_reference}</text>'
        )
        if reference_position and reference_position % 50 == 0:
            chunks.append(
                f'<line x1="9" y1="160" x2="9" y2="168" stroke="#6c757d"/>'
                f'<text class="coord" x="9" y="182">{reference_position}</text>'
            )
        elif reference_position and reference_position % 10 == 0:
            chunks.append('<line x1="9" y1="161" x2="9" y2="166" stroke="#adb5bd"/>')
        chunks.append("</g>")
    chunks.append(
        f'<line x1="0" y1="126" x2="{max(1, alignment_length * column_width)}" '
        'y2="126" stroke="#c7cdd2"/></svg>'
    )
    return "".join(chunks), stats


def mapping_validation_orf_pair_svg(
    aligned_mapping: str,
    aligned_spades: str,
    mapping_cds_start: int,
) -> tuple[str, dict[str, int]]:
    """Render an interactive mapping-CDS versus SPAdes-ORF nucleotide frame."""
    alignment_length = max(len(aligned_mapping), len(aligned_spades), 1)
    mapping = aligned_mapping.ljust(alignment_length, "-")
    spades = aligned_spades.ljust(alignment_length, "-")
    mapping_coords = mapping_validation_aligned_coord_map(mapping)
    spades_coords = mapping_validation_aligned_coord_map(spades)
    column_width = 18
    svg_width = max(1, alignment_length * column_width)
    chunks = [
        f'<svg class="pair-sequence-svg" viewBox="0 0 {svg_width} 116" '
        f'width="{svg_width}" height="116" role="img" '
        'aria-label="Aligned final mapping CDS and selected SPAdes ORF">'
    ]
    stats = {
        "alignment": alignment_length,
        "mismatch": 0,
        "ambiguity": 0,
        "n": 0,
        "gap": 0,
        "issue": 0,
    }

    for column in range(alignment_length):
        x = column * column_width
        mapping_base = mapping[column]
        spades_base = spades[column]
        mapping_relative = mapping_coords[column]
        consensus_position = (
            mapping_cds_start + mapping_relative - 1
            if mapping_relative and mapping_cds_start > 0
            else 0
        )
        spades_position = spades_coords[column]
        ambiguity = (
            mapping_base in MAPPING_VALIDATION_IUPAC
            and len(MAPPING_VALIDATION_IUPAC[mapping_base]) > 1
        )
        n_call = mapping_base == "N"
        mismatch = (
            mapping_base in "ACGT"
            and spades_base in "ACGT"
            and mapping_base != spades_base
        )
        gap = (mapping_base == "-") != (spades_base == "-")
        issue = ambiguity or n_call or mismatch or gap
        stats["ambiguity"] += int(ambiguity)
        stats["n"] += int(n_call)
        stats["mismatch"] += int(mismatch)
        stats["gap"] += int(gap)
        stats["issue"] += int(issue)

        classes = ["base-column", "pair-base-column"]
        if issue:
            classes.append("issue")
        if ambiguity:
            classes.append("ambiguous")
        if n_call:
            classes.append("n-call")

        backgrounds: list[str] = []
        if mismatch:
            backgrounds.append(
                '<rect width="18" height="116" fill="#ffdede" fill-opacity=".78"/>'
            )
        elif gap:
            backgrounds.append(
                '<rect width="18" height="116" fill="#eee5ff" fill-opacity=".78"/>'
            )
        elif n_call:
            backgrounds.append(
                '<rect width="18" height="116" fill="#e9ecef" fill-opacity=".78"/>'
            )
        elif ambiguity:
            backgrounds.append(
                '<rect width="18" height="116" fill="#fff0b8" fill-opacity=".78"/>'
            )

        title = [
            f"alignment column {column + 1}",
            f"mapping CDS {mapping_relative or '-'}"
            + (
                f" (consensus {consensus_position})"
                if consensus_position
                else ""
            )
            + f": {mapping_base}",
            f"SPAdes ORF {spades_position or '-'}: {spades_base}",
        ]
        if ambiguity:
            title.append("mapping CDS contains an IUPAC ambiguity call")
        if n_call:
            title.append("mapping CDS contains N (unknown/masked)")
        if mismatch:
            title.append("nucleotide mismatch")
        if gap:
            title.append("alignment gap")

        chunks.append(
            f'<g class="{" ".join(classes)}" transform="translate({x},0)" '
            f'data-x="{x}" data-cons-pos="{consensus_position or ""}">'
            f"<title>{html.escape(chr(10).join(title))}</title>"
            + "".join(backgrounds)
        )
        if consensus_position and consensus_position % 50 == 0:
            chunks.append(
                '<line x1="9" y1="0" x2="9" y2="7" stroke="#6c757d"/>'
                f'<text class="coord" x="9" y="16">{consensus_position}</text>'
            )
        elif consensus_position and consensus_position % 10 == 0:
            chunks.append(
                '<line x1="9" y1="2" x2="9" y2="7" stroke="#adb5bd"/>'
            )

        shown_mapping = "·" if mapping_base == "-" else mapping_base
        shown_spades = "·" if spades_base == "-" else spades_base
        chunks.append(
            f'<text class="sequence-letter{" ambiguous-call" if ambiguity else ""}" '
            f'x="9" y="39">{shown_mapping}</text>'
            f'<text class="sequence-letter reference{" mismatch-letter" if mismatch else ""}" '
            f'x="9" y="82">{shown_spades}</text>'
        )
        if spades_position and spades_position % 50 == 0:
            chunks.append(
                '<line x1="9" y1="89" x2="9" y2="96" stroke="#6c757d"/>'
                f'<text class="coord" x="9" y="110">{spades_position}</text>'
            )
        elif spades_position and spades_position % 10 == 0:
            chunks.append(
                '<line x1="9" y1="90" x2="9" y2="95" stroke="#adb5bd"/>'
            )
        chunks.append("</g>")

    chunks.append("</svg>")
    return "".join(chunks), stats


def mapping_validation_orf_text(
    label: str,
    orf: dict[str, object] | None,
) -> str:
    if not orf:
        return f"{label}: unavailable"
    details = (
        f"{label}: {orf.get('start')}-{orf.get('end')} "
        f"({orf.get('strand', '+')})"
    )
    annotation = str(orf.get("product") or orf.get("gene") or "")
    return details + (f", {annotation}" if annotation else "")


def mapping_validation_difference_text(row: dict[str, str]) -> str:
    status = (row.get("ORFComparisonStatus", "") or "").strip()
    if status and status != "OK":
        return (
            f"SPAdes versus mapping ORF: not compared ({status}; "
            f"{row.get('ORFComparisonReason', '')})"
        )
    values = []
    for column, label in [
        ("SpadesVsMappingDifferenceColumns", "difference columns"),
        ("SpadesVsMappingMismatches", "mismatches"),
        ("SpadesVsMappingIndelColumns", "indel columns"),
        ("SpadesVsMappingIndelEvents", "indel events"),
    ]:
        value = row.get(column, "")
        if value != "":
            values.append(f"{value} {label}")
    return "SPAdes versus mapping ORF: " + (
        ", ".join(values) if values else "no comparison data"
    )


def build_mapping_validation_html(
    sample_dir: Path,
    args: argparse.Namespace,
    log_path: Path,
) -> Path | None:
    """Build one combined mapped-reference validation HTML for a sample."""
    mapping_dir = sample_dir / "mapping"
    consensus_path = mapping_dir / "rotavirus_mapping_consensus.fasta"
    reference_path = mapping_dir / "rotavirus_mapping_references.fasta"
    summary_path = mapping_dir / "rotavirus_mapping_summary.csv"
    features_path = mapping_dir / "rotavirus_mapping_consensus_features.csv"
    mapping_orf_path = mapping_dir / "rotavirus_mapping_orfs.fasta"
    spades_orf_path = sample_dir / "rotavirus_selected_orfs.fasta"
    spades_blast_path = sample_dir / "rotavirus_blast_evaluated.csv"
    depth_path = mapping_dir / "rotavirus_mapping_depth.tsv"
    vcf_path = mapping_dir / "rotavirus_mapping.variants.vcf.gz"
    output_path = mapping_dir / "rotavirus_mapping_validation.html"
    if not consensus_path.exists() or not reference_path.exists():
        log(
            f"{sample_dir.name}: mapping validation HTML not generated because "
            "consensus/reference FASTA is missing",
            log_path,
        )
        return None

    summary_rows = mapping_validation_load_csv(summary_path)
    summary_by_id = {
        row.get("ConsensusID", ""): row for row in summary_rows
        if row.get("ConsensusID", "")
    }
    feature_rows = mapping_validation_load_csv(features_path)
    consensus_orfs: dict[str, dict[str, object]] = {}
    for row in feature_rows:
        if (row.get("Feature", "") or "").upper() != "CDS":
            continue
        try:
            consensus_orfs[row.get("ConsensusID", "")] = {
                "start": int(row.get("Start", "")),
                "end": int(row.get("End", "")),
                "strand": row.get("Strand", "+") or "+",
                "gene": row.get("Gene", ""),
                "product": row.get("Product", ""),
            }
        except (TypeError, ValueError):
            continue
    depth_by_id = mapping_validation_load_depth(depth_path)
    references = mapping_validation_reference_records(reference_path)
    vcf_depths = mapping_validation_load_allele_depths(vcf_path)
    selected_spades_orfs = mapping_validation_selected_spades_orfs(spades_orf_path)
    spades_hits_by_contig = mapping_validation_best_spades_hits(spades_blast_path)

    mapping_orfs_by_segment_genotype: dict[
        tuple[str, str], list[dict[str, str]]
    ] = defaultdict(list)
    mapping_segment_by_id: dict[str, str] = {}
    for mapping_header, mapping_sequence in parse_fasta(mapping_orf_path):
        mapping_fields = mapping_validation_header_fields(mapping_header)
        mapping_id = mapping_fields["id"]
        mapping_segment = mapping_fields["segment"]
        mapping_genotype = mapping_fields["genotype"]
        mapping_orfs_by_segment_genotype[
            (mapping_segment, mapping_genotype)
        ].append({
            "id": mapping_id,
            "sequence": normalize_dna(mapping_sequence),
        })
        mapping_segment_by_id[mapping_id] = mapping_segment

    global_notes: list[str] = []
    if not vcf_path.exists():
        global_notes.append(
            "Final variant VCF is unavailable; ambiguous stacks use equal component heights."
        )
    email = str(getattr(args, "mapping_validation_email", "") or "")
    api_key = str(getattr(args, "mapping_validation_api_key", "") or "")
    if not email:
        global_notes.append(
            "No mapping-validation email was configured. The exact local mapping "
            "reference is shown, but reference CDS annotations were not downloaded "
            "from GenBank. Set --mapping-validation-email or NCBI_EMAIL."
        )
    genbank_cache: dict[str, tuple[str, dict[str, object] | None] | None] = {}
    blast_db = str(getattr(args, "blast_db", "") or "")
    mapping_cds_hits = run_genotype_orf_blast(
        mapping_dir,
        mapping_orf_path,
        blast_db,
        mapping_segment_by_id,
        log_path,
    )
    if mapping_segment_by_id and not mapping_cds_hits:
        global_notes.append(
            "The separate mapping-CDS genotyping BLAST returned no usable "
            "same-segment hits."
        )
    items: list[dict[str, object]] = []

    for header, consensus_sequence_raw in parse_fasta(consensus_path):
        fields = mapping_validation_header_fields(header)
        consensus_id = fields["id"]
        segment = fields["segment"]
        row = summary_by_id.get(consensus_id, {})
        accession = (
            row.get("MappingReferenceAccession", "") or fields["reference"]
        )
        reference_id = row.get("MappingReferenceID", "")
        local_reference = mapping_validation_choose_reference(
            references, segment, accession, reference_id
        )
        notes: list[str] = []
        if local_reference is None:
            global_notes.append(
                f"{segment}: could not match the consensus to a local mapping reference."
            )
            continue
        accession = accession or local_reference["accession"]
        reference_sequence = local_reference["sequence"]
        reference_source = "local mapping reference"
        reference_orf: dict[str, object] | None = None
        fetched = mapping_validation_fetch_genbank(
            accession, email, api_key, genbank_cache, notes
        )
        if fetched is not None:
            reference_sequence, reference_orf = fetched
            reference_source = "NCBI GenBank record"
        entries = depth_by_id.get(consensus_id, {})
        aligned = mapping_validation_run_alignment(
            normalize_dna(consensus_sequence_raw),
            reference_sequence,
            mapping_validation_depth_sequence(entries),
            log_path,
        )
        if aligned is None:
            global_notes.append(f"{segment}: MAFFT alignment failed.")
            continue
        (
            aligned_consensus,
            aligned_reference,
            aligned_depth,
            reference_reversed,
            depth_reversed,
        ) = aligned
        display_reference_orf = mapping_validation_display_orf(
            reference_orf,
            len(reference_sequence),
            reference_reversed,
        )
        reported_genotype = str(
            row.get("ORFCalledGenotype", "") or fields.get("genotype", "")
        )
        mapping_orf_candidates = mapping_orfs_by_segment_genotype.get(
            (segment, reported_genotype),
            [],
        )
        mapping_orf_record = (
            mapping_orf_candidates[0] if mapping_orf_candidates else None
        )
        mapping_orf_id = (
            mapping_orf_record.get("id", "")
            if isinstance(mapping_orf_record, dict)
            else ""
        )
        mapping_orf_sequence = (
            mapping_orf_record.get("sequence", "")
            if isinstance(mapping_orf_record, dict)
            else ""
        )
        spades_contig_id = str(row.get("OriginalSpadesContigID", "") or "")
        spades_orf_sequence = selected_spades_orfs.get(spades_contig_id, "")
        orf_pair_alignment = mapping_validation_run_orf_pair_alignment(
            mapping_orf_sequence,
            spades_orf_sequence,
            log_path,
        )
        if mapping_orf_sequence and spades_orf_sequence and orf_pair_alignment is None:
            notes.append("SPAdes-ORF versus mapping-CDS nucleotide alignment failed.")
        allele_depths = vcf_depths.get(local_reference["name"], {})
        items.append({
            "fields": fields,
            "summary": row,
            "notes": notes,
            "accession": accession,
            "reference_id": reference_id,
            "reference_source": reference_source,
            "reference_length": len(reference_sequence),
            "consensus_length": len(normalize_dna(consensus_sequence_raw)),
            "aligned_consensus": aligned_consensus,
            "aligned_reference": aligned_reference,
            "aligned_depth": aligned_depth,
            "reference_reversed": reference_reversed,
            "depth_reversed": depth_reversed,
            "depth_entries": entries,
            "allele_depths": allele_depths,
            "consensus_orf": consensus_orfs.get(consensus_id),
            "reference_orf": display_reference_orf,
            "reference_orf_original": reference_orf,
            "spades_contig_id": spades_contig_id,
            "spades_orf_sequence": spades_orf_sequence,
            "mapping_orf_id": mapping_orf_id,
            "mapping_orf_sequence": mapping_orf_sequence,
            "orf_pair_alignment": orf_pair_alignment,
            "spades_genotype_hit": spades_hits_by_contig.get(spades_contig_id),
            "mapping_genotype_hit": mapping_cds_hits.get(mapping_orf_id),
        })

    scale = str(getattr(args, "mapping_validation_depth_scale", "log") or "log")
    cards: list[str] = []
    navigation: list[str] = []
    for item_index, item in enumerate(items):
        fields = item["fields"]
        row = item["summary"]
        if not isinstance(fields, dict) or not isinstance(row, dict):
            continue
        segment = str(fields.get("segment", ""))
        card_id = f"segment-{re.sub(r'[^A-Za-z0-9_-]+', '_', segment)}-{item_index + 1}"
        navigation.append(f'<a href="#{card_id}">{html.escape(segment)}</a>')
        overview_svg = mapping_validation_overview_svg(item)
        svg, stats = mapping_validation_svg(item, scale)
        consensus_orf = item.get("consensus_orf")
        reference_orf = item.get("reference_orf_original")
        pair_alignment = item.get("orf_pair_alignment")
        pair_svg = ""
        pair_stats: dict[str, int] | None = None
        if (
            isinstance(pair_alignment, tuple)
            and len(pair_alignment) == 2
            and isinstance(consensus_orf, dict)
        ):
            try:
                mapping_cds_start = int(consensus_orf.get("start", 0) or 0)
            except (TypeError, ValueError):
                mapping_cds_start = 0
            pair_svg, pair_stats = mapping_validation_orf_pair_svg(
                str(pair_alignment[0]),
                str(pair_alignment[1]),
                mapping_cds_start,
            )
        blast_comparison = mapping_validation_blast_comparison_html(
            str(row.get("ORFCalledGenotype", "") or fields.get("genotype", "")),
            (
                item.get("spades_genotype_hit")
                if isinstance(item.get("spades_genotype_hit"), dict)
                else None
            ),
            (
                item.get("mapping_genotype_hit")
                if isinstance(item.get("mapping_genotype_hit"), dict)
                else None
            ),
        )
        metadata = [
            f"<b>Sample:</b> {html.escape(str(fields.get('sample', '')))}",
            f"<b>Consensus:</b> {int(item['consensus_length']):,} bp; "
            f"{html.escape(str(row.get('ConsensusStatus') or fields.get('status') or 'n/a'))}",
            f"<b>Mapping reference:</b> {html.escape(str(item['accession']))}; "
            f"{int(item['reference_length']):,} bp; "
            f"{html.escape(str(item['reference_source']))}",
            html.escape(mapping_validation_orf_text("Consensus CDS", consensus_orf if isinstance(consensus_orf, dict) else None)),
            html.escape(mapping_validation_orf_text("Reference CDS", reference_orf if isinstance(reference_orf, dict) else None)),
            html.escape(mapping_validation_difference_text(row)),
            f"<b>Reference coverage:</b> {html.escape(str(row.get('CoveragePct', 'n/a')))}%; "
            f"mean depth {html.escape(str(row.get('MeanDepth', 'n/a')))}, "
            f"median depth {html.escape(str(row.get('MedianDepth', 'n/a')))}",
            f"<b>Mapping-reference selection:</b> identity "
            f"{html.escape(str(row.get('MappingReferencePident', 'n/a')))}%; "
            f"bitscore {html.escape(str(row.get('MappingReferenceBitscore', 'n/a')))}",
        ]
        notes = item.get("notes", [])
        notes_html = "".join(
            f"<li>{html.escape(str(note))}</li>"
            for note in notes if str(note)
        ) if isinstance(notes, list) else ""
        cards.append(f"""
<section class="segment-card" id="{card_id}">
  <header>
    <h2>{html.escape(segment)} <span>{html.escape(str(fields.get('genotype', '')))}</span></h2>
    <div class="stats">
      {stats['alignment']:,} columns · {stats['mismatch']:,} mismatch(es) ·
      {stats['ambiguity']:,} ambiguity call(s) · {stats['n']:,} N/masked ·
      {stats['gap']:,} gap column(s) · {stats['low']:,} below depth threshold,
      including {stats['zero']:,} zero-depth · p99 depth {stats['normalizer']:,}
    </div>
  </header>
  <details open><summary>Validation and ORF information</summary>
    <ul>{''.join(f'<li>{value}</li>' for value in metadata)}</ul>
    {f'<ul class="warnings">{notes_html}</ul>' if notes_html else ''}
  </details>
  <details open><summary>Genotyping basis and separate mapping-CDS BLAST</summary>
    {blast_comparison}
  </details>
  <figure class="overview-figure">
    <figcaption>Alignment and CDS overview</figcaption>
    {overview_svg}
  </figure>
  <div class="comparison-frame">
    <h3>Mapping consensus versus mapping reference</h3>
    <div class="interactive-comparison">
      <div class="toolbar">
        <label>Consensus position <input class="jump-input" type="number" min="1"></label>
        <button type="button" class="jump-button">Jump</button>
        <button type="button" class="prev-issue">← Issue</button>
        <button type="button" class="next-issue">Issue →</button>
        <button type="button" class="prev-ambiguity">← Ambiguity</button>
        <button type="button" class="next-ambiguity">Ambiguity →</button>
        <span>Hover a base for coordinates, depth, and AD.</span>
      </div>
      <div class="track-shell">
        <div class="track-labels" aria-hidden="true">
          <div class="track-label consensus-label">Consensus</div>
          <div class="track-label logo-label">Depth logo</div>
          <div class="track-label reference-label">Reference</div>
        </div>
        <div class="sequence-scroll" tabindex="0">{svg}</div>
      </div>
    </div>
  </div>
  <div class="comparison-frame pair-frame">
    <h3>Selected SPAdes ORF versus final mapping CDS</h3>
    {(
      f'''<div class="stats pair-stats">
        {pair_stats['alignment']:,} columns · {pair_stats['mismatch']:,} mismatch(es) ·
        {pair_stats['ambiguity']:,} mapping ambiguity call(s) ·
        {pair_stats['n']:,} mapping N/masked · {pair_stats['gap']:,} gap column(s)
      </div>
      <div class="interactive-comparison">
        <div class="toolbar">
          <label>Consensus position <input class="jump-input" type="number" min="1"></label>
          <button type="button" class="jump-button">Jump</button>
          <button type="button" class="prev-issue">← Issue</button>
          <button type="button" class="next-issue">Issue →</button>
          <button type="button" class="prev-ambiguity">← Ambiguity</button>
          <button type="button" class="next-ambiguity">Ambiguity →</button>
          <span>Coordinates refer to the final mapping consensus.</span>
        </div>
        <div class="track-shell pair-track-shell">
          <div class="track-labels pair-track-labels" aria-hidden="true">
            <div class="track-label pair-mapping-label">Mapping CDS</div>
            <div class="track-label pair-spades-label">SPAdes ORF</div>
          </div>
          <div class="sequence-scroll pair-sequence-scroll" tabindex="0">{pair_svg}</div>
        </div>
      </div>'''
	      if pair_stats is not None
	      else '<div class="comparison-unavailable">SPAdes ORF or mapping CDS alignment is unavailable.</div>'
	    )}
	  </div>
	</section>""")

    notes_html = "".join(f"<li>{html.escape(note)}</li>" for note in global_notes)
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RotaX mapped-reference validation — {html.escape(sample_dir.name)}</title>
<style>
:root{{--ink:#17212b;--muted:#66727e;--line:#d8dee4;--accent:#176b87}}
*{{box-sizing:border-box}} body{{margin:0;background:#f4f7f9;color:var(--ink);
font:14px/1.45 system-ui,-apple-system,Segoe UI,sans-serif}}
main{{width:min(1500px,98vw);margin:auto;padding:22px 0 60px}}
h1{{margin:0}} .subtitle,.stats{{color:var(--muted)}} .intro,.global-notes,
.segment-card{{background:#fff;border:1px solid var(--line);border-radius:9px}}
.intro,.global-notes{{padding:12px 16px;margin:12px 0}} .global-notes{{border-left:5px solid #d99b20}}
nav{{display:flex;gap:7px;flex-wrap:wrap;position:sticky;top:0;z-index:5;
padding:9px;margin:10px 0;background:rgba(244,247,249,.94)}}
nav a{{padding:4px 10px;border:1px solid #bfcbd3;border-radius:16px;background:#fff;
color:var(--accent);text-decoration:none;font-weight:700}}
.segment-card{{margin:0 0 18px;overflow:hidden}} .segment-card header{{padding:12px 16px 8px}}
h2{{margin:0}} h2 span{{font-size:15px;color:var(--accent)}} details{{margin:0 16px 9px}}
.warnings{{color:#8a5400}} .toolbar{{display:flex;align-items:center;flex-wrap:wrap;gap:7px;
padding:8px 12px;border-top:1px solid var(--line);background:#f8fafb}}
.overview-figure{{margin:10px 16px 14px;padding:8px 10px 2px;border:1px solid var(--line);
border-radius:7px;background:#fff}} .overview-figure figcaption{{font-weight:750;margin:0 0 3px}}
.overview-svg{{display:block;width:100%;height:auto;min-height:180px}}
.overview-lane-label{{text-anchor:end;font-size:14px;font-weight:800}}
.overview-lane-detail{{text-anchor:end;font-size:11px;font-weight:650}}
.overview-cds-label{{text-anchor:end;font-size:10px}}
.overview-position{{text-anchor:middle;font-size:9px;fill:#33424f}}
.overview-ruler-label{{text-anchor:middle;font-size:9px;fill:#53616c}}
	.overview-axis-title{{text-anchor:middle;font-size:10px;fill:#53616c}}
	.overview-legend{{text-anchor:end;font-size:10px;fill:#333}}
	.genotype-origin{{margin:7px 0 10px;padding:8px 10px;border-left:4px solid var(--accent);
	background:#f0f8fb}} .blast-comparison{{width:100%;border-collapse:collapse;margin:7px 0}}
	.blast-comparison th,.blast-comparison td{{padding:6px 8px;border:1px solid var(--line);
	text-align:left;overflow-wrap:anywhere}} .blast-comparison thead th{{background:#f2f5f7}}
	.blast-verdict{{display:inline-block;margin:2px 0 7px;padding:5px 9px;border-radius:5px;font-weight:750}}
	.blast-verdict.same{{background:#d9f3e3;color:#165a31}}
	.blast-verdict.same-genotype{{background:#e7f0fb;color:#184f82}}
	.blast-verdict.different{{background:#fde0dc;color:#8b1d14}}
	.blast-verdict.unavailable{{background:#eceff1;color:#53616c}}
	.comparison-frame{{margin:13px 16px 16px;border:1px solid var(--line);border-radius:7px;overflow:hidden}}
	.comparison-frame h3{{margin:0;padding:9px 12px;background:#eef3f6;font-size:15px}}
	.pair-stats{{padding:6px 12px;background:#fbfcfd;border-top:1px solid var(--line)}}
	.comparison-unavailable{{padding:14px;color:var(--muted)}}
	.toolbar input{{width:90px;padding:4px}} button{{padding:5px 9px;border:1px solid #aebbc5;
	border-radius:5px;background:#fff;cursor:pointer}} .toolbar span{{margin-left:auto;color:var(--muted);font-size:12px}}
	.track-shell{{display:grid;grid-template-columns:94px minmax(0,1fr);border-top:1px solid var(--line)}}
	.track-labels{{height:190px;position:relative;background:#f8fafb;border-right:1px solid var(--line);
	box-shadow:3px 0 6px rgba(30,45,60,.09);z-index:2}} .track-label{{position:absolute;right:9px;font-weight:700;font-size:12px}}
	.consensus-label{{top:24px}} .logo-label{{top:82px}} .reference-label{{top:139px}}
	.sequence-scroll{{overflow-x:auto;overflow-y:hidden;height:208px;background:#fff}}
	.sequence-svg{{display:block;max-width:none;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}}
	.pair-track-labels{{height:116px}} .pair-mapping-label{{top:25px}}
	.pair-spades-label{{top:68px}} .pair-sequence-scroll{{height:134px}}
	.pair-sequence-svg{{display:block;max-width:none;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}}
	.sequence-letter{{text-anchor:middle;font-size:15px;font-weight:750;fill:#14212c}}
.sequence-letter.reference{{fill:#53616c}} .ambiguous-call{{fill:#a05a00}}
.mismatch-letter{{fill:#b42318;font-weight:900}} .logo-letter{{text-anchor:middle;font-weight:900}}
.coord{{text-anchor:middle;font-size:8px;fill:#596773}} .base-column.issue:hover rect{{stroke:#212b35;stroke-width:1.5px}}
.base-column.active rect{{stroke:#006d8f;stroke-width:2.5px}}
@media(max-width:700px){{main{{width:100%}}.segment-card{{border-radius:0}}.track-shell{{grid-template-columns:78px minmax(0,1fr)}}}}
</style></head><body><main>
	<h1>RotaX mapped-reference validation</h1>
	<div class="subtitle">Sample {html.escape(sample_dir.name)} · generated {time.strftime('%Y-%m-%d %H:%M:%S')}</div>
	<div class="intro">
	<b>The reported RotaX genotype is based on the selected SPAdes sequence.</b>
	For validation, the final mapping CDS is BLASTed separately against the same genotyping database,
	and both best hits are shown together for each segment.
	The first interactive frame compares the called mapping consensus with its mapping reference above a
	depth-scaled nucleotide logo. The second frame directly compares the selected SPAdes ORF with the final
	mapping CDS, using final-consensus coordinates for navigation.
	The report also includes mapping/reference metadata, CDS coordinates, coverage statistics,
	a compact alignment/CDS bar overview, and the complete per-base alignments.
	The overview reproduces the former PDF figure: thick bars mark CDS regions and black marks N or internal gaps.
Letter-stack height uses a
<b>{html.escape(scale)}</b> depth transform capped at segment p99. IUPAC components use VCF AD proportions.
N is unknown/masked and is not counted as ambiguity. Pale blue marks consensus CDS; pale orange marks reference CDS.
</div>
{f'<div class="global-notes"><b>Notes</b><ul>{notes_html}</ul></div>' if notes_html else ''}
<nav>{''.join(navigation)}</nav>
{''.join(cards) if cards else '<div class="global-notes">No segments could be aligned for validation.</div>'}
</main><script>
	document.querySelectorAll('.interactive-comparison').forEach(panel=>{{
	 const scroll=panel.querySelector('.sequence-scroll');
	 const groups=[...panel.querySelectorAll('.base-column')];
	 const issues=groups.filter(g=>g.classList.contains('issue'));
	 const ambiguities=groups.filter(g=>g.classList.contains('ambiguous'));
	 let issue=-1,ambiguity=-1;
	 const focus=g=>{{if(!g)return;groups.forEach(x=>x.classList.remove('active'));g.classList.add('active');
	  scroll.scrollTo({{left:Math.max(0,Number(g.dataset.x)-scroll.clientWidth/2),behavior:'smooth'}})}};
	 const jumpInput=panel.querySelector('.jump-input');
	 const jumpButton=panel.querySelector('.jump-button');
	 jumpButton.onclick=()=>focus(groups.find(g=>g.dataset.consPos===jumpInput.value));
	 jumpInput.onkeydown=e=>{{if(e.key==='Enter')jumpButton.click()}};
	 panel.querySelector('.next-issue').onclick=()=>{{if(issues.length)focus(issues[issue=(issue+1)%issues.length])}};
	 panel.querySelector('.prev-issue').onclick=()=>{{if(issues.length)focus(issues[issue=(issue-1+issues.length)%issues.length])}};
	 panel.querySelector('.next-ambiguity').onclick=()=>{{if(ambiguities.length)focus(ambiguities[ambiguity=(ambiguity+1)%ambiguities.length])}};
	 panel.querySelector('.prev-ambiguity').onclick=()=>{{if(ambiguities.length)focus(ambiguities[ambiguity=(ambiguity-1+ambiguities.length)%ambiguities.length])}};
	}});
</script></body></html>"""
    output_path.write_text(document, encoding="utf-8")
    log(
        f"{sample_dir.name}: generated combined mapping validation HTML "
        f"{output_path} with {len(items)} segment(s)",
        log_path,
    )
    return output_path


def prepare_fasta_input_sample(fasta_input: Path, parent_dir: Path, log_path: Path) -> Path | None:
    if not fasta_input.exists():
        log(f"FASTA input not found: {fasta_input}", log_path)
        return None

    records = parse_fasta(fasta_input)
    if not records:
        log(f"FASTA input contains no sequences: {fasta_input}", log_path)
        return None

    sample_name = fasta_input.stem
    sample_dir = parent_dir / sample_name
    sample_dir.mkdir(exist_ok=True)

    # For FASTA-only mode, do not apply SPAdes coverage filtering. The BLAST
    # length and percent-identity criteria are still evaluated downstream.
    write_fasta(records, sample_dir / "contigs.fasta")
    write_fasta(records, sample_dir / "contigs500.fasta")
    write_fasta(records, sample_dir / "contigs500_COV.fasta")

    log(f"Prepared FASTA-only sample folder {sample_dir.name} with {len(records)} sequence(s)", log_path)
    return sample_dir


def full_length_low_identity_x_report(sample_dir: Path, log_path: Path) -> None:
    blast_file = sample_dir / "rotavirus_blast_evaluated.csv"
    output_file = sample_dir / "rotavirus_full_length_low_identity_x_hits.csv"

    fields = [
        "Folder Name",
        "FullLengthLowIdentityHitCount",
        "BlastHeader",
        "BlastLines",
    ]

    blast_headers = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "sstrand", "slen",
    ]
    blast_header_line = "\t".join(blast_headers)

    if not blast_file.exists():
        with output_file.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=fields)
            writer.writeheader()
        return

    rows: list[dict[str, str]] = []
    with blast_file.open("r", newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(row)

    accepted_genes: set[str] = set()
    for row in rows:
        evaluation = row.get("Evaluation", "")
        if evaluation and evaluation != "not accepted":
            gene, _number = extract_gene_and_number(evaluation)
            if gene:
                accepted_genes.add(gene)

    low_identity_lines: list[str] = []
    for row in rows:
        if row.get("Evaluation", "") != "not accepted":
            continue

        sseqid = row.get("sseqid", "")
        genotype = sseqid.split("|", 1)[0] if sseqid else ""
        gene = genotype[0] if genotype else ""
        if gene not in BASE_THRESHOLDS or gene in accepted_genes:
            continue

        try:
            pident = float(row.get("pident", ""))
            length = int(row.get("length", ""))
            slen = int(row.get("slen", ""))
        except ValueError:
            continue

        is_full = length >= slen
        identity_failed = pident < BASE_THRESHOLDS[gene]
        if not (is_full and identity_failed):
            continue

        line = "\t".join(row.get(h, "") for h in blast_headers)
        low_identity_lines.append(line)

    with output_file.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        if low_identity_lines:
            writer.writerow({
                "Folder Name": sample_dir.name,
                "FullLengthLowIdentityHitCount": len(low_identity_lines),
                "BlastHeader": blast_header_line,
                "BlastLines": " || ".join(low_identity_lines),
            })

    log(f"Generated {output_file.name}", log_path)


def combine_low_identity_x_reports(parent_dir: Path, log_path: Path) -> Path | None:
    output_file = parent_dir / "rotax_full_length_low_identity_x_report.csv"
    report_files = sorted(parent_dir.glob("*/rotavirus_full_length_low_identity_x_hits.csv"))

    fields = [
        "Folder Name",
        "FullLengthLowIdentityHitCount",
        "BlastHeader",
        "BlastLines",
    ]

    wrote_any = False
    with output_file.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()

        for report in report_files:
            with report.open("r", newline="", encoding="utf-8", errors="replace") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    if row.get("Folder Name"):
                        writer.writerow({field: row.get(field, "") for field in fields})
                        wrote_any = True

    if wrote_any:
        log(f"Generated final full-length low-identity X report: {output_file.name}", log_path)
    else:
        log(f"Generated empty final full-length low-identity X report: {output_file.name}", log_path)

    return output_file


def safe_remove(path: Path, log_path: Path) -> None:
    try:
        if path.is_dir():
            shutil.rmtree(path)
        elif path.exists():
            path.unlink()
    except OSError as exc:
        log(f"Could not remove {path}: {exc}", log_path)


def cleanup_sample(sample_dir: Path, log_path: Path, skip_assembly: bool = False) -> None:
    for info in ASSEMBLY_OPTIONS.values():
        safe_remove(sample_dir / str(info["work_dir"]), log_path)

    names = [
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
    ]
    # In --skip-assembly mode the trimming/dedup stage is not re-run, so these
    # reads cannot be regenerated and mapping depends on them. Keep them so a
    # later --skip-assembly run does not lose its trimmed reads.
    if skip_assembly:
        protected = {
            "R1_paired.fastq.gz",
            "R2_paired.fastq.gz",
            "deduped_R1.fastq.gz",
            "deduped_R2.fastq.gz",
        }
        names = [n for n in names if n not in protected]
    for name in names:
        safe_remove(sample_dir / name, log_path)



def existing_rotax_result_paths(
    sample_dir: Path,
    ignore_names: set[str] | None = None,
) -> list[Path]:
    """Return RotaX-created files/folders already present in a sample folder.

    ``ignore_names`` is used by --skip-assembly so reusable input files such as
    contigs500.fasta are not mistaken for completed downstream results.
    """
    paths: list[Path] = []
    ignored = ignore_names or set()

    for name in ROTAX_SAMPLE_OUTPUT_FILES:
        if name in ignored:
            continue
        path = sample_dir / name
        if path.exists():
            paths.append(path)

    for info in ASSEMBLY_OPTIONS.values():
        path = sample_dir / str(info["work_dir"])
        if path.exists() and path.name not in ignored:
            paths.append(path)

    return paths


def remove_existing_rotax_results(
    sample_dir: Path,
    log_path: Path,
    preserve_names: set[str] | None = None,
) -> None:
    """Remove old RotaX outputs before rerunning a sample in overwrite mode.

    In --skip-assembly mode, contigs500.fasta is an input rather than an output
    and must survive an overwrite rerun.
    """
    preserved = preserve_names or set()
    paths = existing_rotax_result_paths(sample_dir, ignore_names=preserved)
    if not paths:
        return

    if preserved:
        kept = ", ".join(sorted(preserved))
        log(f"{sample_dir.name}: preserving reusable input file(s): {kept}", log_path)
    log(f"{sample_dir.name}: removing {len(paths)} existing downstream RotaX output file(s)/folder(s) before overwrite", log_path)
    for path in paths:
        safe_remove(path, log_path)


def ask_existing_results_action(
    folders_with_results: list[Path],
    args: argparse.Namespace,
    log_path: Path,
    ignore_names: set[str] | None = None,
) -> str:
    """Return 'overwrite' or 'skip' for all folders with existing downstream results."""
    if args.overwrite_existing:
        log("Existing RotaX results found; --overwrite-existing was used", log_path)
        return "overwrite"
    if args.skip_existing:
        log("Existing RotaX results found; --skip-existing was used", log_path)
        return "skip"

    print("\nExisting RotaX results were found in these sample folder(s):", flush=True)
    for folder in folders_with_results:
        found = existing_rotax_result_paths(folder, ignore_names=ignore_names)
        preview = ", ".join(path.name for path in found[:5])
        if len(found) > 5:
            preview += f", ... ({len(found)} total)"
        print(f"  - {folder.name}: {preview}", flush=True)

    print(
        "\nChoose once for all folders with existing results:\n"
        "  O = overwrite old RotaX files/folders and rerun these samples\n"
        "  S = skip folders with existing results\n",
        flush=True,
    )

    while True:
        answer = input("Overwrite existing RotaX results, or skip these folders? [O/S]: ").strip().lower()
        if answer in {"o", "overwrite"}:
            log("User selected: overwrite folders with existing RotaX results", log_path)
            return "overwrite"
        if answer in {"s", "skip"}:
            log("User selected: skip folders with existing RotaX results", log_path)
            return "skip"
        print("Please type O to overwrite, or S to skip.", flush=True)


def apply_existing_results_choice(folders: list[Path], args: argparse.Namespace, log_path: Path) -> list[Path]:
    """Ask once about existing downstream results and return folders to process.

    For --skip-assembly directory mode, contigs500.fasta is reusable input. It is
    excluded from the existing-result prompt and preserved if overwrite is chosen.
    """
    preserve_names: set[str] = set()
    if args.skip_assembly and not args.fasta_input:
        # In --skip-assembly mode the assembly/trimming stage is not re-run, so
        # the trimmed and deduped reads it would produce cannot be regenerated.
        # Preserve them (like contigs500.fasta) so an --overwrite-existing rerun
        # does not delete the very reads mapping depends on and force a raw-read
        # fallback.
        preserve_names.update({
            "contigs500.fasta",
            "R1_paired.fastq.gz",
            "R2_paired.fastq.gz",
            "R1_paired.repaired.fastq.gz",
            "R2_paired.repaired.fastq.gz",
            "deduped_R1.fastq.gz",
            "deduped_R2.fastq.gz",
        })

    folders_with_results = [
        folder
        for folder in folders
        if existing_rotax_result_paths(folder, ignore_names=preserve_names)
    ]
    if not folders_with_results:
        return folders

    action = ask_existing_results_action(
        folders_with_results,
        args,
        log_path,
        ignore_names=preserve_names,
    )

    if action == "skip":
        skip_names = {folder.name for folder in folders_with_results}
        log(f"Skipping {len(skip_names)} folder(s) with existing downstream RotaX results", log_path)
        return [folder for folder in folders if folder.name not in skip_names]

    for folder in folders_with_results:
        remove_existing_rotax_results(folder, log_path, preserve_names=preserve_names)
    return folders


def process_sample(
    sample_dir: Path,
    args: argparse.Namespace,
    log_path: Path,
) -> bool:
    log(f"Processing sample folder: {sample_dir.name}", log_path)

    try:
        if not args.skip_assembly:
            assembled = run_assembly(
                sample_dir,
                args.threads,
                log_path,
                args.assemblies,
                args.fastq_repair,
                args.bbnorm_threads,
            )
            if not assembled and not (sample_dir / "contigs500.fasta").exists():
                log(f"No assembly output for {sample_dir.name}; skipping genotyping", log_path)
                return False
        else:
            reusable_contigs = sample_dir / "contigs500.fasta"
            log(
                f"{sample_dir.name}: --skip-assembly active; reusing {reusable_contigs}",
                log_path,
            )

        if not (sample_dir / "contigs500.fasta").exists():
            log(f"{sample_dir.name}: contigs500.fasta not found; skipping genotyping", log_path)
            return False

        if args.fasta_input:
            log(f"{sample_dir.name}: FASTA-input mode; skipping SPAdes coverage filtering", log_path)
        else:
            filter_by_coverage(sample_dir, args.min_cov, log_path)

        if not blast_rotavirus(sample_dir, args.blast_db, log_path):
            log(f"{sample_dir.name}: rotavirus BLAST failed or produced no usable input", log_path)
            return False

        if not is_nonempty(sample_dir / "rotavirus_blast_raw.tsv"):
            log(f"{sample_dir.name}: rotavirus_blast_raw.tsv is empty; writing empty genotype report", log_path)
            write_empty_genotype_report(sample_dir)
            if bool(getattr(args, "spades_full_segment_qc", False)):
                write_empty_genbank_reports(sample_dir, args, log_path, "empty_rotavirus_blast_raw")
            write_empty_mapping_outputs(sample_dir, "empty_rotavirus_blast_raw", log_path)
            full_length_low_identity_x_report(sample_dir, log_path)
            return True

        evaluate_blast(sample_dir, log_path)
        summarize_blast(sample_dir, log_path)
        genotype_summary(sample_dir, log_path)
        add_vaccine_match_columns(sample_dir, log_path)
        selected_contigs(sample_dir, log_path)
        selected_orfs(sample_dir, log_path)
        if bool(getattr(args, "spades_full_segment_qc", False)):
            evaluate_genbank_full_segments(sample_dir, args, log_path)
        run_mapping_consensus(sample_dir, args, log_path)
        if (
            bool(getattr(args, "mapping_consensus", True))
            and bool(getattr(args, "mapping_validation", True))
        ):
            try:
                build_mapping_validation_html(sample_dir, args, log_path)
            except Exception as exc:  # noqa: BLE001
                # Validation is a review output and must not invalidate otherwise
                # successful genotyping/mapping results.
                log(
                    f"WARNING: {sample_dir.name} mapping validation HTML failed: {exc}",
                    log_path,
                )
        full_length_low_identity_x_report(sample_dir, log_path)

        return True

    except Exception as exc:
        log(f"ERROR in {sample_dir.name}: {exc}", log_path)
        return False

    finally:
        if args.cleanup:
            cleanup_sample(sample_dir, log_path, skip_assembly=bool(getattr(args, "skip_assembly", False)))


def write_empty_genotype_report(sample_dir: Path) -> None:
    output = sample_dir / "rotax_sample_genotyping_report.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_GENOTYPING_FIELDS)
        writer.writeheader()
        writer.writerow({"Folder Name": sample_dir.name, "Extra Information": "None"})


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
        return "UTILSTREKV"

    g_part, p_part = g_p_genotype.split("-", 1)
    if re.fullmatch(r"G\d+", g_part) and (
        re.fullmatch(r"P\d+", p_part) or re.fullmatch(r"P\[\d+\]", p_part)
    ):
        return "RVA"

    return "UTILSTREKV"


def equine_g3_flag(sample_dir: Path) -> str:
    """Flag the G3 hit selected by the existing full/coverage/length rules.

    BLAST sseqid ends at whitespace, so a reference header beginning
    'G3|Equine-like ' is normally represented as exactly 'G3|Equine-like'.
    Accept the full prefix too for evaluated reports retaining full headers.
    Missing/rejected G3 hits and non-equine selected hits remain blank.
    """
    candidate = read_selected_candidate_rows(sample_dir).get("G3", {})
    subject = str(candidate.get("sseqid", ""))
    return "y" if subject == "G3|Equine-like" or subject.startswith("G3|Equine-like ") else ""


def combine_sample_reports(parent_dir: Path, log_path: Path) -> Path | None:
    output_file = parent_dir / "rotax_all_samples_genotyping_report.csv"
    report_files = sorted(parent_dir.glob("*/rotax_sample_genotyping_report.csv"))
    if not report_files:
        log("No sample rotax_sample_genotyping_report.csv files found; no final report created", log_path)
        return None

    with output_file.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=ALL_GENOTYPING_FIELDS)
        writer.writeheader()
        for report in report_files:
            with report.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                for source in csv.DictReader(handle):
                    row = {field: source.get(field, "") for field in SAMPLE_GENOTYPING_FIELDS}
                    constellation = source.get("Genotype Constellation", source.get("Genotype", ""))
                    display, _ = coinfection_display(constellation, source.get("Extra Information", ""))
                    row["Genotype"] = constellation
                    row["LabwareGenotype"] = display
                    # Keep category based on the original resolved primary G/P call.
                    row["Art"] = infer_result_category(g_p_genotype_from_constellation(constellation))
                    row["G3 equine-like"] = equine_g3_flag(report.parent)
                    # Old VP6 values were averages without a closest-match check.
                    # Recompute all RotaTeq fields when combining pre-v32 reports.
                    if "Genotype Constellation" not in source and "LabwareGenotype" not in source:
                        blast_file = report.parent / "rotavirus_blast_evaluated.csv"
                        blast_rows = []
                        if blast_file.exists():
                            with blast_file.open("r", encoding="utf-8", errors="replace", newline="") as blast:
                                blast_rows = list(csv.DictReader(blast))
                        else:
                            log(f"WARNING: {report.parent.name}: no evaluated BLAST evidence; RotaTeq fields left blank", log_path)
                        row.update(rotateq_closest_identities(report.parent, blast_rows))
                    writer.writerow(row)
    log(f"Generated final all-samples genotyping report: {output_file.name}", log_path)
    return output_file



def combine_genbank_outputs(parent_dir: Path, log_path: Path) -> None:
    report_files = sorted(parent_dir.glob("*/rotavirus_genbank_full_segment_report.csv"))
    combined_report = parent_dir / "rotax_genbank_full_segment_report.csv"
    all_full_fasta = parent_dir / "rotax_genbank_clean_all_full_segments.fasta"
    complete_fasta = parent_dir / "rotax_genbank_clean_complete_11seg_samples.fasta"

    all_rows: list[dict[str, str]] = []
    if report_files:
        with combined_report.open("w", newline="", encoding="utf-8") as out:
            writer: csv.DictWriter | None = None
            for report in report_files:
                with report.open("r", newline="", encoding="utf-8", errors="replace") as handle:
                    reader = csv.DictReader(handle)
                    if reader.fieldnames is None:
                        continue
                    if writer is None:
                        writer = csv.DictWriter(out, fieldnames=reader.fieldnames)
                        writer.writeheader()
                    for row in reader:
                        writer.writerow(row)
                        all_rows.append(row)
        log(f"Generated combined GenBank full-segment report: {combined_report.name}", log_path)
    else:
        with combined_report.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=genbank_report_fields())
            writer.writeheader()
        log("No per-sample GenBank reports found; wrote empty combined report", log_path)

    complete_samples = {
        row.get("Sample", "") for row in all_rows
        if row.get("Complete11SegmentSample") == "YES" and row.get("Status") == "GENBANK_READY_FULL"
    }

    all_records: list[tuple[str, str]] = []
    complete_records: list[tuple[str, str]] = []
    for fasta in sorted(parent_dir.glob("*/rotavirus_genbank_clean_full_segments.fasta")):
        sample = fasta.parent.name
        records = parse_fasta(fasta)
        all_records.extend(records)
        if sample in complete_samples:
            complete_records.extend(records)

    write_fasta(all_records, all_full_fasta)
    write_fasta(complete_records, complete_fasta)
    log(f"Generated {all_full_fasta.name} with {len(all_records)} full sequences", log_path)
    log(f"Generated {complete_fasta.name} with {len(complete_records)} sequences from complete 11-segment samples", log_path)


def combine_mapping_outputs(parent_dir: Path, log_path: Path) -> None:
    """Combine only the final, human-facing mapping outputs."""
    consensus_all = parent_dir / "rotax_mapping_consensus_all.fasta"
    consensus_full = parent_dir / "rotax_mapping_consensus_full.fasta"
    combined_summary = parent_dir / "rotax_mapping_summary.csv"
    combined_features = parent_dir / "rotax_mapping_consensus_features.csv"

    # Remove obsolete v10 parent-level mapping reports so they cannot be mistaken
    # for current results after an overwrite run.
    for old_name in [
        "rotax_mapping_consensus_all_segments.fasta",
        "rotax_mapping_consensus_full_segments.fasta",
        "rotax_mapping_genbank_metadata.csv",
        "rotax_mapping_genbank_feature_coordinates.csv",
        "rotax_mapping_orf_comparison_report.csv",
        "rotax_mapping_reference_selection.csv",
        "rotax_mapping_consensus_genbank_ready.fasta",
        "rotax_mapping_consensus_genbank.fasta",
    ]:
        safe_remove(parent_dir / old_name, log_path)

    summary_rows: list[dict[str, str]] = []
    for summary_path in sorted(parent_dir.glob("*/mapping/rotavirus_mapping_summary.csv")):
        with summary_path.open("r", newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                summary_rows.append(row)

    full_ids = {
        row.get("ConsensusID", "") for row in summary_rows
        if row.get("FullConsensusCalled", "").upper() == "YES"
    }
    all_records: list[tuple[str, str]] = []
    full_records: list[tuple[str, str]] = []
    for fasta in sorted(parent_dir.glob("*/mapping/rotavirus_mapping_consensus.fasta")):
        for header, seq in parse_fasta(fasta):
            key = fasta_key(header)
            all_records.append((header, seq))
            if key in full_ids:
                full_records.append((header, seq))
    write_fasta(all_records, consensus_all)
    write_fasta(full_records, consensus_full)

    for output, pattern, fields in [
        (combined_summary, "*/mapping/rotavirus_mapping_summary.csv", mapping_summary_fields()),
        (combined_features, "*/mapping/rotavirus_mapping_consensus_features.csv", mapping_consensus_feature_fields()),
    ]:
        files = sorted(parent_dir.glob(pattern))
        with output.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=fields)
            writer.writeheader()
            for src in files:
                with src.open("r", newline="", encoding="utf-8", errors="replace") as handle:
                    reader = csv.DictReader(handle)
                    for row in reader:
                        writer.writerow({field: row.get(field, "") for field in fields})

    log(f"Generated {consensus_all.name} with {len(all_records)} mapping consensus sequence(s)", log_path)
    log(f"Generated {consensus_full.name} with {len(full_records)} full mapping consensus sequence(s)", log_path)
    log(f"Generated {combined_summary.name} and {combined_features.name}", log_path)


def sample_folders(parent_dir: Path) -> list[Path]:
    folders: list[Path] = []
    for path in sorted(parent_dir.iterdir()):
        if not path.is_dir():
            continue
        if path.name.startswith("."):
            continue
        folders.append(path)
    return folders


def check_required_tools(
    skip_assembly: bool,
    assemblies: list[str],
    log_path: Path,
    mapping_enabled: bool = False,
    mapping_validation_enabled: bool = False,
) -> None:
    required = ["blastn"]
    if not skip_assembly:
        required.append("trimmomatic")
        required.append("repair.sh")
        if any(name in assemblies for name in ["rnaviralspades", "bbnorm100-rnaviralspades", "bbnorm500-rnaviralspades", "clumpify-rnaviralspades"]):
            required.append("rnaviralspades.py")
        if any(name in assemblies for name in ["bbnorm100-rnaviralspades", "bbnorm500-rnaviralspades"]):
            required.append("bbnorm.sh")
            required.append("reformat.sh")
        if "clumpify-rnaviralspades" in assemblies:
            required.append("clumpify.sh")
        if "spades-isolate" in assemblies:
            required.append("spades.py")

    if mapping_enabled:
        required.extend(["bwa", "samtools", "bcftools"])
    if mapping_validation_enabled:
        required.append("mafft")

    missing = sorted({tool for tool in required if shutil.which(tool) is None})
    if missing:
        log("Missing required executable(s) in PATH: " + ", ".join(missing), log_path)
        raise SystemExit(2)


def resolve_blast_db(blast_db_arg: str | None) -> str:
    if blast_db_arg:
        return str(Path(blast_db_arg).expanduser())

    db_dir = Path(__file__).resolve().parent / "blastdb"
    if not db_dir.exists():
        raise SystemExit(
            f"No --blast-db was provided and default BLAST database folder was not found: {db_dir}"
        )

    if db_dir.is_file():
        return str(db_dir)

    # BLAST database files normally share a prefix and end in extensions such as
    # .nhr/.nin/.nsq or .ndb/.not/.ntf/.nto. Auto-detect one unique prefix.
    prefixes: set[Path] = set()
    for path in db_dir.iterdir():
        name = path.name
        if name.endswith((".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto")):
            prefixes.add(path.with_suffix(""))
        elif re.search(r"\.\d{2}\.(nhr|nin|nsq)$", name):
            prefixes.add(Path(str(path).rsplit(".", 2)[0]))

    if len(prefixes) == 1:
        return str(next(iter(prefixes)))

    if len(prefixes) > 1:
        options = "\n  ".join(str(p) for p in sorted(prefixes))
        raise SystemExit(
            "More than one BLAST database prefix was found in the default blastdb folder. "
            "Please choose one with --blast-db. Found:\n  " + options
        )

    return str(db_dir)



def resolve_optional_blast_db(blast_db_arg: str | None, default_folder_name: str) -> str:
    if blast_db_arg:
        return str(Path(blast_db_arg).expanduser())

    db_dir = Path(__file__).resolve().parent / default_folder_name
    if not db_dir.exists():
        return ""
    if db_dir.is_file():
        return str(db_dir)

    prefixes: set[Path] = set()
    for path in db_dir.iterdir():
        name = path.name
        if name.endswith((".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto")):
            prefixes.add(path.with_suffix(""))
        elif re.search(r"\.\d{2}\.(nhr|nin|nsq)$", name):
            prefixes.add(Path(str(path).rsplit(".", 2)[0]))

    if len(prefixes) == 1:
        return str(next(iter(prefixes)))
    if len(prefixes) > 1:
        options = "\n  ".join(str(p) for p in sorted(prefixes))
        raise SystemExit(
            f"More than one BLAST database prefix was found in {db_dir}. "
            "Please choose one with --full-segment-db. Found:\n  " + options
        )
    return ""


def parse_args() -> argparse.Namespace:
    assembly_choices = list(ASSEMBLY_OPTIONS)

    parser = argparse.ArgumentParser(
        description="Single-file RotaX pipeline for rotavirus genotyping from Illumina reads or FASTA sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {ROTAX_VERSION} - {ROTAX_UPDATE_NOTE}",
    )
    parser.add_argument(
        "--blast-db",
        default=None,
        help=(
            "Path/prefix of the rotavirus ORF BLAST database used for genotyping. "
            "Default: auto-detect one BLAST database prefix inside a folder named blastdb next to this script."
        ),
    )
    parser.add_argument(
        "--full-segment-db",
        default=None,
        help=(
            "Optional path/prefix of a full-segment nucleotide BLAST database used only for GenBank/full-segment evaluation. "
            "Default: auto-detect one BLAST database prefix inside fulldb next to this script. "
            "If not found, GenBank reports are still created but marked NOT_EVALUATED_NO_FULL_SEGMENT_DB."
        ),
    )
    parser.add_argument(
        "--genbank-ref-cov-full",
        type=float,
        default=99.0,
        help="Minimum percent of the full-segment reference that must be covered to call GENBANK_READY_FULL.",
    )
    parser.add_argument(
        "--genbank-max-missing-end-bases",
        type=int,
        default=5,
        help="Maximum missing reference bases allowed at each end of the best same-segment full-reference hit.",
    )
    parser.add_argument(
        "--genbank-max-extra-end-bases",
        type=int,
        default=25,
        help="Maximum extra contig bases allowed at each end after terminal artefact trimming.",
    )
    parser.add_argument(
        "--genbank-terminal-homopolymer",
        type=int,
        default=12,
        help="Trim terminal homopolymer runs of at least this length before writing the clean/review FASTA.",
    )
    parser.add_argument(
        "--genbank-terminal-motif-mismatches",
        type=int,
        default=1,
        help=(
            "Maximum mismatches allowed in conserved RVA terminal motifs. "
            "5-prime motif is degenerate GGCWWTWWWW. 3-prime motif is degenerate WTRTGRCC. "
            "For NSP1, either WTRTGRCC or the special 9-mer CTGTGAACC is accepted."
        ),
    )
    parser.add_argument(
        "--genbank-orf-length-tolerance",
        type=int,
        default=90,
        help=(
            "Allowed absolute difference from the segment's expected ORF length before a valid ORF is marked "
            "FULL_BUT_NEEDS_REVIEW. The difference is still reported even within tolerance."
        ),
    )
    parser.add_argument(
        "--spades-full-segment-qc",
        action="store_true",
        help=(
            "Optionally generate the older SPAdes-contig full-segment/terminal QC reports. "
            "Disabled by default because SPAdes contigs are used for genotyping and ORF-guided reference selection, "
            "not as the final GenBank sequences."
        ),
    )
    parser.add_argument(
        "--no-genbank-alignments",
        dest="genbank_make_alignments",
        action="store_false",
        help="When --spades-full-segment-qc is enabled, disable its per-candidate alignment FASTA files.",
    )
    parser.add_argument(
        "--genbank-alignment-full-ref-count",
        type=int,
        default=2,
        help="Number of closest same-segment full-reference sequences to include in each GenBank QC alignment.",
    )
    parser.add_argument(
        "--genbank-alignment-cmd",
        default="mafft",
        help="Alignment command for GenBank QC FASTA files. Use 'none' to write only unaligned FASTA files.",
    )
    parser.add_argument(
        "--genbank-terminal-motif-mode",
        choices=["strict", "warn", "off"],
        default="strict",
        help=(
            "How terminal motif failures affect GenBank classification: strict=call PARTIAL_SEGMENT, "
            "warn=call FULL_BUT_NEEDS_REVIEW if otherwise ready, off=report only."
        ),
    )
    parser.add_argument(
        "--genbank-allow-ambiguity",
        action="store_true",
        help="Allow ambiguous/non-ACGT bases in sequences marked GENBANK_READY_FULL.",
    )
    parser.add_argument(
        "--no-mapping",
        dest="mapping_consensus",
        action="store_false",
        help=(
            "Disable the read-mapping consensus step. By default, reads are mapped once to a shared BWA index "
            "containing one closest same-segment full reference per selected ORF candidate."
        ),
    )
    parser.add_argument(
        "--mapping-r1",
        default=None,
        help=(
            "Optional R1 FASTQ(.gz) file to use for the mapping-consensus step. "
            "This is mainly for --fasta-input mode, where the script otherwise has contigs but no reads."
        ),
    )
    parser.add_argument(
        "--mapping-r2",
        default=None,
        help=(
            "Optional R2 FASTQ(.gz) file to use for the mapping-consensus step. "
            "Must be supplied together with --mapping-r1."
        ),
    )
    parser.add_argument(
        "--allow-raw-fastq-mapping",
        action="store_true",
        help=(
            "Allow the mapping-consensus step to fall back to the raw, untrimmed "
            "*_R1_001.fastq.gz / *_R2_001.fastq.gz reads when no trimmed "
            "(R1_paired.fastq.gz) or deduped (deduped_R1.fastq.gz) reads are present. "
            "Off by default: mapping refuses rather than silently using raw reads, so "
            "a deleted set of trimmed reads is surfaced as an error instead of quietly "
            "producing a consensus from unfiltered reads."
        ),
    )
    parser.add_argument(
        "--mapping-min-depth",
        type=int,
        default=5,
        help="Minimum read depth required for a mapping consensus base; lower-depth positions are masked to N.",
    )
    parser.add_argument(
        "--mapping-polish-rounds",
        type=int,
        default=1,
        help=(
            "Number of iterative remap-to-consensus polish rounds after the initial mapping. "
            "Each round rebuilds a SNP-corrected, coordinate-stable mapping target (indels excluded, "
            "low-depth positions kept as the reference base) and remaps the reads to it so read ends that "
            "soft-clipped against the divergent GenBank reference realign and extend terminal coverage. "
            "Set to 0 to disable and reproduce the single-pass mapping behaviour."
        ),
    )
    parser.add_argument(
        "--mapping-ref-cov-full",
        type=float,
        default=99.0,
        help="Minimum percent of the selected full reference covered at mapping-min-depth to classify a mapping consensus as full.",
    )
    parser.add_argument(
        "--mapping-max-missing-end-bases",
        type=int,
        default=5,
        help="Maximum uncovered bases allowed at each reference end to classify a mapping consensus as full.",
    )
    parser.add_argument(
        "--mapping-orf-alignment-cmd",
        default="mafft",
        help=(
            "Protein alignment command used both to determine mapping-consensus CDS coordinates from the best "
            "genotyping-DB ORF hit and to compare SPAdes versus mapping ORFs. MAFFT is the default. Setting this "
            "to 'none' disables protein-guided CDS annotation; the source feature note will report that no CDS was annotated."
        ),
    )
    parser.add_argument(
        "--mapping-full-ref-count",
        type=int,
        default=1,
        help=(
            "Compatibility option. RotaX v9 mapping consensus always indexes exactly one closest same-segment "
            "full-reference sequence per selected ORF/SPAdes candidate. Values above 1 are ignored for mapping; "
            "use --genbank-alignment-full-ref-count for QC alignments with extra references."
        ),
    )
    parser.add_argument(
        "--no-mapping-validation",
        dest="mapping_validation",
        action="store_false",
        help=(
            "Do not generate mapping/rotavirus_mapping_validation.html. By default "
            "RotaX aligns each mapping consensus to its selected mapping reference "
            "and writes one combined metadata, ORF, depth-logo, and base-alignment HTML."
        ),
    )
    parser.add_argument(
        "--mapping-validation-email",
        default=os.environ.get("NCBI_EMAIL", ""),
        help=(
            "Email supplied to NCBI Entrez when downloading GenBank CDS annotations "
            "for the integrated mapping-validation HTML. Defaults to NCBI_EMAIL. "
            "Without it, the HTML still uses the exact local mapping reference but "
            "reports the reference CDS as unavailable."
        ),
    )
    parser.add_argument(
        "--mapping-validation-api-key",
        default=os.environ.get("NCBI_API_KEY", ""),
        help=(
            "Optional NCBI API key for mapping-validation GenBank retrieval; "
            "defaults to NCBI_API_KEY."
        ),
    )
    parser.add_argument(
        "--mapping-validation-depth-scale",
        choices=["log", "sqrt", "linear"],
        default="log",
        help=(
            "Depth transform for nucleotide-letter height in the integrated "
            "mapping-validation HTML; each segment is capped at p99."
        ),
    )
    parser.add_argument(
        "--allow-mixed-genbank",
        action="store_true",
        help=(
            "Allow samples with Extra Information / possible mixed infection signals to be included in the GenBank-ready "
            "mapping full-consensus FASTA. Default is conservative: mixed samples are mapped and reported but kept in review outputs."
        ),
    )
    parser.add_argument(
        "--fasta-input",
        default=None,
        help=(
            "FASTA file to genotype directly without assembly. A sample folder named after the FASTA file stem "
            "is created in the current working directory."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=12,
        help="Threads for SPAdes/rnaviralSPAdes.",
    )
    parser.add_argument(
        "--bbnorm-threads",
        type=int,
        default=1,
        help=(
            "Threads for BBTools bbnorm.sh. RotaX interleaves paired input to avoid "
            "the BBMap 39.81 dual-stream PairStreamer assertion. Default 1 remains "
            "conservative; increase for large datasets after testing your BBTools version."
        ),
    )
    parser.add_argument(
        "--min-cov",
        type=float,
        default=3.0,
        help=(
            "Minimum SPAdes k-mer coverage parsed from contig headers. "
            "This is not applied in --fasta-input mode."
        ),
    )
    parser.add_argument(
        "--assemblies",
        nargs="+",
        choices=assembly_choices,
        default=DEFAULT_ASSEMBLIES,
        metavar="MODE",
        help=(
            "Assembly mode(s) to run. Default is all modes. "
            "Choices: " + ", ".join(assembly_choices)
        ),
    )
    parser.add_argument(
        "--no-fastq-repair",
        dest="fastq_repair",
        action="store_false",
        help=(
            "Disable repair.sh synchronization of Trimmomatic paired FASTQs before BBnorm/SPAdes. "
            "Default is to repair/synchronise pairs to avoid BBTools PairStreamer list-size mismatch errors."
        ),
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
        help=(
            "Keep intermediate work folders, FASTQ files, mapping BAM/VCFs, "
            "polished references, raw depth, mask BED, and raw consensus."
        ),
    )
    existing_group = parser.add_mutually_exclusive_group()
    existing_group.add_argument(
        "--overwrite-existing",
        action="store_true",
        help=(
            "Do not ask interactively. If existing RotaX outputs are found in a sample folder, "
            "delete those RotaX outputs and rerun the sample."
        ),
    )
    existing_group.add_argument(
        "--skip-existing",
        action="store_true",
        help=(
            "Do not ask interactively. If existing RotaX outputs are found in a sample folder, "
            "skip that folder."
        ),
    )
    parser.set_defaults(
        cleanup=True,
        genbank_make_alignments=True,
        mapping_consensus=True,
        mapping_validation=True,
        fastq_repair=True,
    )
    args = parser.parse_args()

    # FASTA input is a direct genotyping mode and should not try to assemble.
    if args.fasta_input:
        args.skip_assembly = True

    if bool(args.mapping_r1) != bool(args.mapping_r2):
        parser.error("--mapping-r1 and --mapping-r2 must be supplied together")
    if args.mapping_r1:
        args.mapping_r1 = str(Path(args.mapping_r1).expanduser())
        args.mapping_r2 = str(Path(args.mapping_r2).expanduser())

    # Deduplicate while preserving input order.
    seen: set[str] = set()
    args.assemblies = [mode for mode in args.assemblies if not (mode in seen or seen.add(mode))]

    args.blast_db = resolve_blast_db(args.blast_db)
    args.full_segment_db = resolve_optional_blast_db(args.full_segment_db, "fulldb")
    return args


def main() -> int:
    args = parse_args()

    parent_dir = Path.cwd()
    log_path = parent_dir / "rotax_pipeline_log.txt"
    start_time = time.time()
    start_date = time.strftime("%Y-%m-%d %H:%M:%S")

    log(f"Script started in {parent_dir}", log_path)
    log(f"RotaX script version: {ROTAX_VERSION}", log_path)
    log(f"RotaX update note: {ROTAX_UPDATE_NOTE}", log_path)
    log(f"Start time: {start_date}", log_path)
    log(f"Rotavirus ORF BLAST DB: {args.blast_db}", log_path)
    if args.full_segment_db:
        log(f"Rotavirus full-segment BLAST DB: {args.full_segment_db}", log_path)
    else:
        log("Rotavirus full-segment BLAST DB: not found/provided; ORF-guided mapping consensus will be unavailable", log_path)
    if args.skip_assembly:
        log("Assembly step: skipped; existing contigs500.fasta will be reused in each sample folder", log_path)
    else:
        log("Selected assembly modes: " + ", ".join(args.assemblies), log_path)
        log(f"BBnorm threads: {args.bbnorm_threads}", log_path)
    log(f"Mapping consensus step: {'enabled' if args.mapping_consensus else 'disabled'}", log_path)
    log(
        "Integrated mapping-validation HTML: "
        + ("enabled" if args.mapping_validation else "disabled"),
        log_path,
    )
    if args.mapping_validation and args.mapping_validation_email:
        log("Mapping-validation GenBank CDS retrieval: enabled", log_path)
    elif args.mapping_validation:
        log(
            "Mapping-validation GenBank CDS retrieval: disabled (no email); "
            "local mapping references will still be shown",
            log_path,
        )
    log(f"Optional SPAdes full-segment QC: {'enabled' if args.spades_full_segment_qc else 'disabled'}", log_path)
    if args.fasta_input:
        log(f"FASTA input mode: {args.fasta_input}", log_path)
    if args.mapping_r1 and args.mapping_r2:
        log(f"Mapping read override R1: {args.mapping_r1}", log_path)
        log(f"Mapping read override R2: {args.mapping_r2}", log_path)

    mapping_tools_needed = bool(args.mapping_consensus and args.full_segment_db and (not args.fasta_input or (args.mapping_r1 and args.mapping_r2)))
    validation_tools_needed = bool(mapping_tools_needed and args.mapping_validation)
    check_required_tools(
        args.skip_assembly,
        args.assemblies,
        log_path,
        mapping_tools_needed,
        validation_tools_needed,
    )

    if args.fasta_input:
        fasta_input_path = Path(args.fasta_input).expanduser()
        planned_sample_dir = parent_dir / fasta_input_path.stem

        if planned_sample_dir.exists() and existing_rotax_result_paths(planned_sample_dir):
            selected = apply_existing_results_choice([planned_sample_dir], args, log_path)
            if not selected:
                log(f"FASTA-input sample folder {planned_sample_dir.name} was skipped because existing results were found", log_path)
                folders = []
            else:
                fasta_sample = prepare_fasta_input_sample(fasta_input_path, parent_dir, log_path)
                folders = [fasta_sample] if fasta_sample else []
        else:
            fasta_sample = prepare_fasta_input_sample(fasta_input_path, parent_dir, log_path)
            folders = [fasta_sample] if fasta_sample else []
    else:
        folders = sample_folders(parent_dir)
        folders = apply_existing_results_choice(folders, args, log_path)

    if not folders:
        log("No sample folders selected for processing", log_path)
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
    combine_low_identity_x_reports(parent_dir, log_path)
    if bool(getattr(args, "spades_full_segment_qc", False)):
        combine_genbank_outputs(parent_dir, log_path)
    else:
        for legacy_name in [
            "rotax_genbank_full_segment_report.csv",
            "rotax_genbank_clean_all_full_segments.fasta",
            "rotax_genbank_clean_complete_11seg_samples.fasta",
        ]:
            safe_remove(parent_dir / legacy_name, log_path)
    combine_mapping_outputs(parent_dir, log_path)

    elapsed = int(time.time() - start_time)
    log(f"Script ended in {parent_dir}", log_path)
    log(f"Time spent: {elapsed} seconds", log_path)
    log(f"Folders processed successfully: {processed}", log_path)
    log(f"Folders failed/skipped: {failed}", log_path)

    return 0 if processed > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
