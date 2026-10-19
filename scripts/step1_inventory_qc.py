#!/usr/bin/env python3
"""Preliminary streaming FASTA QC for the local Cornus assemblies.

Uses only the Python standard library. Paths are relative to this script's
project root. The script writes results/assembly_qc.tsv, creates results/ only
after input preflight, and skips the run when that output already exists.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
METADATA = ROOT / "genome_metadata.csv"
RESULTS = ROOT / "results"
OUTPUT = RESULTS / "assembly_qc.tsv"
CANONICAL = set("ACGT")
IUPAC_DNA = set("ACGTRYSWKMBDHVN")


def fasta_metrics(path: Path) -> dict:
    lengths: list[int] = []
    counts: Counter[str] = Counter()
    current_length = 0
    records = 0
    invalid = Counter()

    with path.open("r", encoding="ascii", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if records:
                    lengths.append(current_length)
                records += 1
                current_length = 0
                continue
            if not line.strip():
                continue
            if records == 0:
                raise ValueError(f"Sequence data before first FASTA header in {path}:{line_number}")
            seq = line.strip().upper()
            current_length += len(seq)
            counts.update(seq)
            invalid.update(ch for ch in seq if ch not in IUPAC_DNA)

    if records:
        lengths.append(current_length)
    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    total = sum(lengths)
    running = 0
    n50 = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= math.ceil(total / 2):
            n50 = length
            break

    canonical_total = sum(counts[base] for base in CANONICAL)
    gc_percent = (
        100.0 * (counts["G"] + counts["C"]) / canonical_total
        if canonical_total
        else 0.0
    )
    return {
        "fasta_records": records,
        "total_bp": total,
        "n50_bp": n50,
        "gc_percent_acgt": f"{gc_percent:.4f}",
        "n_bases": counts["N"],
        "ambiguous_bases": sum(counts[ch] for ch in counts if ch not in CANONICAL),
        "invalid_symbols": ";".join(f"{char}:{count}" for char, count in sorted(invalid.items())),
        "sequence_lengths": sorted(lengths),
    }


def sequence_report_lengths(path: Path) -> list[int]:
    lengths = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
                lengths.append(int(row["length"]))
            except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"Invalid sequence report row {path}:{line_number}: {exc}") from exc
    return sorted(lengths)


def first_genbank_record_features(path: Path) -> str:
    """Return feature keys from the first GenBank record as a preliminary screen."""
    keys = set()
    in_features = False
    with path.open("r", encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN"):
                break
            if in_features and line.startswith("     ") and len(line) > 5:
                key = line[5:21].strip()
                if key and not key.startswith("/"):
                    keys.add(key)
    return ";".join(sorted(keys))


def present(root: Path, relative_path: str) -> str:
    if not relative_path:
        return ""
    return "yes" if (root / relative_path).is_file() else "MISSING"


def main() -> int:
    if OUTPUT.exists():
        print(f"Output already exists; skipping QC: {OUTPUT}")
        return 0
    if RESULTS.exists() and not RESULTS.is_dir():
        print(f"Results path exists but is not a directory: {RESULTS}", file=sys.stderr)
        return 2
    if not METADATA.is_file():
        print(f"Metadata file not found: {METADATA}", file=sys.stderr)
        return 2

    with METADATA.open(newline="", encoding="utf-8") as handle:
        assemblies = list(csv.DictReader(handle))
    if not assemblies:
        print("No assembly rows found in genome_metadata.csv", file=sys.stderr)
        return 2

    required_fields = ("species", "accession", "genome_fasta", "genbank_flatfile", "assembly_archive", "extracted_directory")
    errors = []
    for assembly in assemblies:
        accession = assembly.get("accession", "<unknown>")
        missing = [field for field in required_fields if not assembly.get(field)]
        if missing:
            errors.append(f"{accession}: blank metadata fields: {', '.join(missing)}")
            continue
        required_paths = (
            assembly["genome_fasta"], assembly["genbank_flatfile"], assembly["assembly_archive"],
            f"{assembly['extracted_directory']}/sequence_report.jsonl",
        )
        absent = [relative for relative in required_paths if not (ROOT / relative).is_file()]
        if absent:
            errors.append(f"{accession}: missing input files: {', '.join(absent)}")
    if errors:
        print("Input metadata/files failed preflight; no QC output was created:", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 2

    RESULTS.mkdir(parents=True, exist_ok=True)

    columns = [
        "species", "tax_id", "accession", "release_date_as_reported", "release_date_note",
        "source_archive", "archive_present", "genome_fasta", "genbank_flatfile",
        "separate_annotation_files", "first_gbff_record_feature_keys",
        "preliminary_annotation_assessment", "sequence_report_present", "fasta_records",
        "sequence_report_records", "total_bp", "sequence_report_total_bp", "n50_bp",
        "gc_percent_acgt", "n_bases", "ambiguous_bases", "invalid_symbols",
        "sequence_report_matches_fasta",
    ]
    result_rows = []

    for assembly in assemblies:
        fasta_path = assembly["genome_fasta"]
        gbff_path = assembly["genbank_flatfile"]
        archive_path = assembly["assembly_archive"]
        species = assembly["species"]
        accession = assembly["accession"]
        if not fasta_path or not (ROOT / fasta_path).is_file():
            raise FileNotFoundError(f"Genome FASTA missing for {accession}: {fasta_path}")
        if not gbff_path or not (ROOT / gbff_path).is_file():
            raise FileNotFoundError(f"GenBank flatfile missing for {accession}: {gbff_path}")

        metrics = fasta_metrics(ROOT / fasta_path)
        report_path = ROOT / assembly["extracted_directory"] / "sequence_report.jsonl"
        report_lengths = sequence_report_lengths(report_path) if report_path.is_file() else []
        report_total = sum(report_lengths) if report_lengths else ""
        matches = bool(report_lengths) and report_lengths == metrics.pop("sequence_lengths")

        feature_keys = first_genbank_record_features(ROOT / gbff_path)
        separate_annotations = [
            (label, assembly[field])
            for label, field in (
                ("GFF", "gff"), ("GTF", "gtf"), ("CDS", "cds_fasta"),
                ("protein", "protein_fasta"), ("RNA", "rna_fasta"),
            )
            if assembly.get(field)
        ]
        present_annotations = [label for label, path in separate_annotations if present(ROOT, path) == "yes"]
        if present_annotations:
            annotation_assessment = "separate annotation bundle present: " + ",".join(present_annotations)
        elif feature_keys:
            annotation_assessment = "no separate bundle; features found in first GBFF record (screen only)"
        else:
            annotation_assessment = "no separate bundle; no features in first GBFF record (screen only)"

        release_date = assembly.get("release_date", "")
        release_note = "verify source report" if release_date.startswith("2026-") else ""
        result_rows.append({
            "species": species,
            "tax_id": assembly.get("tax_id", ""),
            "accession": accession,
            "release_date_as_reported": release_date,
            "release_date_note": release_note,
            "source_archive": archive_path,
            "archive_present": present(ROOT, archive_path),
            "genome_fasta": fasta_path,
            "genbank_flatfile": gbff_path,
            "separate_annotation_files": ";".join(present_annotations),
            "first_gbff_record_feature_keys": feature_keys,
            "preliminary_annotation_assessment": annotation_assessment,
            "sequence_report_present": "yes" if report_lengths else "no",
            "fasta_records": metrics["fasta_records"],
            "sequence_report_records": len(report_lengths) if report_lengths else "",
            "total_bp": metrics["total_bp"],
            "sequence_report_total_bp": report_total,
            "n50_bp": metrics["n50_bp"],
            "gc_percent_acgt": metrics["gc_percent_acgt"],
            "n_bases": metrics["n_bases"],
            "ambiguous_bases": metrics["ambiguous_bases"],
            "invalid_symbols": metrics["invalid_symbols"],
            "sequence_report_matches_fasta": "yes" if matches else "no" if report_lengths else "not available",
        })
        print(f"QC complete: {species} ({accession})", flush=True)

    with OUTPUT.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result_rows)
    print(f"Wrote {OUTPUT.relative_to(ROOT)} ({len(result_rows)} assemblies)")
    print("BUSCO was not run by this script; see README.md for the lineage-data requirement.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
