#!/usr/bin/env python3
"""Validate the 32-sequence reference/context dataset without writing a report."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "Results" / "ugt84_reference_context_panel_32.faa"
MANIFEST = ROOT / "Results" / "ugt84_reference_context_panel_manifest.tsv"
SOURCE_MANIFEST = ROOT / "Inputs" / "reference_panel" / "source_manifest.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    accession = header.split("|", 1)[0].split()[0]
                    if accession in records:
                        raise ValueError(f"Duplicate FASTA accession: {accession}")
                    records[accession] = (header, "".join(parts))
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError("Sequence encountered before FASTA header")
                parts.append(line)
    if header is not None:
        accession = header.split("|", 1)[0].split()[0]
        if accession in records:
            raise ValueError(f"Duplicate FASTA accession: {accession}")
        records[accession] = (header, "".join(parts))
    return records


def main() -> None:
    if not SOURCE_MANIFEST.is_file():
        raise FileNotFoundError("Run 00_stage_reference_sources.py first")
    source_manifest = json.loads(SOURCE_MANIFEST.read_text())
    if source_manifest.get("counts") != {"formal_ugt84": 24, "group_l_ugt74_comparators": 6, "family1_context": 2}:
        raise ValueError("Unexpected staged reference-source counts")
    for item in source_manifest.get("files", []):
        staged = ROOT / item["staged_file"]
        if not staged.is_file() or sha256(staged) != item["sha256"]:
            raise ValueError(f"Staged source file missing or hash mismatch: {staged}")
    if len(source_manifest.get("files", [])) != 4:
        raise ValueError("Expected four staged reference source files")
    if not FASTA.is_file() or not MANIFEST.is_file():
        raise FileNotFoundError("Run 01_build_reference_panel.py first")
    records = read_fasta(FASTA)
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_accession = {row["protein_accession"]: row for row in rows}
    if len(rows) != 32 or len(by_accession) != 32:
        raise ValueError(f"Expected 32 unique manifest rows; got {len(rows)} rows, {len(by_accession)} accessions")
    if set(records) != set(by_accession):
        raise ValueError("FASTA accessions and manifest accessions differ")

    role_counts = Counter(row["panel_role"] for row in rows)
    expected_roles = {
        "UGT84_reference": 24,
        "Group_L_UGT74_comparator": 6,
        "family1_context_outgroup": 2,
    }
    if role_counts != expected_roles:
        raise ValueError(f"Unexpected panel-role counts: {role_counts}")
    group_counts = Counter(row["phylogenetic_group"] for row in rows)
    if group_counts != Counter({"L": 30, "E": 1, "H": 1}):
        raise ValueError(f"Unexpected phylogenetic-group counts: {group_counts}")
    if {"AEW31188.1", "AEW31189.1"} & set(records):
        raise ValueError("VvGT2/VvGT3 are excluded from this panel but were found")

    for accession, row in by_accession.items():
        sequence = records[accession][1]
        if len(sequence) != int(row["sequence_length_aa"]):
            raise ValueError(f"Length mismatch for {accession}")
        if not sequence or any(char not in "ABCDEFGHIKLMNPQRSTVWXYZUO*-" for char in sequence.upper()):
            raise ValueError(f"Invalid or empty amino-acid sequence: {accession}")

    print("PASS: 32 unique protein sequences match 32 manifest records.")
    print("PASS: 24 UGT84 references, 6 UGT74 Group-L comparators, 2 family-1 context sequences.")
    print("PASS: group labels are L=30, E=1 (UGT71G1), H=1 (UGT76G1); VvGT2/3 are absent.")
    print("PASS: sequence lengths and amino-acid characters agree with the manifest.")


if __name__ == "__main__":
    main()
