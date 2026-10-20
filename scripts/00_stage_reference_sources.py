#!/usr/bin/env python3
"""Copy the local curated reference inputs into the project input area."""
from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT / "references" / "ugt84_expanded" / "sequences"
DESTINATION = PROJECT / "Inputs" / "reference_panel"
FILES = (
    "formal_ugt84_reference_seeds.faa",
    "formal_ugt84_reference_manifest.tsv",
    "groupL_ugt74_comparator_seeds.faa",
    "family1_context_outgroups.faa",
)
EXPECTED_OUTGROUPS = {"AAW56092.1", "AAR06912.1"}


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence)))
                header, sequence = line[1:], []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records


def accession(header: str) -> str:
    return header.split("|", 1)[0].split()[0]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    if not SOURCE.is_dir():
        raise FileNotFoundError(f"Included reference source data are missing: {SOURCE}")
    if DESTINATION.exists():
        raise FileExistsError(f"Refusing to overwrite staged inputs: {DESTINATION}")
    for name in FILES:
        if not (SOURCE / name).is_file():
            raise FileNotFoundError(SOURCE / name)

    formal_records = parse_fasta(SOURCE / FILES[0])
    with (SOURCE / FILES[1]).open(newline="") as handle:
        formal_rows = list(csv.DictReader(handle, delimiter="\t"))
    group_records = parse_fasta(SOURCE / FILES[2])
    context_records = parse_fasta(SOURCE / FILES[3])
    formal_ids = {accession(header) for header, _ in formal_records}
    manifest_ids = {row["protein_accession"] for row in formal_rows}
    if len(formal_records) != 24 or manifest_ids != formal_ids or len(formal_rows) != 24:
        raise ValueError("Formal UGT84 FASTA and manifest must contain the same 24 accessions")
    if len(group_records) != 6:
        raise ValueError("The local Group-L UGT74 source FASTA must contain six records")
    if {accession(header) for header, _ in context_records} != EXPECTED_OUTGROUPS:
        raise ValueError("The family-1 context FASTA must contain UGT71G1 and UGT76G1 accessions")

    DESTINATION.mkdir(parents=True, exist_ok=False)
    staged = []
    for name in FILES:
        source_path = SOURCE / name
        target_path = DESTINATION / name
        shutil.copy2(source_path, target_path)
        if sha256(source_path) != sha256(target_path):
            raise RuntimeError(f"Copy hash mismatch: {source_path} -> {target_path}")
        staged.append({
            "source": str(source_path.relative_to(PROJECT)),
            "staged_file": str(target_path.relative_to(PROJECT)),
            "size_bytes": target_path.stat().st_size,
            "sha256": sha256(target_path),
        })

    manifest = {
        "description": "Locally staged inputs for building the 32-sequence reference/context dataset",
        "counts": {"formal_ugt84": 24, "group_l_ugt74_comparators": 6, "family1_context": 2},
        "files": staged,
    }
    manifest_path = DESTINATION / "source_manifest.json"
    with manifest_path.open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    print(f"Staged and hash-verified {len(FILES)} local source files under {DESTINATION}")
    print("Record counts: 24 UGT84 references, six Group-L UGT74 comparators, two family-1 context sequences.")


if __name__ == "__main__":
    main()
