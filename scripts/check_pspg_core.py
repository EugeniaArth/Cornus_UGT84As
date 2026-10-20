#!/usr/bin/env python3
"""Screen proteins for the exact HCGWNS PSPG-box core; this is not a full-box test."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError("Sequence encountered before FASTA header")
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    if not records:
        raise ValueError(f"No sequences in {path}")
    ids = [h.split()[0] for h, _ in records]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate FASTA identifiers")
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True, type=Path, help="Protein FASTA, e.g. PF00201 hit proteins")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if not args.proteins.is_file():
        raise FileNotFoundError(args.proteins)
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    table_path = outdir / "pspg_core_check.tsv"
    manifest_path = outdir / "pspg_core_check_manifest.json"
    if table_path.exists() or manifest_path.exists():
        raise FileExistsError(f"Refusing to overwrite PSPG-core outputs in {outdir}")

    records = read_fasta(args.proteins)
    motif = "HCGWNS"
    fields = ["protein_id", "sequence_length_aa", "exact_HCGWNS_core_present", "positions_1_based", "motif_interpretation"]
    matches = 0
    with table_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for header, sequence in records:
            protein_id = header.split()[0]
            positions = []
            start = 0
            while True:
                index = sequence.upper().find(motif, start)
                if index < 0:
                    break
                positions.append(index + 1)
                start = index + 1
            present = bool(positions)
            matches += int(present)
            writer.writerow({
                "protein_id": protein_id,
                "sequence_length_aa": len(sequence),
                "exact_HCGWNS_core_present": "yes" if present else "no",
                "positions_1_based": ";".join(map(str, positions)),
                "motif_interpretation": "Exact conserved-core screen only; not a complete PSPG-box validation",
            })
    manifest = {
        "input_fasta": str(args.proteins.resolve()),
        "protein_count": len(records),
        "motif": motif,
        "motif_match_count": matches,
        "interpretation": "Exact HCGWNS PSPG-box core screen only; not a full PSPG-box model and not a UGT84 or Group-L classifier",
        "data_output": str(table_path),
    }
    with manifest_path.open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    print(f"Exact HCGWNS core detected in {matches}/{len(records)} proteins.")
    print("This is a conserved-core screen, not a full PSPG-box or subfamily classification.")
    print(f"Per-protein data: {table_path}")


if __name__ == "__main__":
    main()
