#!/usr/bin/env python3
"""Extract HMMER protein hits from tblout into reusable FASTA and accession data."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


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
                    key = header.split()[0]
                    if key in records:
                        raise ValueError(f"Duplicate protein FASTA identifier: {key}")
                    records[key] = (header, "".join(parts))
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError("Protein sequence encountered before FASTA header")
                parts.append(line)
    if header is not None:
        key = header.split()[0]
        if key in records:
            raise ValueError(f"Duplicate protein FASTA identifier: {key}")
        records[key] = (header, "".join(parts))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tblout", required=True, type=Path)
    parser.add_argument("--proteins", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if not args.tblout.is_file() or not args.proteins.is_file():
        raise FileNotFoundError("--tblout and --proteins must exist")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_out = outdir / "pf00201_hit_proteins.faa"
    table_out = outdir / "pf00201_hit_accessions.tsv"
    if fasta_out.exists() or table_out.exists():
        raise FileExistsError(f"Refusing to overwrite outputs in {outdir}")

    hits: dict[str, dict[str, str]] = {}
    with args.tblout.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            columns = line.rstrip("\n").split(maxsplit=18)
            if len(columns) < 18:
                raise ValueError(f"Malformed HMMER tblout row: {line.rstrip()}")
            target = columns[0]
            description = columns[18] if len(columns) > 18 else ""
            row = {
                "protein_id": target,
                "target_accession": columns[1],
                "target_length": columns[2],
                "query_name": columns[3],
                "full_sequence_evalue": columns[6],
                "full_sequence_score": columns[7],
                "best_domain_i_evalue": columns[12],
                "description": description,
            }
            old = hits.get(target)
            if old is None or float(row["full_sequence_evalue"]) < float(old["full_sequence_evalue"]):
                hits[target] = row

    proteins = read_fasta(args.proteins)
    missing = sorted(set(hits) - set(proteins))
    if missing:
        raise KeyError(f"HMMER hit IDs absent from protein FASTA (first 10): {missing[:10]}")
    with fasta_out.open("x") as handle:
        for protein_id in sorted(hits):
            header, sequence = proteins[protein_id]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    fields = ["protein_id", "target_accession", "target_length", "query_name", "full_sequence_evalue",
              "full_sequence_score", "best_domain_i_evalue", "description"]
    with table_out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for protein_id in sorted(hits):
            writer.writerow(hits[protein_id])
    print(f"Extracted {len(hits)} broad PF00201 hit proteins to {fasta_out}")
    print(f"Accession-level hit data: {table_out}")


if __name__ == "__main__":
    main()
