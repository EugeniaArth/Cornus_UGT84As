#!/usr/bin/env python3
"""Choose one longest protein per mapped locus, retaining all selection metadata."""
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
                        raise ValueError(f"Duplicate FASTA identifier: {key}")
                    records[key] = (header, "".join(parts))
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError("Sequence encountered before FASTA header")
                parts.append(line)
    if header is not None:
        key = header.split()[0]
        if key in records:
            raise ValueError(f"Duplicate FASTA identifier: {key}")
        records[key] = (header, "".join(parts))
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path, help="selected_candidate_proteins.faa")
    parser.add_argument("--mapping", required=True, type=Path, help="candidate_gene_map.tsv from script 07")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    proteins = read_fasta(args.candidates)
    if not args.mapping.is_file():
        raise FileNotFoundError(args.mapping)
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_out, table_out = outdir / "locus_representatives.faa", outdir / "locus_representatives.tsv"
    if fasta_out.exists() or table_out.exists():
        raise FileExistsError(f"Refusing to overwrite outputs in {outdir}")

    mapping: dict[str, dict[str, str]] = {}
    with args.mapping.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["protein_accession"]] = row
    groups: dict[str, list[str]] = {}
    for protein_id in proteins:
        row = mapping.get(protein_id, {})
        gene_ids = [x for x in row.get("gene_ids", "").split(";") if x]
        locus_tags = [x for x in row.get("locus_tags", "").split(";") if x]
        if len(gene_ids) == 1:
            locus_key = f"gene:{gene_ids[0]}"
        elif len(locus_tags) == 1:
            locus_key = f"locus_tag:{locus_tags[0]}"
        else:
            locus_key = f"unresolved:{protein_id}"
        groups.setdefault(locus_key, []).append(protein_id)

    selected: list[str] = []
    rows = []
    for locus_key in sorted(groups):
        protein_ids = groups[locus_key]
        representative = sorted(protein_ids, key=lambda pid: (-len(proteins[pid][1]), pid))[0]
        selected.append(representative)
        row = mapping.get(representative, {})
        for protein_id in sorted(protein_ids):
            rows.append({
                "locus_key": locus_key,
                "protein_accession": protein_id,
                "sequence_length_aa": len(proteins[protein_id][1]),
                "selected_representative": "yes" if protein_id == representative else "no",
                "selection_rule": "longest sequence per mapped gene/locus; accession lexical order breaks length ties",
                "gene_ids": row.get("gene_ids", "") if protein_id == representative else mapping.get(protein_id, {}).get("gene_ids", ""),
                "locus_tags": row.get("locus_tags", "") if protein_id == representative else mapping.get(protein_id, {}).get("locus_tags", ""),
                "mapping_status": mapping.get(protein_id, {}).get("mapping_status", "unmapped"),
            })
    with fasta_out.open("x") as handle:
        for protein_id in selected:
            header, sequence = proteins[protein_id]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    fields = ["locus_key", "protein_accession", "sequence_length_aa", "selected_representative", "selection_rule", "gene_ids", "locus_tags", "mapping_status"]
    with table_out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Selected {len(selected)} locus representatives from {len(proteins)} proteins.")
    print(f"Representative FASTA: {fasta_out}; selection data: {table_out}")


if __name__ == "__main__":
    main()
