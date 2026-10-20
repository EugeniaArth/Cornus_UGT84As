#!/usr/bin/env python3
"""Select broad candidates by explicit BLAST thresholds and emit reusable data files."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

BLAST_FIELDS = ["qseqid", "sseqid", "pident", "length", "qcovs", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


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
    parser.add_argument("--blast", required=True, type=Path, help="Raw BLASTP outfmt 6 file from script 05")
    parser.add_argument("--candidates", required=True, type=Path, help="Broad-hit protein FASTA from script 04")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--top-per-reference", type=int, default=10,
                        help="Select this many top hits per reference at the coverage cutoff (default: 10)")
    parser.add_argument("--min-identity", type=float, default=50.0,
                        help="Also select hits at or above this identity percentage (default: 50)")
    parser.add_argument("--min-query-coverage", type=float, default=70.0,
                        help="Minimum query coverage for either selection rule (default: 70)")
    parser.add_argument("--max-evalue", type=float, default=1e-5)
    args = parser.parse_args()
    if not args.blast.is_file() or not args.candidates.is_file():
        raise FileNotFoundError("--blast and --candidates must exist")
    if (args.top_per_reference < 1 or not 0 <= args.min_identity <= 100 or
            not 0 <= args.min_query_coverage <= 100 or args.max_evalue < 0):
        raise ValueError("Invalid selection threshold")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_out, table_out = outdir / "selected_candidate_proteins.faa", outdir / "selected_candidate_hits.tsv"
    if fasta_out.exists() or table_out.exists():
        raise FileExistsError(f"Refusing to overwrite outputs in {outdir}")

    hits_by_reference: dict[str, list[dict[str, str]]] = defaultdict(list)
    with args.blast.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            values = line.rstrip("\n").split("\t")
            if len(values) != len(BLAST_FIELDS):
                raise ValueError(f"Expected {len(BLAST_FIELDS)} BLAST columns, found {len(values)}")
            row = dict(zip(BLAST_FIELDS, values))
            if float(row["qcovs"]) < args.min_query_coverage or float(row["evalue"]) > args.max_evalue:
                continue
            row["qseqid"] = row["qseqid"].split("|", 1)[0]
            hits_by_reference[row["qseqid"]].append(row)

    best: dict[str, dict[str, str]] = {}
    reasons: dict[str, set[str]] = defaultdict(set)
    supporting_references: dict[str, set[str]] = defaultdict(set)

    def retain(row: dict[str, str], reason: str) -> None:
        protein_id = row["sseqid"]
        reasons[protein_id].add(reason)
        supporting_references[protein_id].add(row["qseqid"])
        previous = best.get(protein_id)
        if (previous is None or float(row["bitscore"]) > float(previous["bitscore"]) or
                (float(row["bitscore"]) == float(previous["bitscore"]) and row["qseqid"] < previous["qseqid"])):
            best[protein_id] = row

    for reference, hits in sorted(hits_by_reference.items()):
        ranked = sorted(hits, key=lambda row: (-float(row["bitscore"]), row["sseqid"]))
        for rank, row in enumerate(ranked[:args.top_per_reference], start=1):
            retain(row, f"top{rank}_for_{reference}")
        for row in hits:
            if float(row["pident"]) >= args.min_identity:
                retain(row, f"identity_ge_{args.min_identity:g}_pct_and_query_coverage_ge_{args.min_query_coverage:g}_pct")

    proteins = read_fasta(args.candidates)
    missing = sorted(set(best) - set(proteins))
    if missing:
        raise KeyError(f"Selected BLAST subjects missing from candidate FASTA: {missing[:10]}")
    with fasta_out.open("x") as handle:
        for protein_id in sorted(best):
            header, sequence = proteins[protein_id]
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    fields = ["candidate_protein_id", "best_reference_query", "percent_identity", "query_coverage_percent",
              "evalue", "bitscore", "selection_reason", "supporting_reference_queries",
              "selection_top_n", "selection_min_identity", "selection_min_query_coverage", "selection_max_evalue"]
    with table_out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for protein_id in sorted(best):
            row = best[protein_id]
            writer.writerow({
                "candidate_protein_id": protein_id,
                "best_reference_query": row["qseqid"],
                "percent_identity": row["pident"],
                "query_coverage_percent": row["qcovs"],
                "evalue": row["evalue"],
                "bitscore": row["bitscore"],
                "selection_reason": ";".join(sorted(reasons[protein_id])),
                "supporting_reference_queries": ";".join(sorted(supporting_references[protein_id])),
                "selection_top_n": args.top_per_reference,
                "selection_min_identity": args.min_identity,
                "selection_min_query_coverage": args.min_query_coverage,
                "selection_max_evalue": args.max_evalue,
            })
    print(f"Selected {len(best)} sequence candidates; membership/function remain unresolved.")
    print(f"Candidate FASTA: {fasta_out}; selection data: {table_out}")


if __name__ == "__main__":
    main()
