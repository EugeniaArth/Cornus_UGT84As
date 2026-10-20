#!/usr/bin/env python3
"""Combine the 32-sequence panel with finalized Cornus locus representatives."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_FASTA = ROOT / "Results" / "ugt84_reference_context_panel_32.faa"
REFERENCE_MANIFEST = ROOT / "Results" / "ugt84_reference_context_panel_manifest.tsv"


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
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, type=Path, help="locus_representatives.faa from script 08")
    parser.add_argument("--selection", required=True, type=Path, help="locus_representatives.tsv from script 08")
    parser.add_argument("--taxon", required=True, help="Cornus taxon label, e.g. Cornus_florida")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if not REFERENCE_FASTA.is_file() or not REFERENCE_MANIFEST.is_file():
        raise FileNotFoundError("Build and validate the reference panel first")
    if not args.candidates.is_file() or not args.selection.is_file():
        raise FileNotFoundError("--candidates and --selection must exist")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    fasta_out, manifest_out = outdir / "final_phylogeny_input.faa", outdir / "final_phylogeny_input_manifest.tsv"
    if fasta_out.exists() or manifest_out.exists():
        raise FileExistsError(f"Refusing to overwrite outputs in {outdir}")

    refs = read_fasta(REFERENCE_FASTA)
    candidates = read_fasta(args.candidates)
    with REFERENCE_MANIFEST.open(newline="") as handle:
        ref_rows = list(csv.DictReader(handle, delimiter="\t"))
    with args.selection.open(newline="") as handle:
        selection_rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["selected_representative"] == "yes"]
    candidate_metadata = {row["protein_accession"]: row for row in selection_rows}
    if {header.split("|", 1)[0].split()[0] for header, _ in candidates} != set(candidate_metadata):
        raise ValueError("Candidate representative FASTA and selection table do not match")

    all_records = list(refs)
    all_rows = list(ref_rows)
    reference_ids = {row["protein_accession"] for row in ref_rows}
    fields = list(ref_rows[0])
    for header, sequence in candidates:
        acc = header.split("|", 1)[0].split()[0]
        if acc in reference_ids:
            raise ValueError(f"Candidate accession duplicates reference panel ID: {acc}")
        metadata = candidate_metadata[acc]
        fasta_taxon = args.taxon.replace(" ", "_")
        all_records.append((f"{acc}|taxon={fasta_taxon}|role=Cornus_candidate", sequence))
        all_rows.append({
            "protein_accession": acc,
            "taxon": args.taxon.replace("_", " "),
            "gene_name": metadata.get("locus_tags", "") or metadata.get("gene_ids", ""),
            "family": "unassigned",
            "subfamily": "unassigned",
            "phylogenetic_group": "",
            "panel_role": "Cornus_candidate",
            "sequence_length_aa": str(len(sequence)),
            "membership_evidence": "Candidate recovery only; UGT84 membership not established",
            "function_evidence": "No functional evidence inferred",
            "sequence_source_header": header,
            "source_file": "candidate locus-representative FASTA supplied to script 09",
            "group_classification_source": "",
        })

    accessions = [row["protein_accession"] for row in all_rows]
    if len(accessions) != len(set(accessions)):
        raise ValueError("Duplicate accessions in final phylogeny input")
    with fasta_out.open("x") as handle:
        for header, sequence in all_records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")
    with manifest_out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"Wrote {len(all_records)} final reference-plus-candidate sequences to {fasta_out}")
    print(f"Matching sequence metadata: {manifest_out}")


if __name__ == "__main__":
    main()
