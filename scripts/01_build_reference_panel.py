#!/usr/bin/env python3
"""Build the 32-sequence UGT84 reference/context FASTA and its data manifest."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Inputs" / "reference_panel"
OUTPUT = ROOT / "Results"
FASTA_OUT = OUTPUT / "ugt84_reference_context_panel_32.faa"
MANIFEST_OUT = OUTPUT / "ugt84_reference_context_panel_manifest.tsv"
GROUP_SOURCE = "Hoffmann et al. 2023, Table 2; PMID 37346124; doi:10.3389/fpls.2023.1191625"


def parse_fasta(path: Path) -> list[tuple[str, str]]:
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
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    return records


def accession(header: str) -> str:
    return header.split("|", 1)[0].split()[0]


def parse_pipe_header(header: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for field in header.split("|")[1:]:
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
    return values


def safe(value: str) -> str:
    return value.replace(" ", "_").replace("|", "_")


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("x") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def main() -> None:
    if not SOURCE.is_dir():
        raise FileNotFoundError(f"Run 00_stage_reference_sources.py first: {SOURCE}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if FASTA_OUT.exists() or MANIFEST_OUT.exists():
        raise FileExistsError("Reference-panel output exists; refusing to overwrite")

    manifest_path = SOURCE / "formal_ugt84_reference_manifest.tsv"
    with manifest_path.open(newline="") as handle:
        formal_meta = {row["protein_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
    formal_records = parse_fasta(SOURCE / "formal_ugt84_reference_seeds.faa")
    if len(formal_records) != 24 or set(formal_meta) != {accession(h) for h, _ in formal_records}:
        raise ValueError("Formal UGT84 sequence and metadata sets differ or are not 24 records")

    rows: list[dict[str, str]] = []
    final_records: list[tuple[str, str]] = []
    fields = [
        "protein_accession", "taxon", "gene_name", "family", "subfamily",
        "phylogenetic_group", "panel_role", "sequence_length_aa",
        "membership_evidence", "function_evidence", "sequence_source_header",
        "source_file", "group_classification_source",
    ]

    for original_header, sequence in formal_records:
        acc = accession(original_header)
        meta = formal_meta[acc]
        taxon = meta["taxon"]
        gene = meta["gene_name"]
        # PMC11921145 Table 1 gives the exact alias as VvGT1, not the source-manifest typo VvgGT1.
        if acc == "AEW31187.1":
            gene = "VvGT1 / UGT84A44"
        subfamily = meta["subfamily"]
        fasta_header = f"{acc}|taxon={safe(taxon)}|gene={safe(gene)}|family=UGT84|subfamily={subfamily}|group=L|role=UGT84_reference"
        final_records.append((fasta_header, sequence))
        rows.append({
            "protein_accession": acc,
            "taxon": taxon,
            "gene_name": gene,
            "family": "UGT84",
            "subfamily": subfamily,
            "phylogenetic_group": "L",
            "panel_role": "UGT84_reference",
            "sequence_length_aa": str(len(sequence)),
            "membership_evidence": meta["membership_evidence"],
            "function_evidence": meta["function_evidence"],
            "sequence_source_header": original_header,
            "source_file": "References/ugt84_expanded/sequences/formal_ugt84_reference_seeds.faa",
            "group_classification_source": GROUP_SOURCE,
        })

    group_records = parse_fasta(SOURCE / "groupL_ugt74_comparator_seeds.faa")
    if len(group_records) != 6:
        raise ValueError(f"Expected six UGT74 Group-L comparators, found {len(group_records)}")
    for original_header, sequence in group_records:
        acc = accession(original_header)
        meta = parse_pipe_header(original_header)
        gene = meta.get("gene", "")
        taxon = meta.get("taxon", "").replace("_", " ")
        if not acc or "UGT74" not in gene:
            raise ValueError(f"Unexpected Group-L comparator header: {original_header}")
        fasta_header = f"{acc}|taxon={safe(taxon)}|gene={safe(gene)}|family=UGT74|group=L|role=Group_L_UGT74_comparator"
        final_records.append((fasta_header, sequence))
        rows.append({
            "protein_accession": acc,
            "taxon": taxon,
            "gene_name": gene,
            "family": "UGT74",
            "subfamily": gene,
            "phylogenetic_group": "L",
            "panel_role": "Group_L_UGT74_comparator",
            "sequence_length_aa": str(len(sequence)),
            "membership_evidence": f"NCBI protein accession {acc}; gene label in saved protein record",
            "function_evidence": "Not used as an experimentally demonstrated UGT84 function reference",
            "sequence_source_header": original_header,
            "source_file": "References/ugt84_expanded/sequences/groupL_ugt74_comparator_seeds.faa",
            "group_classification_source": GROUP_SOURCE,
        })

    context_info = {
        "AAW56092.1": ("Medicago truncatula", "UGT71G1", "UGT71", "E"),
        "AAR06912.1": ("Stevia rebaudiana", "UGT76G1", "UGT76", "H"),
    }
    context_records = parse_fasta(SOURCE / "family1_context_outgroups.faa")
    if len(context_records) != 2:
        raise ValueError(f"Expected two context sequences, found {len(context_records)}")
    for original_header, sequence in context_records:
        acc = accession(original_header)
        if acc not in context_info:
            raise ValueError(f"Unexpected context accession: {acc}")
        taxon, gene, family, group = context_info[acc]
        fasta_header = f"{acc}|taxon={safe(taxon)}|gene={gene}|family={family}|group={group}|role=family1_context_outgroup"
        final_records.append((fasta_header, sequence))
        rows.append({
            "protein_accession": acc,
            "taxon": taxon,
            "gene_name": gene,
            "family": family,
            "subfamily": gene,
            "phylogenetic_group": group,
            "panel_role": "family1_context_outgroup",
            "sequence_length_aa": str(len(sequence)),
            "membership_evidence": f"NCBI protein record {acc} product name",
            "function_evidence": "Context sequence only; not a functional UGT84 reference",
            "sequence_source_header": original_header,
            "source_file": "References/ugt84_expanded/sequences/family1_context_outgroups.faa",
            "group_classification_source": GROUP_SOURCE,
        })

    accessions = [row["protein_accession"] for row in rows]
    if len(rows) != 32 or len(set(accessions)) != 32:
        raise ValueError(f"Expected 32 unique sequences; found {len(rows)} rows and {len(set(accessions))} accessions")
    if {"AEW31188.1", "AEW31189.1"} & set(accessions):
        raise ValueError("Related grape VvGT2/VvGT3 must not be in this panel")

    write_fasta(FASTA_OUT, final_records)
    with MANIFEST_OUT.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(final_records)} sequences to {FASTA_OUT}")
    print(f"Wrote sequence-level provenance and classification to {MANIFEST_OUT}")


if __name__ == "__main__":
    main()
