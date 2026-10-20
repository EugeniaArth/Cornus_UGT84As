#!/usr/bin/env python3
"""Combine the curated panel with mapped *C. florida* and *C. suecica* representatives."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
REFERENCE_FASTA = PROJECT / "Results" / "ugt84_reference_context_panel_32.faa"
REFERENCE_MANIFEST = PROJECT / "Results" / "ugt84_reference_context_panel_manifest.tsv"
CANDIDATE_SETS = (
    {
        "taxon": "Cornus florida",
        "assembly": "GCF_030987335.1",
        "fasta": PROJECT / "Results/cornus_florida/locus_representatives/locus_representatives.faa",
        "loci": PROJECT / "Results/cornus_florida/locus_representatives/locus_representatives.tsv",
        "selection": PROJECT / "Results/cornus_florida/candidate_selection/selected_candidate_hits.tsv",
        "core": PROJECT / "Results/cornus_florida/pspg_core/pspg_core_check.tsv",
        "quality": None,
    },
    {
        "taxon": "Cornus suecica",
        "assembly": "GCA_963922005.1",
        "fasta": PROJECT / "Results/cornus_suecica/locus_representatives_canonical_ids/locus_representatives.faa",
        "loci": PROJECT / "Results/cornus_suecica/locus_representatives_canonical_ids/locus_representatives.tsv",
        "selection": PROJECT / "Results/cornus_suecica/candidate_selection_canonical_ids/selected_candidate_hits.tsv",
        "core": PROJECT / "Results/cornus_suecica/pspg_core_canonical_ids/pspg_core_check.tsv",
        "quality": PROJECT / "Results/cornus_suecica/proteome_validated_v3/suecica_protein_transcript_manifest.tsv",
    },
)
DEFAULT_OUTDIR = PROJECT / "Results" / "cornus_comparison" / "final_input"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="ascii", errors="replace") as handle:
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
                    raise ValueError(f"Sequence before header in {path}")
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))
    ids = [header.split()[0] for header, _ in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate FASTA identifiers in {path}")
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR,
                        help="New output directory; existing outputs are never overwritten")
    args = parser.parse_args()
    required = [REFERENCE_FASTA, REFERENCE_MANIFEST]
    for item in CANDIDATE_SETS:
        required.extend([item["fasta"], item["loci"], item["selection"], item["core"]])
        if item["quality"] is not None:
            required.append(item["quality"])
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing combined-input source files: {missing}")
    outdir = args.outdir.resolve()
    fasta_out = outdir / "combined_cornus_reference_candidates_65.faa"
    manifest_out = outdir / "combined_cornus_reference_candidates_65_manifest.tsv"
    if outdir.exists() or fasta_out.exists() or manifest_out.exists():
        raise FileExistsError(f"Refusing to overwrite combined-input outputs under {outdir}")

    reference_records = read_fasta(REFERENCE_FASTA)
    reference_rows = read_tsv(REFERENCE_MANIFEST)
    ref_by_accession = {row["protein_accession"]: row for row in reference_rows}
    if len(reference_records) != 32 or len(ref_by_accession) != 32:
        raise ValueError("Expected exactly 32 reference/context proteins")
    reference_ids = {header.split("|", 1)[0].split()[0] for header, _ in reference_records}
    if reference_ids != set(ref_by_accession):
        raise ValueError("Reference FASTA IDs and manifest accessions differ")

    final_records: list[tuple[str, str]] = []
    final_rows: list[dict[str, str]] = []
    tree_tips: set[str] = set()
    protein_accessions: set[str] = set()
    fields = [
        "tree_tip", "protein_accession", "taxon", "assembly_accession", "gene_name",
        "family", "subfamily", "phylogenetic_group", "panel_role", "sequence_length_aa",
        "membership_evidence", "function_evidence", "sequence_source_header", "source_file",
        "group_classification_source", "locus_key", "mapping_status", "best_reference_query",
        "percent_identity", "query_coverage_percent", "evalue", "selection_reason",
        "HCGWNS_core", "model_quality_flags",
    ]

    for header, sequence in reference_records:
        tree_tip = header.split()[0]
        accession = tree_tip.split("|", 1)[0]
        row = ref_by_accession[accession]
        if tree_tip in tree_tips or accession in protein_accessions:
            raise ValueError(f"Duplicate reference identifier: {tree_tip}")
        tree_tips.add(tree_tip)
        protein_accessions.add(accession)
        if len(sequence) != int(row["sequence_length_aa"]):
            raise ValueError(f"Reference sequence length mismatch for {accession}")
        final_records.append((tree_tip, sequence))
        final_rows.append({
            "tree_tip": tree_tip,
            "protein_accession": accession,
            "taxon": row["taxon"],
            "assembly_accession": "",
            "gene_name": row["gene_name"],
            "family": row["family"],
            "subfamily": row["subfamily"],
            "phylogenetic_group": row["phylogenetic_group"],
            "panel_role": row["panel_role"],
            "sequence_length_aa": str(len(sequence)),
            "membership_evidence": row["membership_evidence"],
            "function_evidence": row["function_evidence"],
            "sequence_source_header": row["sequence_source_header"],
            "source_file": row["source_file"],
            "group_classification_source": row["group_classification_source"],
            "locus_key": "", "mapping_status": "reference",
            "best_reference_query": "", "percent_identity": "", "query_coverage_percent": "",
            "evalue": "", "selection_reason": "", "HCGWNS_core": "", "model_quality_flags": "",
        })

    for item in CANDIDATE_SETS:
        fasta_records = read_fasta(item["fasta"])
        locus_rows = [row for row in read_tsv(item["loci"]) if row["selected_representative"] == "yes"]
        selection_rows = {row["candidate_protein_id"]: row for row in read_tsv(item["selection"])}
        core_rows = {row["protein_id"]: row for row in read_tsv(item["core"])}
        quality_rows = {}
        if item["quality"] is not None:
            quality_rows = {row["protein_accession"]: row for row in read_tsv(item["quality"])}
        fasta_by_id = {header.split()[0]: (header, sequence) for header, sequence in fasta_records}
        if set(fasta_by_id) != {row["protein_accession"] for row in locus_rows}:
            raise ValueError(f"Representative FASTA and locus table disagree for {item['taxon']}")
        if len({row["locus_key"] for row in locus_rows}) != len(locus_rows):
            raise ValueError(f"More than one representative per locus for {item['taxon']}")
        for locus_row in locus_rows:
            accession = locus_row["protein_accession"]
            if accession in protein_accessions:
                raise ValueError(f"Candidate accession duplicates an existing reference/candidate: {accession}")
            header, sequence = fasta_by_id[accession]
            if len(sequence) != int(locus_row["sequence_length_aa"]):
                raise ValueError(f"Candidate sequence length mismatch for {accession}")
            selection = selection_rows.get(accession)
            if selection is None:
                raise KeyError(f"No selection evidence for candidate {accession}")
            core = core_rows.get(accession, {})
            quality = quality_rows.get(accession, {})
            if locus_row["mapping_status"] != "mapped":
                raise ValueError(f"Candidate is not mapped to a locus: {accession}")
            taxon_label = item["taxon"].replace(" ", "_")
            tree_tip = f"{accession}|taxon={taxon_label}|assembly={item['assembly']}|role=Cornus_candidate"
            if tree_tip in tree_tips:
                raise ValueError(f"Duplicate tree tip: {tree_tip}")
            tree_tips.add(tree_tip)
            protein_accessions.add(accession)
            gene_name = locus_row["gene_ids"] or locus_row["locus_tags"] or locus_row["locus_key"]
            final_records.append((tree_tip, sequence))
            final_rows.append({
                "tree_tip": tree_tip,
                "protein_accession": accession,
                "taxon": item["taxon"],
                "assembly_accession": item["assembly"],
                "gene_name": gene_name,
                "family": "unassigned",
                "subfamily": "unassigned",
                "phylogenetic_group": "",
                "panel_role": "Cornus_candidate",
                "sequence_length_aa": str(len(sequence)),
                "membership_evidence": "Candidate recovery only; UGT84 membership not established",
                "function_evidence": "No functional evidence inferred",
                "sequence_source_header": header,
                "source_file": str(item["fasta"].relative_to(PROJECT)),
                "group_classification_source": "",
                "locus_key": locus_row["locus_key"],
                "mapping_status": locus_row["mapping_status"],
                "best_reference_query": selection["best_reference_query"],
                "percent_identity": selection["percent_identity"],
                "query_coverage_percent": selection["query_coverage_percent"],
                "evalue": selection["evalue"],
                "selection_reason": selection["selection_reason"],
                "HCGWNS_core": core.get("exact_HCGWNS_core_present", "not_checked"),
                "model_quality_flags": quality.get("model_flags", "not_available"),
            })

    if len(final_records) != 65:
        raise ValueError(f"Expected 65 combined sequences (32 references + 20 Florida + 13 Suecica), got {len(final_records)}")
    if len(tree_tips) != len(final_records) or len(protein_accessions) != len(final_records):
        raise ValueError("Tree-tip or accession duplication detected")
    outdir.mkdir(parents=True, exist_ok=False)
    with fasta_out.open("x", encoding="ascii") as handle:
        for tree_tip, sequence in final_records:
            handle.write(f">{tree_tip}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    with manifest_out.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(final_rows)
    print(f"Combined {len(final_records)} sequences: 32 reference/context, 20 *C. florida*, 13 *C. suecica*.")
    print(f"Final combined input: {fasta_out}")
    print(f"Matching manifest: {manifest_out}")


if __name__ == "__main__":
    main()
