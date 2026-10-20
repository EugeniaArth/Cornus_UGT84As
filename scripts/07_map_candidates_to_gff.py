#!/usr/bin/env python3
"""Map selected candidate protein IDs to GFF3 gene/transcript features."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


def parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip('"')
        else:
            continue
        attrs[unquote(key)] = unquote(value)
    return attrs


def versionless(value: str) -> str:
    return value.rsplit(".", 1)[0] if value.rsplit(".", 1)[-1].isdigit() else value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selected", required=True, type=Path, help="selected_candidate_hits.tsv from script 06")
    parser.add_argument("--gff", required=True, type=Path, help="GFF3 annotation")
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    if not args.selected.is_file() or not args.gff.is_file():
        raise FileNotFoundError("--selected and --gff must exist")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    output = outdir / "candidate_gene_map.tsv"
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite {output}")

    candidates: list[str] = []
    with args.selected.open(newline="") as handle:
        candidates = [row["candidate_protein_id"] for row in csv.DictReader(handle, delimiter="\t")]
    candidate_keys = set(candidates) | {versionless(x) for x in candidates}
    id_attrs: dict[str, dict[str, str]] = {}
    parent_by_id: dict[str, list[str]] = defaultdict(list)
    protein_to_ids: dict[str, list[str]] = defaultdict(list)
    feature_type: dict[str, str] = {}

    with args.gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                continue
            kind, attrs = cols[2], parse_attributes(cols[8])
            identifier = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("gene_id")
            parents = [p for p in attrs.get("Parent", "").split(",") if p]
            if identifier:
                id_attrs[identifier] = attrs
                feature_type[identifier] = kind
                parent_by_id[identifier].extend(parents)
            if attrs.get("protein_id"):
                for protein_id in attrs["protein_id"].split(","):
                    linked_features = ([identifier] if identifier else []) + parents
                    protein_to_ids[protein_id].extend(linked_features)
                    protein_to_ids[versionless(protein_id)].extend(linked_features)
            # Some GFFs identify proteins in ID or Name rather than protein_id.
            for value in (identifier, attrs.get("Name"), attrs.get("proteinId")):
                if value:
                    protein_to_ids[value].append(identifier or "")
                    protein_to_ids[versionless(value)].append(identifier or "")

    def ancestors(identifier: str) -> set[str]:
        found: set[str] = set()
        frontier = [identifier]
        seen: set[str] = set()
        while frontier and len(seen) < 10000:
            current = frontier.pop()
            if current in seen:
                continue
            seen.add(current)
            attrs = id_attrs.get(current, {})
            if feature_type.get(current, "").lower() == "gene" or current.startswith("gene-"):
                found.add(current)
            for key in ("gene_id", "locus_tag"):
                if attrs.get(key):
                    found.add(attrs[key])
            frontier.extend(parent_by_id.get(current, []))
        return {gene_id.removeprefix("gene:") for gene_id in found}

    rows = []
    for protein_id in candidates:
        feature_ids = set(protein_to_ids.get(protein_id, [])) | set(protein_to_ids.get(versionless(protein_id), []))
        feature_ids.discard("")
        gene_ids: set[str] = set()
        transcripts: set[str] = set()
        locus_tags: set[str] = set()
        gene_names: set[str] = set()
        for feature_id in feature_ids:
            if feature_type.get(feature_id, "").lower() in {"mrna", "transcript", "rna"}:
                transcripts.add(feature_id)
            for gene_id in ancestors(feature_id):
                gene_ids.add(gene_id)
            attrs = id_attrs.get(feature_id, {})
            if attrs.get("locus_tag"):
                locus_tags.add(attrs["locus_tag"])
            if attrs.get("Name"):
                gene_names.add(attrs["Name"])
        for gene_id in gene_ids:
            attrs = id_attrs.get(gene_id) or id_attrs.get(f"gene:{gene_id}", {})
            if attrs.get("locus_tag"):
                locus_tags.add(attrs["locus_tag"])
            if attrs.get("Name"):
                gene_names.add(attrs["Name"])
        rows.append({
            "protein_accession": protein_id,
            "gene_ids": ";".join(sorted(gene_ids)),
            "transcript_ids": ";".join(sorted(transcripts)),
            "locus_tags": ";".join(sorted(locus_tags)),
            "gene_names": ";".join(sorted(gene_names)),
            "mapping_status": "mapped" if gene_ids or locus_tags else ("feature_found_no_gene" if feature_ids else "unmapped"),
            "source_gff": str(args.gff.resolve()),
        })
    fields = ["protein_accession", "gene_ids", "transcript_ids", "locus_tags", "gene_names", "mapping_status", "source_gff"]
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} candidate-to-GFF mappings to {output}")


if __name__ == "__main__":
    main()
