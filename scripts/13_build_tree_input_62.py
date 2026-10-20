#!/usr/bin/env python3
"""Build a tree-only 62-tip input, retaining all candidates in source data."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in sequences:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
                sequences[current] = []
            else:
                if current is None:
                    raise ValueError("Sequence data encountered before the first FASTA header")
                sequences[current].append(line)
    return {key: "".join(parts) for key, parts in sequences.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path,
                        help="New directory; existing directories are not overwritten")
    args = parser.parse_args()

    fasta = args.input_fasta.resolve()
    manifest = args.input_manifest.resolve()
    outdir = args.outdir.resolve()
    if not fasta.is_file() or not manifest.is_file():
        raise FileNotFoundError("Both the combined FASTA and manifest must exist")
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite tree-input output directory: {outdir}")

    sequences = read_fasta(fasta)
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or "tree_tip" not in fields or "taxon" not in fields or "model_quality_flags" not in fields:
            raise ValueError("Input manifest is missing required columns")
        rows = list(reader)
    by_tip = {row["tree_tip"]: row for row in rows}
    if len(by_tip) != len(rows) or set(by_tip) != set(sequences):
        raise ValueError("Combined FASTA and manifest identifiers do not match uniquely")
    for row in rows:
        if len(sequences[row["tree_tip"]]) != int(row["sequence_length_aa"]):
            raise ValueError(f"Length mismatch for {row['tree_tip']}")

    excluded = [row for row in rows
                if row["taxon"] == "Cornus suecica" and row["model_quality_flags"].strip()]
    if len(excluded) != 3 or any(row["model_quality_flags"].strip() != "start_codon_feature_missing"
                                 for row in excluded):
        raise ValueError("Expected exactly the three Suecica start-codon-flagged candidates")
    included = [row for row in rows if row not in excluded]
    if len(included) != 62:
        raise ValueError(f"Expected 62 tree tips after exclusions; found {len(included)}")

    outdir.mkdir(parents=True, exist_ok=False)
    fasta_out = outdir / "cornus_ugt84_tree_input_62.faa"
    manifest_out = outdir / "cornus_ugt84_tree_input_62_manifest.tsv"
    excluded_out = outdir / "excluded_suecica_start_flagged_candidates.tsv"
    build_manifest_out = outdir / "tree_input_build_manifest.json"

    with fasta_out.open("x") as handle:
        for row in included:
            tip = row["tree_tip"]
            sequence = sequences[tip]
            handle.write(f">{tip}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")
    with manifest_out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(included)
    with excluded_out.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=[*fields, "tree_exclusion_reason"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in excluded:
            writer.writerow({**row, "tree_exclusion_reason":
                             "User-directed holdout: missing annotated start-codon feature; retained in candidate data"})

    role_counts: dict[str, int] = {}
    for row in included:
        role_counts[row["panel_role"]] = role_counts.get(row["panel_role"], 0) + 1
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_fasta": str(fasta),
        "source_fasta_sha256": sha256(fasta),
        "source_manifest": str(manifest),
        "source_manifest_sha256": sha256(manifest),
        "tree_input_fasta": str(fasta_out),
        "tree_input_fasta_sha256": sha256(fasta_out),
        "tree_input_manifest": str(manifest_out),
        "tree_input_manifest_sha256": sha256(manifest_out),
        "excluded_candidates_table": str(excluded_out),
        "excluded_candidates_sha256": sha256(excluded_out),
        "input_tips": len(sequences),
        "tree_tips": len(included),
        "excluded_suecica_candidates": [row["tree_tip"] for row in excluded],
        "excluded_reason": "missing annotated start-codon feature; user directed holdout from tree only",
        "role_counts": role_counts,
    }
    with build_manifest_out.open("x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")

    print(f"Wrote {len(included)} tree tips to {fasta_out}")
    print(f"Retained all {len(excluded)} held-out Suecica candidates in the source candidate dataset; exclusions recorded at {excluded_out}")


if __name__ == "__main__":
    main()
