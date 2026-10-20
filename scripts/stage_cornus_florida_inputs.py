#!/usr/bin/env python3
"""Copy the *C. florida* protein, annotation, and PF00201 profile inputs locally."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
DESTINATION = PROJECT / "Inputs" / "cornus_florida"
INPUTS = (
    ("protein_fasta", "GCF_030987335.1_protein.faa"),
    ("gff3", "GCF_030987335.1_genomic.gff"),
    ("sequence_report", "GCF_030987335.1_sequence_report.jsonl"),
    ("hmm_profile", "PF00201.25.hmm"),
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def recorded_path(path: Path) -> str:
    try:
        return str(path.relative_to(PROJECT.parent))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proteins", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--sequence-report", required=True, type=Path)
    parser.add_argument("--hmm", required=True, type=Path, help="PF00201 profile; expected accession PF00201.25")
    args = parser.parse_args()
    sources = {
        "protein_fasta": args.proteins.resolve(),
        "gff3": args.gff.resolve(),
        "sequence_report": args.sequence_report.resolve(),
        "hmm_profile": args.hmm.resolve(),
    }
    for key, path in sources.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {key}: {path}")
    hmm_header = args.hmm.read_text(errors="replace")[:1000]
    if "ACC   PF00201.25" not in hmm_header:
        raise ValueError("The supplied HMM profile does not declare accession PF00201.25")
    if DESTINATION.exists():
        raise FileExistsError(f"Refusing to overwrite staged inputs: {DESTINATION}")

    DESTINATION.mkdir(parents=True, exist_ok=False)
    entries = []
    for key, output_name in INPUTS:
        source = sources[key]
        target = DESTINATION / output_name
        shutil.copy2(source, target)
        source_hash, target_hash = sha256(source), sha256(target)
        if source_hash != target_hash:
            raise RuntimeError(f"Copy hash mismatch: {source} -> {target}")
        entries.append({
            "input_type": key,
            "source_path_at_staging": recorded_path(source),
            "project_path": str(target.relative_to(PROJECT)),
            "size_bytes": target.stat().st_size,
            "sha256": target_hash,
        })
    manifest = {
        "taxon": "Cornus florida",
        "assembly": "GCF_030987335.1",
        "profile_accession": "PF00201.25",
        "profile_role": "Broad UDP-glycosyltransferase-family screen; not a UGT84 classifier",
        "inputs": entries,
    }
    with (DESTINATION / "source_manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    print(f"Copied and hash-verified {len(entries)} *C. florida* inputs under {DESTINATION}")


if __name__ == "__main__":
    main()
