#!/usr/bin/env python3
"""Stage local *C. suecica* assembly and annotation inputs without overwriting."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = PROJECT.parent / "GCA_963922005.1"
DESTINATION = PROJECT / "Inputs" / "cornus_suecica"
ASSEMBLY = "GCA_963922005.1"
INPUTS = (
    ("genome_fasta", "GCA_963922005.1_drCorSuec2.1_genomic.fna"),
    ("gff3", "cornus_suecica_gca963922005v1.gff3"),
    ("gtf", "cornus_suecica_gca963922005v1.gtf"),
    ("sequence_report", "sequence_report.jsonl"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE,
                        help="Directory containing the local GCA_963922005.1 files")
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()
    if DESTINATION.exists():
        raise FileExistsError(f"Refusing to overwrite staged inputs: {DESTINATION}")
    sources = [(kind, name, source_dir / name) for kind, name in INPUTS]
    missing = [(kind, path) for kind, _, path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Suecica source inputs: {missing}")

    DESTINATION.mkdir(parents=True, exist_ok=False)
    entries = []
    for kind, name, source in sources:
        target = DESTINATION / name
        shutil.copy2(source, target)
        source_hash = sha256(source)
        target_hash = sha256(target)
        if source_hash != target_hash:
            raise RuntimeError(f"Copy hash mismatch for {source}")
        try:
            source_record = str(source.relative_to(PROJECT.parent))
        except ValueError:
            source_record = source.name
        entries.append({
            "input_type": kind,
            "source_path_at_staging": source_record,
            "project_path": str(target.relative_to(PROJECT)),
            "size_bytes": target.stat().st_size,
            "sha256": target_hash,
        })

    manifest = {
        "taxon": "Cornus suecica",
        "assembly_accession": ASSEMBLY,
        "annotation_header_build": "drCorSuec2.1",
        "staged_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": entries,
        "note": "Copies of local source files; original assembly and annotations were not modified.",
    }
    with (DESTINATION / "source_manifest.json").open("x", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    print(f"Copied and hash-verified {len(entries)} Suecica inputs under {DESTINATION}")


if __name__ == "__main__":
    main()
