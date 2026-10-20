#!/usr/bin/env python3
"""Align the finalized combined protein FASTA with MAFFT L-INS-i."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="final_phylogeny_input.faa from script 09")
    parser.add_argument("--outdir", required=True, type=Path, help="New output directory")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    input_fasta = args.input.resolve()
    if not input_fasta.is_file():
        raise FileNotFoundError(input_fasta)
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    alignment = outdir / "final_alignment.faa"
    log_path = outdir / "mafft.log"
    manifest_path = outdir / "mafft_run_manifest.json"
    if any(path.exists() for path in (alignment, log_path, manifest_path)):
        raise FileExistsError(f"Refusing to overwrite MAFFT outputs in {outdir}")
    executable = shutil.which("mafft")
    if not executable:
        raise FileNotFoundError("mafft is not available on PATH")
    command = [executable, "--localpair", "--maxiterate", "1000", "--thread", str(args.threads), str(input_fasta)]
    started = datetime.now(timezone.utc).isoformat()
    with alignment.open("x") as output, log_path.open("x") as log:
        result = subprocess.run(command, stdout=output, stderr=log, text=True, check=False)
    manifest = {
        "tool": "MAFFT",
        "command": command,
        "exit_status": result.returncode,
        "started_utc": started,
        "input_fasta": str(input_fasta),
        "input_sha256": sha256(input_fasta),
        "alignment": str(alignment),
        "alignment_sha256": sha256(alignment),
        "log": str(log_path),
    }
    with manifest_path.open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    print(f"Alignment and run provenance written under {outdir}")


if __name__ == "__main__":
    main()
