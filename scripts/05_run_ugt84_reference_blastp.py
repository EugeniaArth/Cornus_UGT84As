#!/usr/bin/env python3
"""Run reference-guided BLASTP against broad HMMER-hit proteins; retain raw hits."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

OUTFMT = "6 qseqid sseqid pident length qcovs qstart qend sstart send evalue bitscore"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query", type=Path, default=Path(__file__).resolve().parents[1] / "references/ugt84_expanded/sequences/formal_ugt84_reference_seeds.faa")
    parser.add_argument("--subject", required=True, type=Path, help="Broad PF00201-hit protein FASTA")
    parser.add_argument("--outdir", required=True, type=Path, help="New, non-existing output directory")
    parser.add_argument("--evalue", default="1e-5")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    query, subject = args.query.resolve(), args.subject.resolve()
    if not query.is_file() or not subject.is_file():
        raise FileNotFoundError("--query and --subject must be existing FASTA files")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    outputs = {"hits": outdir / "ugt84_reference_blastp.tsv", "log": outdir / "blastp.log", "manifest": outdir / "run_manifest.json"}
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError(f"Refusing to overwrite outputs in {outdir}")
    executable = shutil.which("blastp")
    if not executable:
        raise FileNotFoundError("blastp is not available on PATH")
    command = [executable, "-query", str(query), "-subject", str(subject), "-evalue", args.evalue,
               "-max_target_seqs", "1000", "-max_hsps", "1", "-seg", "yes", "-num_threads", str(args.threads), "-outfmt", OUTFMT,
               "-out", str(outputs["hits"])]
    with outputs["log"].open("x") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    manifest = {
        "tool": "blastp",
        "command": command,
        "exit_status": result.returncode,
        "recorded_utc": datetime.now(timezone.utc).isoformat(),
        "query_fasta": str(query),
        "query_sha256": sha256(query),
        "subject_fasta": str(subject),
        "subject_sha256": sha256(subject),
        "evalue_cutoff": args.evalue,
        "outfmt": OUTFMT,
        "outputs": {key: str(path) for key, path in outputs.items()},
        "interpretation": "Similarity evidence for candidate recovery only; not a UGT84 or function assignment",
    }
    with outputs["manifest"].open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    print(f"Raw BLASTP hits and run manifest written under {outdir}")


if __name__ == "__main__":
    main()
