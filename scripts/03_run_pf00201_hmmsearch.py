#!/usr/bin/env python3
"""Run a broad PF00201 HMMER screen and save raw outputs plus provenance."""
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
    parser.add_argument("--hmm", required=True, type=Path, help="Local PF00201 HMM profile; not distributed here")
    parser.add_argument("--proteins", required=True, type=Path, help="Annotated protein FASTA for one taxon")
    parser.add_argument("--outdir", required=True, type=Path, help="New, non-existing output directory")
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()
    hmm, proteins = args.hmm.resolve(), args.proteins.resolve()
    if not hmm.is_file() or not proteins.is_file():
        raise FileNotFoundError("Both --hmm and --proteins must be files")
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    names = ["pf00201.tblout", "pf00201.domtblout", "hmmsearch.log", "run_manifest.json"]
    if any((outdir / name).exists() for name in names):
        raise FileExistsError(f"Refusing to overwrite output in {outdir}")
    executable = shutil.which("hmmsearch")
    if not executable:
        raise FileNotFoundError("hmmsearch is not available on PATH")
    command = [executable, "--cut_ga", "--cpu", str(args.threads), "--noali",
               "--tblout", str(outdir / names[0]), "--domtblout", str(outdir / names[1]),
               str(hmm), str(proteins)]
    started = datetime.now(timezone.utc).isoformat()
    with (outdir / names[2]).open("x") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    manifest = {
        "tool": "hmmsearch",
        "command": command,
        "exit_status": result.returncode,
        "started_utc": started,
        "profile": str(hmm),
        "profile_sha256": sha256(hmm),
        "proteins": str(proteins),
        "proteins_sha256": sha256(proteins),
        "outputs": {name: str(outdir / name) for name in names[:3]},
        "profile_source": "Pfam PF00201 / InterPro; broad UDP-glycosyltransferase-family screen, not a UGT84 classifier",
    }
    with (outdir / names[3]).open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    print(f"HMMER raw outputs and manifest written under {outdir}")


if __name__ == "__main__":
    main()
