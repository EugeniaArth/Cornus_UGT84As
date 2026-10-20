#!/usr/bin/env python3
"""Infer an exploratory unrooted IQ-TREE phylogeny from a final alignment."""
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
    parser.add_argument("--alignment", required=True, type=Path, help="final_alignment.faa from script 10")
    parser.add_argument("--outdir", required=True, type=Path, help="New, empty or non-existing output directory")
    parser.add_argument("--prefix", default="ugt84_final")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--ufboot", type=int, default=1000)
    parser.add_argument("--alrt", type=int, default=1000)
    args = parser.parse_args()
    alignment = args.alignment.resolve()
    if not alignment.is_file():
        raise FileNotFoundError(alignment)
    if args.threads < 1 or args.seed < 1 or args.ufboot < 1 or args.alrt < 1:
        raise ValueError("threads, seed, and bootstrap replicate counts must be positive")
    outdir = args.outdir.resolve()
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f"IQ-TREE output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    prefix_path = outdir / args.prefix
    if list(outdir.glob(args.prefix + ".*")):
        raise FileExistsError(f"IQ-TREE prefix already has outputs: {prefix_path}")
    manifest_path = outdir / "iqtree_run_manifest.json"
    if manifest_path.exists():
        raise FileExistsError(manifest_path)
    executable = shutil.which("iqtree2") or shutil.which("iqtree")
    if not executable:
        raise FileNotFoundError("iqtree2/iqtree is not available on PATH")
    command = [executable, "-s", str(alignment), "-m", "MFP", "-B", str(args.ufboot),
               "-alrt", str(args.alrt), "-bnni", "-T", str(args.threads), "-seed", str(args.seed),
               "-pre", str(prefix_path)]
    started = datetime.now(timezone.utc).isoformat()
    log_path = outdir / "iqtree_console.log"
    with log_path.open("x") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, check=False)
    treefile = Path(str(prefix_path) + ".treefile")
    manifest = {
        "tool": "IQ-TREE",
        "command": command,
        "exit_status": result.returncode,
        "started_utc": started,
        "alignment": str(alignment),
        "alignment_sha256": sha256(alignment),
        "treefile": str(treefile) if treefile.exists() else None,
        "treefile_sha256": sha256(treefile) if treefile.exists() else None,
        "console_log": str(log_path),
        "rooting": "unrooted inference; no root specified",
        "support_label_order": "SH-aLRT/UFBoot",
    }
    with manifest_path.open("x") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    print(f"IQ-TREE output and run provenance written under {outdir}")


if __name__ == "__main__":
    main()
