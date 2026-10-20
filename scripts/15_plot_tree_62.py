#!/usr/bin/env python3
"""Plot the 62-tip Cornus UGT84-reference tree with minimal annotations."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

from Bio import Phylo

SPECIES_COLORS = {
    "Cornus florida": "#0072B2",  # colorblind-safe blue
    "Cornus suecica": "#CC79A7",  # colorblind-safe magenta
}
REFERENCE_COLOR = "#333333"
BRANCH_COLOR = "#68717A"
SUBFAMILY_COLORS = {"A": "#6A3D9A", "B": "#E69F00", "J": "#009E73"}

TAXON_SHORT = {
    "Quercus robur": "Q. robur",
    "Camellia sinensis": "C. sinensis",
    "Punica granatum": "P. granatum",
    "Eucalyptus camaldulensis": "E. camaldulensis",
    "Vitis vinifera": "V. vinifera",
    "Fragaria x ananassa": "F. × ananassa",
    "Canarium album": "C. album",
    "Litsea cubeba": "L. cubeba",
    "Citrus unshiu": "C. unshiu",
    "Brassica napus": "B. napus",
    "Populus fremontii x P. angustifolia": "P. fremontii × P. angustifolia",
    "Arabidopsis thaliana": "A. thaliana",
    "Stevia rebaudiana": "S. rebaudiana",
    "Medicago truncatula": "M. truncatula",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iqtree_label(identifier: str) -> str:
    return identifier.replace("=", "_")


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def display_label(row: dict[str, str]) -> str:
    taxon = TAXON_SHORT.get(row["taxon"], row["taxon"])
    role = row["panel_role"]
    if row["taxon"] == "Cornus florida":
        gene = row["gene_name"].removeprefix("gene-")
        return f"C. florida · {gene}"
    if row["taxon"] == "Cornus suecica":
        return f"C. suecica · {row['protein_accession']}"
    gene = row["gene_name"]
    if role == "family1_context_outgroup":
        group = row["phylogenetic_group"]
        return f"{gene} · {taxon} (Group {group})"
    return f"{gene} · {taxon}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path, help="62-tip unrooted IQ-TREE .treefile")
    parser.add_argument("--manifest", required=True, type=Path, help="Manifest corresponding to tree tips")
    parser.add_argument("--outdir", required=True, type=Path, help="New figure output directory")
    args = parser.parse_args()

    tree_path = args.tree.resolve()
    manifest_path = args.manifest.resolve()
    outdir = args.outdir.resolve()
    if not tree_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError("Tree and matching manifest must exist")
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite figure directory: {outdir}")

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 62:
        raise ValueError(f"Expected 62 manifest rows; found {len(rows)}")
    by_label = {iqtree_label(row["tree_tip"]): row for row in rows}
    if len(by_label) != len(rows):
        raise ValueError("Tip-label normalization creates duplicate labels")
    tree = Phylo.read(str(tree_path), "newick")
    tips = tree.get_terminals()
    if len(tips) != 62 or {tip.name for tip in tips} != set(by_label):
        raise ValueError("Tree tips do not match the 62-tip manifest")

    # Build concise display labels and prevent accidental collisions.
    label_by_tip = {tip: display_label(by_label[tip.name]) for tip in tips}
    if len(set(label_by_tip.values())) != len(label_by_tip):
        raise ValueError("Display labels are not unique; add accession suffixes before plotting")
    outdir.mkdir(parents=True, exist_ok=False)
    cache = outdir / ".mplconfig"
    cache.mkdir()
    os.environ["MPLCONFIGDIR"] = str(cache)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
    })

    # Ladderizing changes only display order. The input tree remains unrooted;
    # the Newick trifurcation is a display anchor, not an inferred root.
    tree.ladderize(reverse=True)
    tips = tree.get_terminals()
    y_by_node: dict[object, float] = {tip: float(len(tips) - i - 1) for i, tip in enumerate(tips)}
    x_by_node: dict[object, float] = {tree.root: 0.0}

    def assign_positions(node) -> float:
        if node.is_terminal():
            return y_by_node[node]
        child_y = []
        for child in node.clades:
            x_by_node[child] = x_by_node[node] + float(child.branch_length or 0.0)
            child_y.append(assign_positions(child))
        y_by_node[node] = sum(child_y) / len(child_y)
        return y_by_node[node]

    assign_positions(tree.root)
    max_x = max(x_by_node.values())
    if max_x <= 0:
        raise ValueError("Tree has no positive branch lengths")

    fig = plt.figure(figsize=(15.5, 14.5), facecolor="white")
    ax_tree = fig.add_axes([0.045, 0.13, 0.55, 0.78])
    ax_labels = fig.add_axes([0.61, 0.13, 0.285, 0.78], sharey=ax_tree)
    ax_subfamily = fig.add_axes([0.91, 0.13, 0.045, 0.78], sharey=ax_tree)

    for parent in tree.find_clades(order="level"):
        children = parent.clades
        if not children:
            continue
        parent_x = x_by_node[parent]
        child_ys = [y_by_node[child] for child in children]
        if len(children) > 1:
            ax_tree.plot([parent_x, parent_x], [min(child_ys), max(child_ys)],
                         color=BRANCH_COLOR, linewidth=0.55, solid_capstyle="round", zorder=1)
        for child in children:
            child_x = x_by_node[child]
            child_y = y_by_node[child]
            ax_tree.plot([parent_x, child_x], [child_y, child_y], color=BRANCH_COLOR,
                         linewidth=0.55, solid_capstyle="round", zorder=1)

    for tip in tips:
        row = by_label[tip.name]
        y = y_by_node[tip]
        species_color = SPECIES_COLORS.get(row["taxon"], REFERENCE_COLOR)
        ax_tree.scatter([x_by_node[tip]], [y], s=9, color=species_color, linewidths=0, zorder=3)
        ax_labels.text(0.0, y, label_by_tip[tip], color=species_color,
                       fontsize=7.7, va="center", ha="left", clip_on=False)
        if row["panel_role"] == "UGT84_reference" and row["subfamily"] in SUBFAMILY_COLORS:
            subfamily = row["subfamily"]
            ax_subfamily.text(0.13, y, subfamily, color=SUBFAMILY_COLORS[subfamily],
                              fontsize=8.0, fontweight="bold", va="center", ha="left")

    ax_tree.set_xlim(-0.015 * max_x, 1.025 * max_x)
    ax_tree.set_ylim(-2.5, len(tips) - 0.5)
    ax_tree.invert_yaxis()
    ax_tree.axis("off")
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(ax_tree.get_ylim())
    ax_labels.axis("off")
    ax_subfamily.set_xlim(0, 1)
    ax_subfamily.set_ylim(ax_tree.get_ylim())
    ax_subfamily.axis("off")
    ax_subfamily.text(0.1, 1.012, "UGT84\nsubfamily", transform=ax_subfamily.transAxes,
                      fontsize=7.2, color="#333333", ha="left", va="bottom")

    # Small branch-length scale bar; the Newick root orientation is arbitrary.
    scale = 0.1
    if scale < max_x * 0.35:
        scale_y = 0.035  # axes fraction: keep the scale bar at the bottom regardless of y orientation
        ax_tree.plot([0.02 * max_x, 0.02 * max_x + scale], [scale_y, scale_y],
                     transform=ax_tree.get_xaxis_transform(), color="#333333",
                     linewidth=1.0, clip_on=False)
        ax_tree.text(0.02 * max_x + scale / 2, scale_y - 0.012, "0.1 substitutions/site",
                     transform=ax_tree.get_xaxis_transform(), fontsize=6.8,
                     color="#333333", ha="center", va="top", clip_on=False)

    fig.suptitle("Cornus candidates with UGT84 reference panel", y=0.982,
                 fontsize=14, fontweight="bold", color="#202124")
    fig.text(0.5, 0.955,
             "Unrooted maximum-likelihood topology · 62 proteins · display orientation arbitrary",
             ha="center", va="center", fontsize=9, color="#4B5563")

    species_handles = [
        Line2D([0], [0], color=SPECIES_COLORS["Cornus florida"], marker="o", linestyle="none",
               markersize=5, label="C. florida candidates"),
        Line2D([0], [0], color=SPECIES_COLORS["Cornus suecica"], marker="o", linestyle="none",
               markersize=5, label="C. suecica candidates"),
        Line2D([0], [0], color=REFERENCE_COLOR, marker="o", linestyle="none",
               markersize=5, label="References / context"),
    ]
    fig.legend(handles=species_handles, loc="lower center", bbox_to_anchor=(0.33, 0.055),
               ncol=3, frameon=False, fontsize=8)
    subfamily_handles = [
        Line2D([0], [0], color=SUBFAMILY_COLORS[key], marker="s", linestyle="none",
               markersize=5, label=key) for key in ("A", "B", "J")
    ]
    fig.legend(handles=subfamily_handles, title="UGT84 reference subfamily",
               loc="lower center", bbox_to_anchor=(0.72, 0.052), ncol=4,
               frameon=False, fontsize=8, title_fontsize=8)
    fig.text(0.5, 0.018,
             "Three flagged Suecica candidates are held out; A/B/J markers apply only to named UGT84 references.",
             ha="center", va="center", fontsize=7.5, color="#555555")

    outputs = {
        "png": outdir / "cornus_ugt84_62_tree.png",
        "svg": outdir / "cornus_ugt84_62_tree.svg",
        "pdf": outdir / "cornus_ugt84_62_tree.pdf",
    }
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(path)
    fig.savefig(outputs["png"], dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(outputs["svg"], facecolor="white", bbox_inches="tight")
    fig.savefig(outputs["pdf"], facecolor="white", bbox_inches="tight")
    plt.close(fig)

    map_path = outdir / "tree_plot_tip_label_map.tsv"
    with map_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["iqtree_tip_label", "input_tree_tip", "protein_accession", "taxon",
                         "panel_role", "gene_name", "UGT84_subfamily", "plot_label"])
        for tip in tips:
            row = by_label[tip.name]
            writer.writerow([tip.name, row["tree_tip"], row["protein_accession"], row["taxon"],
                             row["panel_role"], row["gene_name"], row["subfamily"], label_by_tip[tip]])

    plot_manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_tree": str(tree_path),
        "source_tree_sha256": sha256(tree_path),
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": sha256(manifest_path),
        "tip_count": len(tips),
        "plot_description": "Unrooted rectangular phylogram; Newick trifurcation/display orientation is arbitrary, not a biological root.",
        "branch_lengths": "IQ-TREE amino-acid substitutions per site; displayed to scale",
        "node_support_labels": "Omitted for minimal presentation; present in the source Newick/IQ-TREE outputs (SH-aLRT/UFBoot).",
        "species_label_colors": SPECIES_COLORS,
        "reference_context_color": REFERENCE_COLOR,
        "UGT84_reference_subfamily_colors": SUBFAMILY_COLORS,
        "candidate_subfamily_assignment": "None; candidates are not assigned UGT84 subfamilies by this plot.",
        "output_files": {key: {"path": str(value), "sha256": sha256(value), "bytes": value.stat().st_size}
                         for key, value in outputs.items()},
        "tip_label_map": str(map_path),
        "tip_label_map_sha256": sha256(map_path),
        "matplotlib_version": matplotlib.__version__,
        "biopython_version": __import__("Bio").__version__,
    }
    manifest_path_out = outdir / "tree_plot_manifest.json"
    with manifest_path_out.open("x") as handle:
        json.dump(plot_manifest, handle, indent=2)
        handle.write("\n")

    print("Saved minimal tree plots:")
    for key, path in outputs.items():
        print(f"  {key}: {path}")
    print(f"Tip label mapping: {map_path}")


if __name__ == "__main__":
    main()
