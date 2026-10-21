#!/usr/bin/env python3
"""Render detailed and PowerPoint-friendly variants of the 62-tip tree."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from Bio import Phylo

SPECIES_COLORS = {"Cornus florida": "#0072B2", "Cornus suecica": "#CC79A7"}
REFERENCE_COLOR = "#333333"
BRANCH_COLOR = "#68717A"
SUBFAMILY_COLORS = {"A": "#6A3D9A", "B": "#E69F00", "J": "#009E73"}
TAXON_SHORT = {
    "Quercus robur": "Q. robur", "Camellia sinensis": "C. sinensis",
    "Punica granatum": "P. granatum", "Eucalyptus camaldulensis": "E. camaldulensis",
    "Vitis vinifera": "V. vinifera", "Fragaria x ananassa": "F. × ananassa",
    "Canarium album": "C. album", "Litsea cubeba": "L. cubeba",
    "Citrus unshiu": "C. unshiu", "Brassica napus": "B. napus",
    "Populus fremontii x P. angustifolia": "P. fremontii × P. angustifolia",
    "Arabidopsis thaliana": "A. thaliana", "Stevia rebaudiana": "S. rebaudiana",
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


def support_value(clade) -> float | None:
    """Return the UFBoot part of IQ-TREE's SH-aLRT/UFBoot node label."""
    label = str(clade.name or "")
    if "/" not in label:
        return None
    try:
        return float(label.rsplit("/", 1)[1])
    except ValueError:
        return None


def short_gene(row: dict[str, str]) -> str:
    gene = row["gene_name"]
    if row["panel_role"] == "UGT84_reference" and " / " in gene:
        return gene.rsplit(" / ", 1)[1]
    return gene


def label_detailed(row: dict[str, str]) -> str:
    role, taxon = row["panel_role"], TAXON_SHORT.get(row["taxon"], row["taxon"])
    if row["taxon"] == "Cornus florida":
        return f"C. florida · {row['gene_name'].removeprefix('gene-')}"
    if row["taxon"] == "Cornus suecica":
        return f"C. suecica · {row['protein_accession']}"
    gene = short_gene(row)
    if role == "Group_L_UGT74_comparator":
        return f"{gene} · {taxon} — Group L UGT74 comparator"
    if role == "family1_context_outgroup":
        return f"{gene} · {taxon} (Group {row['phylogenetic_group']})"
    return f"{gene} · {taxon}"


def label_ppt(row: dict[str, str]) -> str:
    role = row["panel_role"]
    if row["taxon"] == "Cornus florida":
        return f"Cfl · {row['gene_name'].removeprefix('gene-')}"
    if row["taxon"] == "Cornus suecica":
        accession = row["protein_accession"]
        compact = accession[:7] + accession[-5:] if accession.startswith("ENSLUWP") else accession
        return f"Csue · {compact}"
    gene = short_gene(row)
    if role == "Group_L_UGT74_comparator":
        return f"{gene} · Group L UGT74 comparator"
    if role == "family1_context_outgroup":
        return f"{gene} (Group {row['phylogenetic_group']})"
    return gene


def build_positions(tree):
    tree.ladderize(reverse=True)
    tips = tree.get_terminals()
    y_by_node = {tip: float(len(tips) - index - 1) for index, tip in enumerate(tips)}
    x_by_node = {tree.root: 0.0}

    def assign(node):
        if node.is_terminal():
            return y_by_node[node]
        child_y = []
        for child in node.clades:
            x_by_node[child] = x_by_node[node] + float(child.branch_length or 0.0)
            child_y.append(assign(child))
        y_by_node[node] = sum(child_y) / len(child_y)
        return y_by_node[node]

    assign(tree.root)
    max_x = max(x_by_node.values())
    if max_x <= 0:
        raise ValueError("Tree has no positive branch lengths")
    return tips, x_by_node, y_by_node, max_x


def render_variant(*, tree, tips, x_by_node, y_by_node, max_x, rows_by_label,
                   label_func, output_base: Path, variant: str,
                   support_threshold: float, figsize: tuple[float, float],
                   label_font: float, support_font: float) -> dict:
    """Render a radial (unrooted-style) phylogram; orientation is arbitrary."""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    is_ppt = variant == "ppt"
    title = ("UGT84-reference placement of Cornus candidates" if is_ppt
             else "Cornus candidates with UGT84 reference panel")
    subtitle = ("Radial unrooted-style ML tree · 62 proteins · orientation arbitrary"
                if is_ppt else
                "Radial unrooted-style maximum-likelihood tree · 62 proteins · orientation arbitrary")
    png_dpi = 300
    marker_size = 6.5 if is_ppt else 9

    # Assign evenly spaced angles to tips, then place internal nodes at the
    # circular mean of their descendants. Branch radii use cumulative lengths.
    ordered_tips = list(tips)
    angle_by_node = {tip: 2 * np.pi * i / len(ordered_tips)
                     for i, tip in enumerate(ordered_tips)}
    radius_by_node = {tree.root: 0.0}

    def place(node):
        if node.is_terminal():
            return angle_by_node[node]
        # Assign each child's radius BEFORE recursing into it. Internal child
        # nodes need their own radius when their descendants are positioned.
        child_angles = []
        for child in node.clades:
            radius_by_node[child] = radius_by_node[node] + float(child.branch_length or 0.0)
            child_angles.append(place(child))
        z = sum(np.exp(1j * a) for a in child_angles)
        angle_by_node[node] = float(np.angle(z) % (2 * np.pi)) if abs(z) > 1e-12 else child_angles[0]
        return angle_by_node[node]

    place(tree.root)
    max_radius = max(radius_by_node.values())
    if max_radius <= 0:
        raise ValueError("Tree has no positive branch lengths")

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    ax.set_aspect("equal")
    ax.axis("off")

    # Draw each edge as a radial branch plus a short circular connector at
    # the parent's radius. This keeps branch length encoded by radial distance.
    for parent in tree.find_clades(order="level"):
        for child in parent.clades:
            rp = radius_by_node[parent]
            rc = radius_by_node[child]
            ap = angle_by_node[parent]
            ac = angle_by_node[child]
            delta = (ac - ap + np.pi) % (2 * np.pi) - np.pi
            if rp > 0 and abs(delta) > 1e-8:
                arc_angles = np.linspace(ap, ap + delta, max(8, int(abs(delta) * 30)))
                ax.plot(rp * np.cos(arc_angles), rp * np.sin(arc_angles),
                        color=BRANCH_COLOR, linewidth=0.55 if not is_ppt else 0.45,
                        solid_capstyle="round", zorder=1)
            ax.plot([rp * np.cos(ac), rc * np.cos(ac)],
                    [rp * np.sin(ac), rc * np.sin(ac)],
                    color=BRANCH_COLOR, linewidth=0.55 if not is_ppt else 0.45,
                    solid_capstyle="round", zorder=1)

    support_count = 0
    for clade in tree.find_clades(order="level"):
        if clade is tree.root or clade.is_terminal():
            continue
        value = support_value(clade)
        if value is not None and value >= support_threshold:
            r, a = radius_by_node[clade], angle_by_node[clade]
            ax.text((r + 0.012 * max_radius) * np.cos(a),
                    (r + 0.012 * max_radius) * np.sin(a),
                    str(int(round(value))), fontsize=support_font,
                    color="#4B5563", ha="center", va="center", zorder=4)
            support_count += 1

    # Anchor every label to its own terminal node, then nudge it outward in
    # screen points. This keeps labels attached to the actual branch tips even
    # when terminal branch lengths differ.
    max_label_offset_pt = 8.0 if is_ppt else 9.0
    for tip in tips:
        row = rows_by_label[tip.name]
        a = angle_by_node[tip]
        r = radius_by_node[tip]
        x, y = r * np.cos(a), r * np.sin(a)
        species_color = SPECIES_COLORS.get(row["taxon"], REFERENCE_COLOR)
        ax.scatter([x], [y], s=marker_size, color=species_color, linewidths=0, zorder=3)

        right_side = np.cos(a) >= 0
        # Offset the text radially outward from this tip. Matplotlib interprets
        # xytext in points, so the label remains a small, consistent distance
        # from the tip regardless of the branch length or plot scale.
        dx = max_label_offset_pt * np.cos(a)
        dy = max_label_offset_pt * np.sin(a)
        ax.annotate(
            label_func(row),
            xy=(x, y),
            xycoords="data",
            xytext=(dx, dy),
            textcoords="offset points",
            color=species_color,
            fontsize=label_font,
            va="center",
            ha="left" if right_side else "right",
            annotation_clip=False,
            clip_on=False,
            zorder=5,
        )
        if row["panel_role"] == "UGT84_reference" and row["subfamily"] in SUBFAMILY_COLORS:
            # Subfamily identity is retained in the manifest and label maps.
            # Candidate Cornus subfamilies are intentionally not inferred here.
            pass

    # Leave room for labels outside the longest tip. This is a data-coordinate
    # margin; bbox_inches="tight" will include text that extends beyond the axes.
    extent = max_radius * (1.35 if is_ppt else 1.45)
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    fig.suptitle(title, y=0.985, fontsize=13 if is_ppt else 14,
                 fontweight="bold", color="#202124")
    fig.text(0.5, 0.955, subtitle, ha="center", va="center",
             fontsize=8 if is_ppt else 9, color="#4B5563")

    species_handles = [
        Line2D([0], [0], color=SPECIES_COLORS["Cornus florida"], marker="o",
               linestyle="none", markersize=4.5, label="C. florida candidates"),
        Line2D([0], [0], color=SPECIES_COLORS["Cornus suecica"], marker="o",
               linestyle="none", markersize=4.5, label="C. suecica candidates"),
        Line2D([0], [0], color=REFERENCE_COLOR, marker="o",
               linestyle="none", markersize=4.5, label="References / context"),
    ]
    fig.legend(handles=species_handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.035), ncol=3, frameon=False,
               fontsize=6.2 if is_ppt else 7)
    footer = ("Three start-codon-flagged Suecica candidates are held out; "
              "UGT74 labels identify Group L UGT74 comparators; Cornus subfamilies remain unassigned. "
              f"UFBoot values shown only for ≥{int(support_threshold)}%.")
    fig.text(0.5, 0.012, footer, ha="center", va="center",
             fontsize=5.8 if is_ppt else 7, color="#555555")

    outputs = {"png": output_base.with_suffix(".png"),
               "svg": output_base.with_suffix(".svg"),
               "pdf": output_base.with_suffix(".pdf")}
    for path in outputs.values():
        if path.exists():
            raise FileExistsError(path)
    fig.savefig(outputs["png"], dpi=png_dpi, facecolor="white", bbox_inches="tight")
    fig.savefig(outputs["svg"], facecolor="white", bbox_inches="tight")
    fig.savefig(outputs["pdf"], facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return {"files": outputs, "support_threshold": support_threshold,
            "support_values_rendered": support_count, "figure_inches": list(figsize),
            "tip_label_font_pt": label_font, "support_font_pt": support_font}

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path, help="New output directory")
    args = parser.parse_args()
    tree_path, input_manifest, outdir = args.tree.resolve(), args.manifest.resolve(), args.outdir.resolve()
    if not tree_path.is_file() or not input_manifest.is_file():
        raise FileNotFoundError("Tree and matching manifest must exist")
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite figure output directory: {outdir}")

    with input_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 62:
        raise ValueError(f"Expected 62 manifest rows, found {len(rows)}")
    rows_by_label = {iqtree_label(row["tree_tip"]): row for row in rows}
    if len(rows_by_label) != len(rows):
        raise ValueError("IQ-TREE label normalization creates duplicate tips")
    tree = Phylo.read(str(tree_path), "newick")
    if {tip.name for tip in tree.get_terminals()} != set(rows_by_label):
        raise ValueError("Tree tips do not match the 62-tip manifest")

    detailed_labels = {tip.name: label_detailed(rows_by_label[tip.name]) for tip in tree.get_terminals()}
    ppt_labels = {tip.name: label_ppt(rows_by_label[tip.name]) for tip in tree.get_terminals()}
    if len(set(detailed_labels.values())) != 62 or len(set(ppt_labels.values())) != 62:
        raise ValueError("One or more display-label sets are not unique")
    for path in (outdir,):
        if path.exists():
            raise FileExistsError(path)
    outdir.mkdir(parents=True, exist_ok=False)
    cache = outdir / ".mplconfig"
    cache.mkdir()
    os.environ["MPLCONFIGDIR"] = str(cache)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "pdf.fonttype": 42})

    tree.ladderize(reverse=True)
    tips, x_by_node, y_by_node, max_x = build_positions(tree)
    detailed = render_variant(
        tree=tree, tips=tips, x_by_node=x_by_node, y_by_node=y_by_node, max_x=max_x,
        rows_by_label=rows_by_label, label_func=lambda row: detailed_labels[iqtree_label(row["tree_tip"])],
        output_base=outdir / "cornus_ugt84_62_tree_detailed", variant="detailed",
        support_threshold=95, figsize=(15.5, 14.5), label_font=7.4, support_font=5.8)
    ppt = render_variant(
        tree=tree, tips=tips, x_by_node=x_by_node, y_by_node=y_by_node, max_x=max_x,
        rows_by_label=rows_by_label, label_func=lambda row: ppt_labels[iqtree_label(row["tree_tip"])],
        output_base=outdir / "cornus_ugt84_62_tree_ppt", variant="ppt",
        support_threshold=99, figsize=(13.333, 7.5), label_font=5.2, support_font=4.2)

    detailed_map = outdir / "tree_plot_tip_label_map_detailed.tsv"
    ppt_map = outdir / "tree_plot_tip_label_map_ppt.tsv"
    for map_path, label_map in ((detailed_map, detailed_labels), (ppt_map, ppt_labels)):
        with map_path.open("x", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["iqtree_tip_label", "input_tree_tip", "protein_accession", "taxon",
                             "panel_role", "gene_name", "UGT84_subfamily", "plot_label"])
            for tip in tips:
                row = rows_by_label[tip.name]
                writer.writerow([tip.name, row["tree_tip"], row["protein_accession"], row["taxon"],
                                 row["panel_role"], row["gene_name"], row["subfamily"], label_map[tip.name]])

    variants = {"detailed": detailed, "ppt": ppt}
    output_files = {}
    for variant, info in variants.items():
        output_files[variant] = {
            "figure_inches": info["figure_inches"],
            "support_threshold_percent": info["support_threshold"],
            "support_values_rendered": info["support_values_rendered"],
            "tip_label_font_pt": info["tip_label_font_pt"],
            "support_font_pt": info["support_font_pt"],
            "files": {key: {"path": str(path), "bytes": path.stat().st_size,
                            "sha256": sha256(path)} for key, path in info["files"].items()},
        }
    plot_manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_tree": str(tree_path),
        "source_tree_sha256": sha256(tree_path),
        "source_manifest": str(input_manifest),
        "source_manifest_sha256": sha256(input_manifest),
        "tip_count": 62,
        "tree_is_unrooted": True,
        "display_orientation": "arbitrary Newick trifurcation; not a biological root",
        "bootstrap_display": "UFBoot values parsed from the second component of IQ-TREE SH-aLRT/UFBoot labels; only values above each variant's threshold are shown.",
        "detailed_label_role": "UGT74 comparator tips explicitly labeled Group L UGT74 comparator.",
        "ugt84_subfamily_badges": "A/B/J shown only for formal UGT84 references; no Cornus candidate subfamily assigned.",
        "candidate_species_colors": SPECIES_COLORS,
        "ugt84_subfamily_colors": SUBFAMILY_COLORS,
        "variants": output_files,
        "tip_label_maps": {
            "detailed": {"path": str(detailed_map), "sha256": sha256(detailed_map)},
            "ppt": {"path": str(ppt_map), "sha256": sha256(ppt_map)},
        },
        "matplotlib_version": matplotlib.__version__,
        "biopython_version": __import__("Bio").__version__,
    }
    provenance = outdir / "tree_plot_manifest.json"
    with provenance.open("x") as handle:
        json.dump(plot_manifest, handle, indent=2)
        handle.write("\n")
    print(f"Detailed figure support values (UFBoot >=95): {detailed['support_values_rendered']}")
    print(f"PPT figure support values (UFBoot >=99): {ppt['support_values_rendered']}")
    print(f"Saved detailed and 16:9 presentation variants under {outdir}")


if __name__ == "__main__":
    main()
