#!/usr/bin/env python3
"""Summarize supported reference associations in the unrooted Cornus tree.

The classification is descriptive only. It uses the smallest side of each
supported unrooted split (UFBoot >= threshold) that contains a Cornus tip and
one or more reference/context tips. It does not assign UGT84 membership,
subfamily, or biochemical function.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TREE = ROOT / "Results/tree_62/cornus_ugt84_62.treefile"
DEFAULT_TIP_MAP = ROOT / "Results/tree_62/tree_plot_tip_label_map_detailed.tsv"
DEFAULT_EXCLUDED = ROOT / "Results/tree_62/excluded_suecica_start_flagged_candidates.tsv"
DEFAULT_SUMMARY = ROOT / "Results/cornus_candidate_tree_placements.tsv"
DEFAULT_COUNTS = ROOT / "Results/cornus_candidate_tree_placement_counts.tsv"
DEFAULT_PLOT = ROOT / "figures/cornus_candidate_placement_counts.png"

REFERENCE_CLASSES = {
    "UGT84_reference": "UGT84 reference-associated",
    "Group_L_UGT74_comparator": "UGT74 Group-L comparator-associated",
    "family1_context_outgroup": "Group E/H context-associated",
}
NO_ASSOCIATION = "No >=95 UFBoot reference association"
MIXED = "Mixed/ambiguous supported reference associations"
HELD_OUT = "Held out: missing annotated start-codon feature"
CLASS_ORDER = [
    "UGT84 reference-associated",
    "UGT74 Group-L comparator-associated",
    "Group E/H context-associated",
    NO_ASSOCIATION,
    MIXED,
    HELD_OUT,
]

@dataclass
class Node:
    name: str = ""
    ufboot: float | None = None
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    leaf_names: frozenset[str] | None = None


def parse_newick(text: str) -> Node:
    """Parse the simple unquoted IQ-TREE Newick labels used in this project."""
    i = 0

    def token(stoppers: str) -> str:
        nonlocal i
        start = i
        while i < len(text) and text[i] not in stoppers:
            i += 1
        return text[start:i].strip()

    def parse_node() -> Node:
        nonlocal i
        if i >= len(text):
            raise ValueError("Unexpected end of Newick string")
        if text[i] == "(":
            i += 1
            node = Node()
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if i >= len(text):
                    raise ValueError("Unclosed Newick group")
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
                raise ValueError(f"Expected ',' or ')' at Newick offset {i}")
            label = token(":,);")
            if "/" in label:
                try:
                    node.ufboot = float(label.rsplit("/", 1)[1])
                except ValueError:
                    node.ufboot = None
        else:
            node = Node(name=token(":,);"))
            if not node.name:
                raise ValueError(f"Empty leaf label at Newick offset {i}")
        if i < len(text) and text[i] == ":":
            i += 1
            token(",);" )  # Branch lengths are not used in this summary.
        return node

    root = parse_node()
    while i < len(text) and text[i].isspace():
        i += 1
    if i < len(text) and text[i] == ";":
        i += 1
    if text[i:].strip():
        raise ValueError(f"Unparsed text after Newick tree: {text[i:i+40]!r}")
    return root


def leaves(node: Node) -> frozenset[str]:
    if node.leaf_names is None:
        if not node.children:
            node.leaf_names = frozenset([node.name])
        else:
            node.leaf_names = frozenset().union(*(leaves(child) for child in node.children))
    return node.leaf_names


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def supported_small_sides(root: Node, all_tips: frozenset[str], threshold: float):
    """Yield smaller sides of internal unrooted splits meeting the support cutoff."""
    sides = []

    def visit(node: Node) -> None:
        if node.parent is not None and node.children and node.ufboot is not None and node.ufboot >= threshold:
            desc = leaves(node)
            other = all_tips - desc
            if len(desc) != len(other):
                side = desc if len(desc) < len(other) else other
                sides.append((side, node.ufboot))
        for child in node.children:
            visit(child)

    visit(root)
    return sides


def classify_tip(tip: str, splits, tip_map: dict[str, dict[str, str]]) -> dict[str, str]:
    associations = []
    for side, support in splits:
        if tip not in side:
            continue
        refs = [name for name in side if tip_map[name]["panel_role"] in REFERENCE_CLASSES]
        groups = {REFERENCE_CLASSES[tip_map[name]["panel_role"]] for name in refs}
        if groups:
            associations.append((len(side), groups, support, refs))

    if not associations:
        return {
            "placement_category": NO_ASSOCIATION,
            "smallest_supported_split_tips": "",
            "support_ufboot": "",
            "reference_categories_in_split": "",
            "reference_accessions_in_split": "",
        }

    min_size = min(item[0] for item in associations)
    nearest = [item for item in associations if item[0] == min_size]
    signatures = {tuple(sorted(item[1])) for item in nearest}
    if len(signatures) == 1 and len(next(iter(signatures))) == 1:
        category = next(iter(signatures))[0]
    else:
        category = MIXED
    groups = sorted({group for item in nearest for group in item[1]})
    accessions = sorted({tip_map[name]["protein_accession"] for item in nearest for name in item[3]})
    supports = sorted({item[2] for item in nearest})
    return {
        "placement_category": category,
        "smallest_supported_split_tips": str(min_size),
        "support_ufboot": ";".join(f"{value:g}" for value in supports),
        "reference_categories_in_split": ";".join(groups),
        "reference_accessions_in_split": ";".join(accessions),
    }


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_plot(count_rows: list[dict[str, str]], path: Path, threshold: float) -> None:
    taxa = ["Cornus florida", "Cornus suecica"]
    categories = CLASS_ORDER
    colors = {
        "UGT84 reference-associated": "#6A3D9A",
        "UGT74 Group-L comparator-associated": "#65717C",
        "Group E/H context-associated": "#E69F00",
        NO_ASSOCIATION: "#C4CCD2",
        MIXED: "#9A6FB0",
        HELD_OUT: "#CC79A7",
    }
    lookup = {(row["taxon"], row["placement_category"]): int(row["count"]) for row in count_rows}
    fig, ax = plt.subplots(figsize=(12.2, 6.4), constrained_layout=False)
    lefts = [0, 0]
    y_positions = [1, 0]
    for category in categories:
        values = [lookup.get((taxon, category), 0) for taxon in taxa]
        bars = ax.barh(y_positions, values, left=lefts, height=0.55,
                       color=colors[category], edgecolor="white", linewidth=1.2,
                       label=category)
        for bar, value in zip(bars, values):
            if value:
                text_color = "white" if category in {
                    "UGT84 reference-associated", "UGT74 Group-L comparator-associated",
                    "Group E/H context-associated", MIXED,
                } else "#26323B"
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_y() + bar.get_height() / 2,
                        str(value), ha="center", va="center", fontsize=11,
                        fontweight="bold", color=text_color)
        lefts = [lefts[i] + values[i] for i in range(len(taxa))]

    selected = [sum(lookup.get((taxon, cat), 0) for cat in categories) for taxon in taxa]
    for y, total, taxon in zip(y_positions, selected, taxa):
        ax.text(total + 0.35, y, f"n={total} selected", va="center", ha="left",
                fontsize=10.5, color="#303B44")
    ax.set_yticks(y_positions, ["C. florida", "C. suecica"], fontsize=12)
    ax.set_xlim(0, max(selected) + 5)
    ax.set_xlabel("Candidate proteins (locus representatives)", fontsize=11)
    ax.set_title("Cornus candidates by supported reference association in the unrooted tree",
                 fontsize=14, fontweight="bold", pad=14)
    ax.grid(axis="x", color="#E6EAED", linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              frameon=False, fontsize=9.5, columnspacing=1.3, handlelength=1.4)
    fig.text(0.5, 0.015,
             f"Association = unique smallest side of a supported unrooted split (UFBoot ≥ {threshold:g}); "
             "not a subfamily or function assignment. Three start-flagged Suecica candidates are held out.",
             ha="center", va="bottom", fontsize=9, color="#4C5963", wrap=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite figure: {path}")
    fig.savefig(path, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE)
    parser.add_argument("--tip-map", type=Path, default=DEFAULT_TIP_MAP)
    parser.add_argument("--excluded", type=Path, default=DEFAULT_EXCLUDED)
    parser.add_argument("--summary-out", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--counts-out", type=Path, default=DEFAULT_COUNTS)
    parser.add_argument("--plot-out", type=Path, default=DEFAULT_PLOT)
    parser.add_argument("--min-ufboot", type=float, default=95.0)
    args = parser.parse_args()
    if not 0 <= args.min_ufboot <= 100:
        parser.error("--min-ufboot must be between 0 and 100")
    paths = [args.tree, args.tip_map, args.excluded]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing input files: {', '.join(missing)}")
    outputs = [args.summary_out, args.counts_out, args.plot_out]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError(f"Refusing to overwrite outputs: {', '.join(existing)}")

    tree = parse_newick(args.tree.read_text(encoding="utf-8").strip())
    all_tips = leaves(tree)
    map_rows = read_tsv(args.tip_map)
    tip_map = {row["iqtree_tip_label"]: row for row in map_rows}
    if len(tip_map) != len(map_rows):
        raise ValueError("Duplicate IQ-TREE tip labels in tip map")
    if set(tip_map) != set(all_tips):
        raise ValueError(f"Tree/tip-map mismatch: tree={len(all_tips)}, map={len(tip_map)}")

    excluded_rows = read_tsv(args.excluded)
    excluded_by_id = {row["protein_accession"]: row for row in excluded_rows}
    if len(excluded_by_id) != len(excluded_rows):
        raise ValueError("Duplicate protein accessions in exclusion table")
    if set(excluded_by_id) & {row["protein_accession"] for row in map_rows if row["panel_role"] == "Cornus_candidate"}:
        raise ValueError("A held-out candidate unexpectedly occurs in the tree tip map")

    splits = supported_small_sides(tree, all_tips, args.min_ufboot)
    result_rows = []
    for tip in sorted(all_tips):
        row = tip_map[tip]
        if row["panel_role"] != "Cornus_candidate":
            continue
        placement = classify_tip(tip, splits, tip_map)
        result_rows.append({
            "protein_accession": row["protein_accession"],
            "taxon": row["taxon"],
            "assembly_accession": next((part.split("=", 1)[1] for part in row["input_tree_tip"].split("|") if part.startswith("assembly=")), ""),
            "gene_name": row["gene_name"],
            "tree_status": "included",
            **placement,
            "support_threshold_ufboot": f">={args.min_ufboot:g}",
            "interpretation": "descriptive topology association only; not a UGT84/subfamily/function assignment",
        })
    if len(result_rows) != 30:
        raise ValueError(f"Expected 30 Cornus tree tips; found {len(result_rows)}")

    for accession, row in sorted(excluded_by_id.items()):
        result_rows.append({
            "protein_accession": accession,
            "taxon": row.get("taxon", ""),
            "assembly_accession": row.get("assembly_accession", ""),
            "gene_name": row.get("gene_name", ""),
            "tree_status": "held_out",
            "placement_category": HELD_OUT,
            "smallest_supported_split_tips": "",
            "support_ufboot": "",
            "reference_categories_in_split": "",
            "reference_accessions_in_split": "",
            "support_threshold_ufboot": f">={args.min_ufboot:g}",
            "interpretation": "not included in the 62-tip tree; retained in candidate data",
        })
    if len(result_rows) != 33:
        raise ValueError(f"Expected 33 selected Cornus representatives including holdouts; found {len(result_rows)}")

    count_map = Counter((row["taxon"], row["placement_category"]) for row in result_rows)
    taxa = ["Cornus florida", "Cornus suecica"]
    count_rows = [
        {"taxon": taxon, "placement_category": category,
         "count": str(count_map.get((taxon, category), 0))}
        for taxon in taxa for category in CLASS_ORDER
    ]
    fields = [
        "protein_accession", "taxon", "assembly_accession", "gene_name", "tree_status",
        "placement_category", "smallest_supported_split_tips", "support_ufboot",
        "reference_categories_in_split", "reference_accessions_in_split",
        "support_threshold_ufboot", "interpretation",
    ]
    write_tsv(args.summary_out, result_rows, fields)
    write_tsv(args.counts_out, count_rows, ["taxon", "placement_category", "count"])
    make_plot(count_rows, args.plot_out, args.min_ufboot)
    print(f"Wrote {len(result_rows)} placement rows: {args.summary_out}")
    print(f"Wrote count table: {args.counts_out}")
    print(f"Wrote figure: {args.plot_out}")


if __name__ == "__main__":
    main()
