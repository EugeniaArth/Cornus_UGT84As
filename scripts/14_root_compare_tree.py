#!/usr/bin/env python3
"""Create a strict two-outgroup rooted tree and compare it with the unrooted tree."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

from Bio import Phylo


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def iqtree_label(identifier: str) -> str:
    # IQ-TREE converts '=' to '_' in Newick tip labels.
    return identifier.replace("=", "_")


def terminals(clade) -> set[str]:
    names = {tip.name for tip in clade.get_terminals()}
    if None in names:
        raise ValueError("Tree contains an unnamed terminal")
    return names


def canonical_split(side: set[str], all_tips: set[str]) -> tuple[str, ...]:
    other = all_tips - side
    left, right = tuple(sorted(side)), tuple(sorted(other))
    return min((left, right), key=lambda group: (len(group), group))


def nontrivial_splits(tree) -> set[tuple[str, ...]]:
    all_tips = terminals(tree.root)
    splits = set()
    for clade in tree.find_clades(order="level"):
        if clade is tree.root:
            continue
        side = terminals(clade)
        if 1 < len(side) < len(all_tips) - 1:
            splits.add(canonical_split(side, all_tips))
    return splits


def split_labels(tree) -> dict[tuple[str, ...], str]:
    all_tips = terminals(tree.root)
    labels: dict[tuple[str, ...], str] = {}
    for clade in tree.find_clades(order="level"):
        if clade is tree.root or not clade.name:
            continue
        side = terminals(clade)
        if 1 < len(side) < len(all_tips) - 1:
            labels[canonical_split(side, all_tips)] = str(clade.name)
    return labels


def pairwise_distances(tree) -> dict[tuple[str, str], float]:
    names = sorted(terminals(tree.root))
    return {(left, right): float(tree.distance(left, right))
            for i, left in enumerate(names) for right in names[i + 1:]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", required=True, type=Path, help="Unrooted IQ-TREE .treefile")
    parser.add_argument("--manifest", required=True, type=Path, help="Manifest for the exact tree input")
    parser.add_argument("--outdir", required=True, type=Path, help="New output directory")
    args = parser.parse_args()

    tree_path = args.tree.resolve()
    manifest_path = args.manifest.resolve()
    outdir = args.outdir.resolve()
    if not tree_path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError("Tree and input manifest must exist")
    if outdir.exists():
        raise FileExistsError(f"Refusing to overwrite tree-variant directory: {outdir}")

    with manifest_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    outgroup_rows = [row for row in rows if row["panel_role"] == "family1_context_outgroup"]
    if len(rows) != 62 or len(outgroup_rows) != 2:
        raise ValueError(f"Expected 62 tips and two context outgroups; got {len(rows)} and {len(outgroup_rows)}")

    outgroup_accessions = {row["protein_accession"] for row in outgroup_rows}
    outgroup_original_tips = {row["tree_tip"] for row in outgroup_rows}
    outgroup_tips = {iqtree_label(tip) for tip in outgroup_original_tips}
    if outgroup_accessions != {"AAR06912.1", "AAW56092.1"}:
        raise ValueError(f"Unexpected context outgroups: {sorted(outgroup_accessions)}")

    tree = Phylo.read(str(tree_path), "newick")
    all_tips = terminals(tree.root)
    expected_tips = {iqtree_label(row["tree_tip"]) for row in rows}
    if len(expected_tips) != len(rows):
        raise ValueError("IQ-TREE-safe tip-label conversion causes collisions")
    if all_tips != expected_tips:
        raise ValueError("IQ-TREE tips do not exactly match the normalized 62-tip input manifest")

    # Find the unrooted edge that separates exactly the two context sequences
    # from all other taxa. Refuse an unsupported or ambiguous rooting.
    edge_candidates = []
    for clade in tree.find_clades(order="level"):
        if clade is tree.root:
            continue
        side = terminals(clade)
        if side == outgroup_tips:
            edge_candidates.append((0, clade, str(clade.name or "")))
        elif all_tips - side == outgroup_tips:
            edge_candidates.append((1, clade, str(clade.name or "")))
    if not edge_candidates:
        raise ValueError("The two designated context sequences do not define an exclusive outgroup split; no rooted tree written")
    edge_candidates.sort(key=lambda item: item[0])
    target = edge_candidates[0][1]
    root_edge_label = edge_candidates[0][2]

    original_splits = nontrivial_splits(tree)
    original_distances = pairwise_distances(tree)
    labels_before = split_labels(tree)
    rooted = tree
    rooted.root_with_outgroup(target)
    rooted.rooted = True
    rooted_root_children = [terminals(child) for child in rooted.root.clades]
    if len(rooted_root_children) != 2 or outgroup_tips not in rooted_root_children:
        raise ValueError("Rerooting did not place the two designated outgroups on one root side")
    if nontrivial_splits(rooted) != original_splits:
        raise ValueError("Rooted and unrooted trees differ in unrooted split topology")
    rooted_distances = pairwise_distances(rooted)
    max_distance_difference = max(abs(original_distances[key] - rooted_distances[key])
                                  for key in original_distances)
    if max_distance_difference > 1e-8:
        raise ValueError(f"Rerooting changed patristic distances (max difference {max_distance_difference})")

    labels_after = split_labels(rooted)
    label_differences = [
        {"split": list(split), "unrooted_label": labels_before.get(split),
         "rooted_label": labels_after.get(split)}
        for split in sorted(set(labels_before) | set(labels_after))
        if labels_before.get(split) != labels_after.get(split)
    ]

    outdir.mkdir(parents=True, exist_ok=False)
    unrooted_out = outdir / "cornus_ugt84_62_unrooted.treefile"
    rooted_out = outdir / "cornus_ugt84_62_rooted_two_outgroups.treefile"
    mapping_out = outdir / "tree_tip_label_map.tsv"
    comparison_out = outdir / "rooting_comparison.tsv"
    manifest_out = outdir / "rooting_manifest.json"
    shutil.copy2(tree_path, unrooted_out)
    Phylo.write(rooted, str(rooted_out), "newick", format_branch_length="%1.10f")
    with mapping_out.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["iqtree_newick_tip", "input_tree_tip", "protein_accession", "taxon"])
        for row in rows:
            writer.writerow([iqtree_label(row["tree_tip"]), row["tree_tip"],
                             row["protein_accession"], row["taxon"]])

    measures = [
        ("tree_tips", len(all_tips)),
        ("designated_outgroups", ";".join(sorted(outgroup_accessions))),
        ("outgroup_split_present", "yes"),
        ("root_edge_support_label", root_edge_label or "not_labeled"),
        ("unrooted_nontrivial_splits", len(original_splits)),
        ("rooted_nontrivial_splits", len(nontrivial_splits(rooted))),
        ("unrooted_topology_identical_after_rooting", "yes"),
        ("max_pairwise_patristic_distance_difference", f"{max_distance_difference:.12g}"),
        ("internal_support_label_differences", len(label_differences)),
    ]
    with comparison_out.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["measure", "value"])
        writer.writerows(measures)
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "method": "Post-hoc reroot of the unrooted IQ-TREE topology using the unique edge separating both designated context taxa from the other 60 tips",
        "rooting_is_for_orientation_only": True,
        "outgroup_accessions": sorted(outgroup_accessions),
        "outgroup_tree_tips_original": sorted(outgroup_original_tips),
        "outgroup_tree_tips_iqtree_newick": sorted(outgroup_tips),
        "tip_label_map": str(mapping_out),
        "tip_label_map_sha256": sha256(mapping_out),
        "iqtree_newick_label_normalization": "replace '=' with '_'",
        "outgroup_edge_support_label": root_edge_label or None,
        "tree_input_manifest": str(manifest_path),
        "tree_input_manifest_sha256": sha256(manifest_path),
        "unrooted_iqtree_tree": str(tree_path),
        "unrooted_iqtree_tree_sha256": sha256(tree_path),
        "saved_unrooted_tree": str(unrooted_out),
        "saved_unrooted_tree_sha256": sha256(unrooted_out),
        "saved_rooted_tree": str(rooted_out),
        "saved_rooted_tree_sha256": sha256(rooted_out),
        "rooted_tree_written_with_biopython": __import__("Bio").__version__,
        "same_unrooted_splits": True,
        "max_pairwise_patristic_distance_difference": max_distance_difference,
        "support_label_differences": label_differences,
        "no_ancestral_state_or_directionality_analysis": True,
    }
    with manifest_out.open("x") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")

    print(f"Saved unrooted and rooted variants under {outdir}")
    print(f"Outgroup split present; maximum pairwise-distance difference: {max_distance_difference:.3g}")
    if label_differences:
        print(f"Warning: {len(label_differences)} internal support-label mapping(s) changed; inspect {manifest_out}")


if __name__ == "__main__":
    main()
