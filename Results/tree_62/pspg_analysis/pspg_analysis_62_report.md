# PSPG-box analysis of the 62-sequence tree alignment

## Scope and inputs

- Alignment: `Revised_analysis/Results/tree_62/mafft/final_alignment.faa` (62 sequences; 685 alignment columns; existing MAFFT output; no re-alignment).
- Tip metadata: `Revised_analysis/Results/tree_62/figures_v3/tree_plot_tip_label_map_detailed.tsv`.
- Visualization/extraction: R packages `ggmsa`, `Biostrings`, and `IRanges`; profiles and figures use `ggplot2`.
- Included groups:
- UGT84 reference: n = 24
- UGT74 Group-L comparator: n = 6
- Group E/H: n = 2
- C. florida candidate: n = 20
- C. suecica candidate: n = 10

## PSPG-box interval and caveats

- Boundary evidence: Li et al. (2024), Communications Chemistry, Figure 2a (https://doi.org/10.1038/s42004-024-01231-1), describe a 44-aa PSPG box and UGT84 landmarks at positions 2 (C/S), 5 (E), 15 (C), 40 (Q), and 44 (Q).
- In UGT84A13 reference AHA54051.1, alignment columns 450-493 give the 44-aa sequence `WSPQEQVLAHPSVACFVTHCGWNSTMESLTSGMPVVAFPQWGDQ`. Residues at positions 2, 5, 15, 40, and 44 are S, E, C, Q, Q, matching those published landmarks. The HCGWNS-like pattern begins at alignment column 468 and occupies PSPG positions 19-24; exact `HCGWNS` occurs in 56/62, with six aligned variants.
- The earlier 453-496 interval was offset three residues downstream and did not satisfy these landmarks; its profiles and candidate scores are superseded. The corrected interval is columns 450-493. Gaps are retained and reported in the corrected tables.

## Provisional UGT84-enriched positions

Threshold: the major UGT84 residue occurs in at least 80% of non-gap UGT84 references (n=24), and occurs in no more than 2 of the 6 Group-L UGT74 comparators. This is a descriptive filter for this panel, not a validated family classifier.
 5 (alignment column 454): UGT84 E at 83%; same residue in UGT74 Group-L at 0%
11 (alignment column 460): UGT84 P at 92%; same residue in UGT74 Group-L at 0%
14 (alignment column 463): UGT84 A at 83%; same residue in UGT74 Group-L at 17%
31 (alignment column 480): UGT84 S at 83%; same residue in UGT74 Group-L at 0%

## Candidate-profile comparison

- The candidate table scores the 30 Cornus sequences in the tree against 4 selected positions. Gaps are counted as non-matches in the all-position score; coverage-adjusted scores are also reported.
- Candidate scores are comparisons to the UGT84-major-residue profile and, in parallel, the Group-L UGT74-major-residue profile. They do not assign Cornus candidates to UGT84, UGT74, a UGT84 subfamily, or a biochemical function.

## Leave-one-reference-out check

- For each of the 24 UGT84 references, positions were reselected using the other 23 references and the six Group-L UGT74 comparators. In 23 of 24 held-out references, the score against the training UGT84 profile exceeded the score against the UGT74 profile. See `ugt84_reference_leave_one_out.tsv` for fold-specific position counts, coverage, and scores.
- This check is exploratory, uses the same small comparator panel, and does not substitute for broader non-UGT84 controls or independent experimental validation.

## Interpretation

The Nature article's Figure 2a compares selected UGT72 and UGT84 sequences; this panel instead contains 24 UGT84 references, six Group-L UGT74 comparators, two Group E/H context proteins, and 30 Cornus candidates. Different group composition and profile summaries need not produce identical residue frequencies. The corrected UGT84A13 landmark check supports the box mapping, but this limited panel cannot establish a general UGT84-specific signature. PSPG similarity, candidate scores, tree placement, and screening results do not demonstrate enzyme function.

## Outputs

- `msa_62_full_alignment_overview.png` — full-alignment conservation overview.
- `pspg_box_62_alignment_detail.png` — explicitly labeled complete PSPG box; x-axis shows positions 1-44 and the HCGWNS-like core is highlighted.
- `pspg_box_boundary_validation.tsv` — UGT84A13 landmark-residue check supporting the 450-493 interval.
- `pspg_box_group_residue_profiles.png` — group-wise dominant residue and frequency by position.
- `cornus_pspg_profile_match.png` — candidate-level descriptive profile scores.
- `pspg_box_site_comparison.tsv` and `pspg_box_group_residue_frequencies.tsv` — full site-level profiles and residue frequencies.
- `cornus_pspg_profile_matches.tsv` and `cornus_pspg_position_matches.tsv` — candidate scores and per-position details.
- `provisional_ugt84_enriched_positions.tsv` and `ugt84_reference_leave_one_out.tsv` — exploratory positions and reference holdout check.
- `tree_62_alignment_grouped_short_headers.faa`, `tree_62_pspg_box_44_alignment_columns.faa`, and `alignment_label_map.tsv` — extracted alignment data and readable label mapping.

## Reproducibility

- Script: `Revised_analysis/Scripts/18_visualize_pspg_62.R`.
- Run from the project root with: `Rscript Revised_analysis/Scripts/18_visualize_pspg_62.R`.
- R versions: ggmsa 1.17.0; IRanges 2.44.0; Biostrings 2.78.0; ggplot2 4.0.3.

