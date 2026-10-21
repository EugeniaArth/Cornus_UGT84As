#!/usr/bin/env Rscript

# Visualize the existing 62-tip MAFFT alignment and compare the PSPG box
# descriptively across UGT84 references, Group-L UGT74 comparators, context
# proteins, and Cornus candidates. No re-alignment or functional assignment is
# performed.

                                                                                                                            
                                                                     
                                                                                                                           
                                                         
                                                                         
     

suppressPackageStartupMessages({
  library(Biostrings)
  library(IRanges)
  library(ggmsa)
  library(ggplot2)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
alignment_file <- file.path(project_root, "Revised_analysis/Results/tree_62/mafft/final_alignment.faa")
tip_map_file <- file.path(project_root, "Revised_analysis/Results/tree_62/figures_v3/tree_plot_tip_label_map_detailed.tsv")
out_dir <- file.path(project_root, "Revised_analysis/Results/tree_62/pspg_analysis_62_group_EH")

required_files <- c(alignment_file, tip_map_file)
if (!all(file.exists(required_files))) {
  stop("Missing required input(s): ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}
if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Output directory is non-empty; refusing to overwrite: ", out_dir)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The 2024 Communications Chemistry article reports UGT84 PSPG landmarks at
# positions 2 (C/S), 5 (E), 15 (C), 40 (Q), and 44 (Q). In UGT84A13
# (AHA54051.1), these map to alignment columns 450-493, a 44-residue interval.
# The aligned HCGWNS-like core H[CS]G[WF]NS begins at column 468, i.e. PSPG
# position 19. These independent residue landmarks validate the box mapping.
pspg_box <- IRanges::IRanges(start = 450L, end = 493L)
pspg_length <- IRanges::width(pspg_box)
core_pattern <- "H[CS]G[WF]NS"
core_motif_position <- 19L

# Read and validate the exact tree alignment.
seqs <- Biostrings::readAAStringSet(alignment_file, format = "fasta", use.names = TRUE)
if (length(seqs) != 62L) stop("Expected 62 tree sequences; found ", length(seqs))
seq_lengths <- Biostrings::width(seqs)
if (length(unique(seq_lengths)) != 1L) stop("Alignment sequences do not have a common length")
if (max(seq_lengths) < IRanges::end(pspg_box)) stop("PSPG interval exceeds alignment length")

core_matches <- regexpr(core_pattern, as.character(seqs), perl = TRUE)
if (any(core_matches < 0L) || length(unique(core_matches)) != 1L) {
  stop("The HCGWNS-like PSPG core pattern is absent or not aligned at a shared column in all 62 sequences")
}
core_start <- as.integer(unique(core_matches))
if (core_start != IRanges::start(pspg_box) + core_motif_position - 1L) {
  stop("PSPG core-pattern coordinate does not match the configured PSPG-box anchor")
}

a13_accessions <- sub("\\|.*$", "", names(seqs))
a13_index <- match("AHA54051.1", a13_accessions)
if (is.na(a13_index)) stop("UGT84A13 reference AHA54051.1 is required to validate PSPG boundaries")
a13_box <- strsplit(as.character(Biostrings::subseq(seqs[a13_index], start = IRanges::start(pspg_box), end = IRanges::end(pspg_box))), "", fixed = TRUE)[[1]]
if (length(a13_box) != 44L) stop("UGT84A13 validation interval is not 44 residues")
landmark_positions <- c(2L, 5L, 15L, 40L, 44L)
landmark_residues <- a13_box[landmark_positions]
landmark_pass <- c(landmark_residues[1] %in% c("C", "S"), landmark_residues[2] == "E", landmark_residues[3] == "C", landmark_residues[4] == "Q", landmark_residues[5] == "Q")
if (!all(landmark_pass)) stop("UGT84A13 PSPG landmarks do not match the published residue positions; refusing to label this interval as the PSPG box")
boundary_validation <- data.frame(
  protein = "UGT84A13", protein_accession = "AHA54051.1",
  alignment_start = IRanges::start(pspg_box), alignment_end = IRanges::end(pspg_box),
  interval_width = pspg_length, box_sequence = paste0(a13_box, collapse = ""),
  pspg_position = landmark_positions, observed_residue = landmark_residues,
  published_ugt84_landmark = c("C or S", "E", "C", "Q", "Q"),
  landmark_match = landmark_pass, stringsAsFactors = FALSE
)

meta <- read.delim(tip_map_file, check.names = FALSE, stringsAsFactors = FALSE)
accessions <- sub("\\|.*$", "", names(seqs))
meta_index <- match(accessions, meta$protein_accession)
if (anyNA(meta_index) || anyDuplicated(meta$protein_accession[meta_index]) > 0L) {
  stop("Could not map every alignment sequence uniquely to the tree tip map")
}
meta <- meta[meta_index, , drop = FALSE]
meta$alignment_id <- names(seqs)
meta$protein_accession <- accessions
meta$sequence_group <- ifelse(
  meta$panel_role == "UGT84_reference", "UGT84 reference",
  ifelse(meta$panel_role == "Group_L_UGT74_comparator", "UGT74 Group-L comparator",
         ifelse(meta$panel_role == "family1_context_outgroup", "Group E/H",
                ifelse(meta$taxon == "Cornus florida", "C. florida candidate",
                       ifelse(meta$taxon == "Cornus suecica", "C. suecica candidate", "Other"))))
)
if (any(meta$sequence_group == "Other")) stop("Unexpected sequence role/taxon in the 62-tip map")

group_order <- c("UGT84 reference", "UGT74 Group-L comparator", "Group E/H", "C. florida candidate", "C. suecica candidate")
meta$sequence_group <- factor(meta$sequence_group, levels = group_order)
expected_group_counts <- c(24L, 6L, 2L, 20L, 10L)
observed_group_counts <- as.integer(table(factor(meta$sequence_group, levels = group_order)))
if (!identical(observed_group_counts, expected_group_counts)) {
  stop("Unexpected sequence-group counts; expected UGT84=24, UGT74=6, Group E/H=2, Florida=20, Suecica=10")
}
subfamily_order <- c(A = 1L, B = 2L, J = 3L, unassigned = 4L)
meta$subfamily_order <- unname(subfamily_order[ifelse(meta$UGT84_subfamily %in% names(subfamily_order), meta$UGT84_subfamily, "unassigned")])
meta$subfamily_order[is.na(meta$subfamily_order)] <- 4L
meta <- meta[order(meta$sequence_group, meta$subfamily_order, meta$gene_name, meta$protein_accession), , drop = FALSE]

short_label <- ifelse(
  meta$panel_role == "UGT84_reference",
  paste0("UGT84_", meta$UGT84_subfamily, "_", meta$gene_name),
  ifelse(meta$panel_role == "Group_L_UGT74_comparator",
         paste0("UGT74_L_", meta$gene_name),
         ifelse(meta$panel_role == "family1_context_outgroup",
                paste0("Context_", meta$gene_name),
                ifelse(meta$taxon == "Cornus florida",
                       paste0("Cfl_", meta$protein_accession),
                       paste0("Csue_", meta$protein_accession))))
)
short_label <- gsub("[^A-Za-z0-9_.-]", "_", short_label)
if (anyDuplicated(short_label)) short_label <- paste0(short_label, "_", meta$protein_accession)
meta$short_label <- short_label

# Reorder the same aligned sequences for readable group-wise plots.
seq_index <- match(meta$protein_accession, accessions)
ordered_seqs <- seqs[seq_index]
names(ordered_seqs) <- meta$short_label
grouped_fasta <- file.path(out_dir, "tree_62_alignment_grouped_short_headers.faa")
Biostrings::writeXStringSet(ordered_seqs, filepath = grouped_fasta, format = "fasta", width = 80L)

label_map <- meta[, c("short_label", "protein_accession", "alignment_id", "taxon", "panel_role", "UGT84_subfamily", "gene_name", "sequence_group")]
label_map$sequence_group <- as.character(label_map$sequence_group)
write.table(label_map, file.path(out_dir, "alignment_label_map.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(boundary_validation, file.path(out_dir, "pspg_box_boundary_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Extract the validated 44-column PSPG interval via IRanges/Biostrings.
pspg_seqs <- Biostrings::subseq(ordered_seqs, start = IRanges::start(pspg_box), end = IRanges::end(pspg_box))
pspg_fasta <- file.path(out_dir, "tree_62_pspg_box_44_alignment_columns.faa")
Biostrings::writeXStringSet(pspg_seqs, filepath = pspg_fasta, format = "fasta", width = 80L)

# Alignment overview and PSPG-box detail using ggmsa.
full_plot <- ggmsa::ggmsa(
  grouped_fasta, start = 1L, end = as.integer(seq_lengths[1]),
  font = NULL, seq_name = TRUE, by_conservation = TRUE, char_width = 0.55
) + ggplot2::labs(
  title = "MAFFT alignment used for the 62-tip tree",
  subtitle = "Rows ordered by panel role; color intensity reflects alignment conservation; no residues displayed"
) + ggplot2::theme(plot.title = element_text(face = "bold"))
ggplot2::ggsave(file.path(out_dir, "msa_62_full_alignment_overview.png"), full_plot, width = 22, height = 15, dpi = 160, bg = "white")

detail_plot <- ggmsa::ggmsa(
  grouped_fasta, start = IRanges::start(pspg_box), end = IRanges::end(pspg_box),
  font = "mono", color = "Chemistry_AA", seq_name = TRUE, char_width = 0.75,
  position_highlight = core_start:(core_start + 5L)
) + ggplot2::labs(
  title = "PSPG BOX — complete 44-residue interval",
  subtitle = "Alignment columns 450-493 mapped to PSPG positions 1-44 using UGT84A13 landmarks; HCGWNS-like core at positions 19-24 is highlighted",
  x = "PSPG-box residue position (alignment columns 450-493)",
  caption = "The displayed interval is the full PSPG box; highlighted columns show H[CS]G[WF]NS."
) + ggplot2::theme(plot.title = element_text(face = "bold"))
detail_plot <- suppressMessages(detail_plot + ggplot2::scale_x_continuous(
  breaks = c(450L, 454L, 458L, 462L, 466L, 468L, 473L, 478L, 482L, 486L, 490L, 493L),
  labels = c(1L, 5L, 9L, 13L, 17L, 19L, 24L, 29L, 33L, 37L, 41L, 44L),
  expand = c(0, 0)
))
ggplot2::ggsave(file.path(out_dir, "pspg_box_62_alignment_detail.png"), detail_plot, width = 21, height = 16, dpi = 190, bg = "white")

# Convert the box alignment to a character matrix; gaps are retained for
# coverage calculations and excluded only when estimating residue frequencies.
box_matrix <- do.call(rbind, lapply(as.character(pspg_seqs), function(s) strsplit(s, "", fixed = TRUE)[[1]]))
if (!identical(dim(box_matrix), c(62L, as.integer(pspg_length)))) stop("Unexpected PSPG-box matrix dimensions")
rownames(box_matrix) <- meta$short_label
colnames(box_matrix) <- as.character(seq_len(pspg_length))
meta$sequence_group <- factor(meta$sequence_group, levels = group_order)

profile_rows <- list()
frequency_rows <- list()
for (group_name in group_order) {
  row_idx <- which(as.character(meta$sequence_group) == group_name)
  for (pos in seq_len(pspg_length)) {
    vals <- box_matrix[row_idx, pos]
    nongap <- vals[!vals %in% c("-", ".", "?")]
    counts <- table(nongap)
    if (length(counts)) {
      max_n <- max(counts)
      major_candidates <- sort(names(counts)[counts == max_n])
      major <- major_candidates[1]
      major_n <- as.integer(counts[major])
      major_freq <- major_n / length(nongap)
      probs <- as.numeric(counts) / sum(counts)
      entropy <- -sum(probs * log2(probs))
      for (residue in names(counts)) {
        frequency_rows[[length(frequency_rows) + 1L]] <- data.frame(
          sequence_group = group_name, motif_position = pos,
          alignment_column = IRanges::start(pspg_box) + pos - 1L,
          residue = residue, count = as.integer(counts[residue]),
          n_non_gap = length(nongap), frequency = as.integer(counts[residue]) / length(nongap),
          stringsAsFactors = FALSE
        )
      }
    } else {
      major <- NA_character_; major_n <- 0L; major_freq <- NA_real_; entropy <- NA_real_
    }
    profile_rows[[length(profile_rows) + 1L]] <- data.frame(
      sequence_group = group_name, motif_position = pos,
      alignment_column = IRanges::start(pspg_box) + pos - 1L,
      n_sequences = length(row_idx), n_non_gap = length(nongap),
      n_gap = length(row_idx) - length(nongap),
      gap_fraction = (length(row_idx) - length(nongap)) / length(row_idx),
      major_residue = major, major_count = major_n,
      major_frequency = major_freq, entropy_bits = entropy,
      stringsAsFactors = FALSE
    )
  }
}
profiles <- do.call(rbind, profile_rows)
residue_frequencies <- do.call(rbind, frequency_rows)
write.table(profiles, file.path(out_dir, "pspg_box_group_site_profiles.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(residue_frequencies, file.path(out_dir, "pspg_box_group_residue_frequencies.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Wide comparison table: UGT84 reference profile against the sampled Group-L
# UGT74 comparators, context proteins, and each Cornus candidate set.
group_keys <- c(
  "UGT84 reference" = "ugt84",
  "UGT74 Group-L comparator" = "ugt74_groupL",
  "Group E/H" = "group_EH",
  "C. florida candidate" = "cornus_florida",
  "C. suecica candidate" = "cornus_suecica"
)
if (!identical(names(group_keys), group_order)) {
  stop("Group-profile keys must exactly match the sequence-group labels and order")
}
site_summary <- data.frame(
  motif_position = seq_len(pspg_length),
  alignment_column = IRanges::start(pspg_box) + seq_len(pspg_length) - 1L,
  stringsAsFactors = FALSE
)
for (group_name in names(group_keys)) {
  key <- group_keys[[group_name]]
  sub <- profiles[profiles$sequence_group == group_name, , drop = FALSE]
  sub <- sub[order(sub$motif_position), , drop = FALSE]
  site_summary[[paste0(key, "_n_sequences")]] <- sub$n_sequences
  site_summary[[paste0(key, "_n_non_gap")]] <- sub$n_non_gap
  site_summary[[paste0(key, "_gap_fraction")]] <- sub$gap_fraction
  site_summary[[paste0(key, "_major_residue")]] <- sub$major_residue
  site_summary[[paste0(key, "_major_frequency")]] <- sub$major_frequency
}

freq_of_residue <- function(group_name, pos, residue) {
  x <- residue_frequencies[residue_frequencies$sequence_group == group_name & residue_frequencies$motif_position == pos & residue_frequencies$residue == residue, , drop = FALSE]
  if (nrow(x) == 0L) return(0)
  x$frequency[1]
}
site_summary$ugt84_residue_frequency_in_ugt74 <- mapply(
  function(pos, aa) if (is.na(aa)) NA_real_ else freq_of_residue("UGT74 Group-L comparator", pos, aa),
  site_summary$motif_position, site_summary$ugt84_major_residue
)
site_summary$ugt84_minus_ugt74_frequency <- site_summary$ugt84_major_frequency - site_summary$ugt84_residue_frequency_in_ugt74
site_summary$provisional_ugt84_enriched <- !is.na(site_summary$ugt84_major_frequency) &
  site_summary$ugt84_n_non_gap >= 20L & site_summary$ugt74_groupL_n_non_gap >= 5L &
  site_summary$ugt84_major_frequency >= 0.80 & site_summary$ugt84_residue_frequency_in_ugt74 <= (2 / 6)
write.table(site_summary, file.path(out_dir, "pspg_box_site_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

selected_positions <- which(site_summary$provisional_ugt84_enriched)
selected_summary <- site_summary[site_summary$provisional_ugt84_enriched, , drop = FALSE]
write.table(selected_summary, file.path(out_dir, "provisional_ugt84_enriched_positions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Candidate-level profile scores are descriptive matches to the positions
# enriched in this small sampled reference/comparator panel.
candidate_rows <- which(meta$panel_role == "Cornus_candidate")
candidate_scores <- list()
candidate_position_rows <- list()
for (i in candidate_rows) {
  group_label <- as.character(meta$sequence_group[i])
  if (length(selected_positions) == 0L) {
    candidate_scores[[length(candidate_scores) + 1L]] <- data.frame(
      short_label = meta$short_label[i], protein_accession = meta$protein_accession[i],
      taxon = meta$taxon[i], gene_name = meta$gene_name[i], sequence_group = group_label,
      n_selected_positions = 0L, n_non_gap = 0L, coverage_pct = NA_real_,
      ugt84_match_n = NA_integer_, ugt84_match_pct_all = NA_real_, ugt84_match_pct_covered = NA_real_,
      ugt74_groupL_match_n = NA_integer_, ugt74_groupL_match_pct_all = NA_real_,
      ugt84_minus_ugt74_match_pct_points = NA_real_, stringsAsFactors = FALSE
    )
    next
  }
  residues <- box_matrix[i, selected_positions]
  nongap_idx <- which(!residues %in% c("-", ".", "?"))
  ugt84_consensus <- site_summary$ugt84_major_residue[selected_positions]
  ugt74_consensus <- site_summary$ugt74_groupL_major_residue[selected_positions]
  match84 <- !residues %in% c("-", ".", "?") & residues == ugt84_consensus
  match74 <- !residues %in% c("-", ".", "?") & residues == ugt74_consensus
  n84 <- sum(match84, na.rm = TRUE)
  n74 <- sum(match74, na.rm = TRUE)
  candidate_scores[[length(candidate_scores) + 1L]] <- data.frame(
    short_label = meta$short_label[i], protein_accession = meta$protein_accession[i],
    taxon = meta$taxon[i], gene_name = meta$gene_name[i], sequence_group = group_label,
    n_selected_positions = length(selected_positions), n_non_gap = length(nongap_idx),
    coverage_pct = 100 * length(nongap_idx) / length(selected_positions),
    ugt84_match_n = n84, ugt84_match_pct_all = 100 * n84 / length(selected_positions),
    ugt84_match_pct_covered = if (length(nongap_idx)) 100 * n84 / length(nongap_idx) else NA_real_,
    ugt74_groupL_match_n = n74, ugt74_groupL_match_pct_all = 100 * n74 / length(selected_positions),
    ugt84_minus_ugt74_match_pct_points = 100 * (n84 - n74) / length(selected_positions),
    stringsAsFactors = FALSE
  )
  for (j in seq_along(selected_positions)) {
    pos <- selected_positions[j]
    candidate_position_rows[[length(candidate_position_rows) + 1L]] <- data.frame(
      short_label = meta$short_label[i], protein_accession = meta$protein_accession[i],
      taxon = meta$taxon[i], gene_name = meta$gene_name[i], sequence_group = group_label,
      motif_position = pos, alignment_column = site_summary$alignment_column[pos],
      candidate_residue = residues[j], ugt84_major_residue = ugt84_consensus[j],
      ugt84_major_frequency = site_summary$ugt84_major_frequency[pos],
      ugt74_groupL_major_residue = ugt74_consensus[j],
      ugt74_frequency_of_ugt84_residue = site_summary$ugt84_residue_frequency_in_ugt74[pos],
      matches_ugt84_major = if (residues[j] %in% c("-", ".", "?")) NA else residues[j] == ugt84_consensus[j],
      matches_ugt74_major = if (residues[j] %in% c("-", ".", "?")) NA else residues[j] == ugt74_consensus[j],
      stringsAsFactors = FALSE
    )
  }
}
candidate_scores <- do.call(rbind, candidate_scores)
candidate_positions <- if (length(candidate_position_rows)) do.call(rbind, candidate_position_rows) else data.frame()
write.table(candidate_scores, file.path(out_dir, "cornus_pspg_profile_matches.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(candidate_positions, file.path(out_dir, "cornus_pspg_position_matches.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Leave-one-UGT84-reference-out stability check. The comparator profile remains
# the six sampled Group-L UGT74 proteins in each fold.
ugt84_rows <- which(meta$panel_role == "UGT84_reference")
ugt74_rows <- which(meta$panel_role == "Group_L_UGT74_comparator")
loo_rows <- list()
for (heldout in ugt84_rows) {
  train84 <- setdiff(ugt84_rows, heldout)
  selected <- integer(0)
  consensus84 <- character(0)
  for (pos in seq_len(pspg_length)) {
    vals84 <- box_matrix[train84, pos]
    vals84 <- vals84[!vals84 %in% c("-", ".", "?")]
    vals74 <- box_matrix[ugt74_rows, pos]
    vals74 <- vals74[!vals74 %in% c("-", ".", "?")]
    if (length(vals84) < 18L || length(vals74) < 5L) next
    tab84 <- table(vals84)
    major <- sort(names(tab84)[tab84 == max(tab84)])[1]
    f84 <- as.integer(tab84[major]) / length(vals84)
    f74 <- sum(vals74 == major) / length(vals74)
    if (f84 >= 0.80 && f74 <= (2 / 6)) {
      selected <- c(selected, pos)
      consensus84 <- c(consensus84, major)
    }
  }
  heldout_residues <- box_matrix[heldout, selected]
  covered <- !heldout_residues %in% c("-", ".", "?")
  match84 <- sum(covered & heldout_residues == consensus84)
  cons74 <- vapply(selected, function(pos) {
    vals <- box_matrix[ugt74_rows, pos]
    vals <- vals[!vals %in% c("-", ".", "?")]
    if (!length(vals)) return(NA_character_)
    tab <- table(vals); sort(names(tab)[tab == max(tab)])[1]
  }, character(1))
  match74 <- sum(covered & heldout_residues == cons74, na.rm = TRUE)
  loo_rows[[length(loo_rows) + 1L]] <- data.frame(
    heldout_accession = meta$protein_accession[heldout], heldout_gene = meta$gene_name[heldout],
    heldout_subfamily = meta$UGT84_subfamily[heldout],
    n_training_ugt84 = length(train84), n_ugt74_comparators = length(ugt74_rows),
    n_selected_positions = length(selected), n_covered_positions = sum(covered),
    ugt84_match_n = match84,
    ugt84_match_pct_covered = if (sum(covered)) 100 * match84 / sum(covered) else NA_real_,
    ugt74_match_n = match74,
    ugt74_match_pct_covered = if (sum(covered)) 100 * match74 / sum(covered) else NA_real_,
    stringsAsFactors = FALSE
  )
}
loo <- do.call(rbind, loo_rows)
write.table(loo, file.path(out_dir, "ugt84_reference_leave_one_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Group profile plot: each cell is the most common amino acid and its count;
# background intensity is its frequency among non-gap residues.
plot_profiles <- profiles
plot_profiles$sequence_group <- factor(plot_profiles$sequence_group, levels = group_order)
plot_profiles$cell_label <- ifelse(
  is.na(plot_profiles$major_residue), "-",
  paste0(plot_profiles$major_residue, "\n", plot_profiles$major_count, "/", plot_profiles$n_non_gap)
)
profile_plot <- ggplot2::ggplot(plot_profiles, ggplot2::aes(x = motif_position, y = sequence_group, fill = major_frequency)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.25) +
  ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.5, lineheight = 0.85) +
  ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90", name = "Major-residue\nfrequency") +
  ggplot2::scale_x_continuous(breaks = seq(1, pspg_length, by = 2), expand = c(0, 0)) +
  ggplot2::labs(
    title = "PSPG-box residue profiles in the 62-sequence panel",
    subtitle = "Cell text: dominant residue and count/non-gap sequences; group sizes are 24, 6, 2, 20, and 10",
    x = "PSPG-box position (44-column alignment interval)", y = NULL,
    caption = "UGT74 proteins are Group-L comparators. Frequencies describe this sampled panel only."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(panel.grid = element_blank(), axis.text.x = element_text(size = 8),
                 axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold"),
                 legend.position = "bottom")
ggplot2::ggsave(file.path(out_dir, "pspg_box_group_residue_profiles.png"), profile_plot, width = 17, height = 6.2, dpi = 220, bg = "white")

if (length(selected_positions) > 0L && any(!is.na(candidate_scores$ugt84_match_pct_all))) {
  candidate_scores$plot_label <- paste0(candidate_scores$short_label, " [", ifelse(candidate_scores$sequence_group == "C. florida candidate", "Florida", "Suecica"), "]")
  candidate_scores$plot_label <- factor(candidate_scores$plot_label, levels = candidate_scores$plot_label[order(candidate_scores$ugt84_match_pct_all)])
  candidate_plot <- ggplot2::ggplot(candidate_scores, ggplot2::aes(x = plot_label, y = ugt84_match_pct_all, fill = sequence_group)) +
    ggplot2::geom_col(width = 0.78) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = c("C. florida candidate" = "#2D7790", "C. suecica candidate" = "#C5773E"), name = "Candidate set") +
    ggplot2::labs(
      title = "Cornus candidates: match to provisional UGT84-enriched PSPG residues",
      subtitle = paste0("", length(selected_positions), " positions selected from 24 UGT84 references vs 6 Group-L UGT74 comparators; gaps count as non-matches"),
      x = NULL, y = "Matches to UGT84 major residue (%)",
      caption = "Profile match is sequence evidence only; it does not establish UGT84 membership or enzyme activity."
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold"),
                   legend.position = "bottom", axis.text.y = element_text(size = 7))
  ggplot2::ggsave(file.path(out_dir, "cornus_pspg_profile_match.png"), candidate_plot, width = 15, height = 11, dpi = 220, bg = "white")
} else {
  no_profile_plot <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = "No positions met the stated exploratory UGT84-enrichment threshold.") +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) + ggplot2::theme_void() +
    ggplot2::labs(title = "Cornus candidate PSPG comparison")
  ggplot2::ggsave(file.path(out_dir, "cornus_pspg_profile_match.png"), no_profile_plot, width = 9, height = 3.5, dpi = 180, bg = "white")
}

# Run manifest and concise human-readable report.
role_counts <- as.data.frame(table(meta$sequence_group), stringsAsFactors = FALSE)
colnames(role_counts) <- c("sequence_group", "n")
role_lines <- paste0("- ", role_counts$sequence_group, ": n = ", role_counts$n)
selected_text <- if (nrow(selected_summary)) {
  paste(apply(selected_summary[, c("motif_position", "alignment_column", "ugt84_major_residue", "ugt84_major_frequency", "ugt84_residue_frequency_in_ugt74")], 1, function(x) {
    sprintf("%s (alignment column %s): UGT84 %s at %.0f%%; same residue in UGT74 Group-L at %.0f%%", x[1], x[2], x[3], 100 * as.numeric(x[4]), 100 * as.numeric(x[5]))
  }), collapse = "\n")
} else {
  "- None met the prespecified descriptive threshold (UGT84 major-residue frequency >= 0.80 and frequency of that residue in Group-L UGT74 <= 2/6)."
}
loo_positive <- sum(loo$ugt84_match_pct_covered > loo$ugt74_match_pct_covered, na.rm = TRUE)
report_lines <- c(
  "# PSPG-box analysis of the 62-sequence tree alignment",
  "",
  "## Scope and inputs",
  "",
  paste0("- Alignment: `", file.path("Revised_analysis/Results/tree_62/mafft/final_alignment.faa"), "` (62 sequences; ", seq_lengths[1], " alignment columns; existing MAFFT output; no re-alignment)."),
  paste0("- Tip metadata: `", file.path("Revised_analysis/Results/tree_62/figures_v3/tree_plot_tip_label_map_detailed.tsv"), "`."),
  "- Visualization/extraction: R packages `ggmsa`, `Biostrings`, and `IRanges`; profiles and figures use `ggplot2`.",
  "- Included groups:", role_lines,
  "",
  "## PSPG-box interval and caveats",
  "",
  paste0("- Boundary evidence: Li et al. (2024), Communications Chemistry, Figure 2a (https://doi.org/10.1038/s42004-024-01231-1), describe a 44-aa PSPG box and UGT84 landmarks at positions 2 (C/S), 5 (E), 15 (C), 40 (Q), and 44 (Q)."),
  paste0("- In UGT84A13 reference AHA54051.1, alignment columns ", IRanges::start(pspg_box), "-", IRanges::end(pspg_box), " give the 44-aa sequence `", paste0(a13_box, collapse = ""), "`. Residues at positions 2, 5, 15, 40, and 44 are ", paste(landmark_residues, collapse = ", "), ", matching those published landmarks. The HCGWNS-like pattern begins at alignment column ", core_start, " and occupies PSPG positions 19-24; exact `HCGWNS` occurs in ", sum(grepl("HCGWNS", as.character(seqs), fixed = TRUE)), "/62, with six aligned variants."),
  "- The earlier 453-496 interval was offset three residues downstream and did not satisfy these landmarks; its profiles and candidate scores are superseded. The corrected interval is columns 450-493. Gaps are retained and reported in the corrected tables.",
  "",
  "## Provisional UGT84-enriched positions",
  "",
  "Threshold: the major UGT84 residue occurs in at least 80% of non-gap UGT84 references (n=24), and occurs in no more than 2 of the 6 Group-L UGT74 comparators. This is a descriptive filter for this panel, not a validated family classifier.",
  selected_text,
  "",
  "## Candidate-profile comparison",
  "",
  if (length(selected_positions)) {
    paste0("- The candidate table scores the 30 Cornus sequences in the tree against ", length(selected_positions), " selected positions. Gaps are counted as non-matches in the all-position score; coverage-adjusted scores are also reported." )
  } else {
    "- No signature positions passed the threshold; candidate matching scores are not estimable under this rule."
  },
  "- Candidate scores are comparisons to the UGT84-major-residue profile and, in parallel, the Group-L UGT74-major-residue profile. They do not assign Cornus candidates to UGT84, UGT74, a UGT84 subfamily, or a biochemical function.",
  "",
  "## Leave-one-reference-out check",
  "",
  paste0("- For each of the 24 UGT84 references, positions were reselected using the other 23 references and the six Group-L UGT74 comparators. In ", loo_positive, " of 24 held-out references, the score against the training UGT84 profile exceeded the score against the UGT74 profile. See `ugt84_reference_leave_one_out.tsv` for fold-specific position counts, coverage, and scores."),
  "- This check is exploratory, uses the same small comparator panel, and does not substitute for broader non-UGT84 controls or independent experimental validation.",
  "",
  "## Interpretation",
  "",
  "The Nature article's Figure 2a compares selected UGT72 and UGT84 sequences; this panel instead contains 24 UGT84 references, six Group-L UGT74 comparators, two Group E/H context proteins, and 30 Cornus candidates. Different group composition and profile summaries need not produce identical residue frequencies. The corrected UGT84A13 landmark check supports the box mapping, but this limited panel cannot establish a general UGT84-specific signature. PSPG similarity, candidate scores, tree placement, and screening results do not demonstrate enzyme function.",
  "",
  "## Outputs",
  "",
  "- `msa_62_full_alignment_overview.png` — full-alignment conservation overview.",
  "- `pspg_box_62_alignment_detail.png` — explicitly labeled complete PSPG box; x-axis shows positions 1-44 and the HCGWNS-like core is highlighted.",
  "- `pspg_box_boundary_validation.tsv` — UGT84A13 landmark-residue check supporting the 450-493 interval.",
  "- `pspg_box_group_residue_profiles.png` — group-wise dominant residue and frequency by position.",
  "- `cornus_pspg_profile_match.png` — candidate-level descriptive profile scores.",
  "- `pspg_box_site_comparison.tsv` and `pspg_box_group_residue_frequencies.tsv` — full site-level profiles and residue frequencies.",
  "- `cornus_pspg_profile_matches.tsv` and `cornus_pspg_position_matches.tsv` — candidate scores and per-position details.",
  "- `provisional_ugt84_enriched_positions.tsv` and `ugt84_reference_leave_one_out.tsv` — exploratory positions and reference holdout check.",
  "- `tree_62_alignment_grouped_short_headers.faa`, `tree_62_pspg_box_44_alignment_columns.faa`, and `alignment_label_map.tsv` — extracted alignment data and readable label mapping.",
  "",
  "## Reproducibility",
  "",
  paste0("- Script: `Revised_analysis/Scripts/18_visualize_pspg_62.R`."),
  "- Run from the project root with: `Rscript Revised_analysis/Scripts/18_visualize_pspg_62.R`.",
  paste0("- R versions: ggmsa ", as.character(packageVersion("ggmsa")), "; IRanges ", as.character(packageVersion("IRanges")), "; Biostrings ", as.character(packageVersion("Biostrings")), "; ggplot2 ", as.character(packageVersion("ggplot2")), "."),
  ""
)
writeLines(report_lines, file.path(out_dir, "pspg_analysis_62_report.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

input_hashes <- tools::md5sum(c(alignment_file, tip_map_file))
manifest <- data.frame(
  item = c("alignment", "tip_map"),
  path = substring(names(input_hashes), nchar(project_root) + 2L),
  md5 = unname(input_hashes),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(out_dir, "input_manifest_md5.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("PSPG analysis complete.\n")
cat("Sequences:", length(seqs), " Alignment columns:", seq_lengths[1], " PSPG interval:", IRanges::start(pspg_box), "-", IRanges::end(pspg_box), "\n")
cat("Core pattern:", core_pattern, "at column", core_start, " Exact HCGWNS:", sum(grepl("HCGWNS", as.character(seqs), fixed = TRUE)), "/", length(seqs), "\n")
cat("Provisional UGT84-enriched positions:", length(selected_positions), "\n")
cat("Output directory:", out_dir, "\n")
