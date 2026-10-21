#!/usr/bin/env Rscript
# Create README count charts and a 16:9 pipeline diagram using base R graphics.
# Run from Revised_analysis/: Rscript Scripts/17_create_readme_figures.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  outdir <- file.path("Results", "readme_figures_v2")
} else if (length(args) == 2 && args[[1]] == "--outdir") {
  outdir <- args[[2]]
} else {
  stop("Usage: Rscript Scripts/17_create_readme_figures.R [--outdir NEW_DIRECTORY]")
}

if (dir.exists(outdir) || file.exists(outdir)) {
  stop(sprintf("Refusing to overwrite existing output path: %s", outdir))
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing required input: %s", path))
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE,
             quote = "", comment.char = "")
}

panel <- read_tsv(file.path("Results", "ugt84_reference_context_panel_manifest.tsv"))
candidates <- read_tsv(file.path("Results", "cornus_comparison", "final_input",
                                  "combined_cornus_reference_candidates_65_manifest.tsv"))
tree_input <- read_tsv(file.path("Results", "tree_62", "input",
                                 "cornus_ugt84_tree_input_62_manifest.tsv"))
held_out <- read_tsv(file.path("Results", "tree_62", "input",
                              "excluded_suecica_start_flagged_candidates.tsv"))
fl_core <- read_tsv(file.path("Results", "cornus_florida", "pspg_core",
                              "pspg_core_check.tsv"))
su_core <- read_tsv(file.path("Results", "cornus_suecica", "pspg_core_canonical_ids",
                              "pspg_core_check.tsv"))

stopifnot(nrow(panel) == 32, nrow(candidates) == 65, nrow(tree_input) == 62,
          nrow(held_out) == 3)

ref_counts <- c(
  "UGT84-A" = sum(panel$panel_role == "UGT84_reference" & panel$subfamily == "A"),
  "UGT84-B" = sum(panel$panel_role == "UGT84_reference" & panel$subfamily == "B"),
  "UGT84-J" = sum(panel$panel_role == "UGT84_reference" & panel$subfamily == "J"),
  "UGT74 / Group-L comparator" = sum(panel$panel_role == "Group_L_UGT74_comparator"),
  "Group E context" = sum(panel$panel_role == "family1_context_outgroup" & panel$phylogenetic_group == "E"),
  "Group H context" = sum(panel$panel_role == "family1_context_outgroup" & panel$phylogenetic_group == "H")
)
if (sum(ref_counts) != 32) stop("Reference/context composition does not sum to 32")

species <- c("Cornus florida", "Cornus suecica")
selected <- vapply(species, function(s) sum(candidates$taxon == s & candidates$panel_role == "Cornus_candidate"), numeric(1))
in_tree <- vapply(species, function(s) sum(tree_input$taxon == s & tree_input$panel_role == "Cornus_candidate"), numeric(1))
held <- vapply(species, function(s) sum(held_out$taxon == s & held_out$panel_role == "Cornus_candidate"), numeric(1))
if (any(selected != in_tree + held)) stop("Selected candidate counts do not equal tree plus holdout counts")

core_yes <- function(x) sum(x$exact_HCGWNS_core_present == "yes")
screen_counts <- data.frame(
  dataset = c("C. florida", "C. suecica"),
  PF00201_hits = c(nrow(fl_core), nrow(su_core)),
  exact_HCGWNS_core = c(core_yes(fl_core), core_yes(su_core)),
  selected_loci = as.numeric(selected),
  tree_tips = as.numeric(in_tree),
  held_out_from_tree = as.numeric(held),
  stringsAsFactors = FALSE
)

# Save the exact values used to build the figures alongside them.
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figure_counts <- rbind(
  data.frame(dataset = "Reference/context panel", category = names(ref_counts), count = as.integer(ref_counts)),
  data.frame(dataset = "Cornus candidates", category = paste0(species, " selected loci"), count = as.integer(selected)),
  data.frame(dataset = "Cornus candidates", category = paste0(species, " tree tips"), count = as.integer(in_tree)),
  data.frame(dataset = "Cornus candidates", category = paste0(species, " held out from tree"), count = as.integer(held)),
  data.frame(dataset = "PF00201/core screen", category = paste0(screen_counts$dataset, " PF00201 hits"), count = screen_counts$PF00201_hits),
  data.frame(dataset = "PF00201/core screen", category = paste0(screen_counts$dataset, " exact HCGWNS"), count = screen_counts$exact_HCGWNS_core)
)
write.table(figure_counts, file.path(outdir, "readme_figure_counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Figure 1: composition of the 32-protein reference/context panel.
panel_colors <- c("#6A3D9A", "#E69F00", "#009E73", "#65717C", "#9AA7B1", "#C4CCD2")
png(file.path(outdir, "reference_panel_composition.png"), width = 1800, height = 1150,
    res = 180, type = "quartz")
par(mar = c(5, 13, 4, 2), family = "sans", xaxs = "i")
mid <- barplot(rev(ref_counts), horiz = TRUE, names.arg = rev(names(ref_counts)),
               las = 1, col = rev(panel_colors), border = NA,
               xlim = c(0, max(ref_counts) + 4), xlab = "Protein sequences",
               main = "Reference/context panel composition (n = 32)",
               cex.names = 0.9, cex.axis = 0.85, cex.lab = 0.95)
text(rev(ref_counts) + 0.35, mid, labels = rev(as.integer(ref_counts)),
     adj = c(0, 0.5), cex = 0.9, col = "#333333")
box(col = "#D7DCE0")
dev.off()

# Figure 2: selected Cornus loci partitioned by tree inclusion/holdout status.
status_counts <- rbind("Included in 62-tip tree" = in_tree,
                       "Held out from tree" = held)
png(file.path(outdir, "cornus_candidate_tree_counts.png"), width = 1500, height = 1050,
    res = 180, type = "quartz")
par(mar = c(5, 5, 4, 2), family = "sans", yaxs = "i")
bp <- barplot(status_counts, beside = FALSE, names.arg = c("C. florida", "C. suecica"),
              col = c("#4C78A8", "#E7A6C9"), border = "white",
              ylim = c(0, max(selected) + 5), ylab = "Locus representatives",
              main = "Cornus candidate counts and tree inclusion",
              cex.names = 1, cex.axis = 0.9, cex.lab = 0.95)
text(bp, selected + 0.8, labels = paste0("selected n = ", selected), cex = 0.9)
legend("topright", legend = rownames(status_counts), fill = c("#4C78A8", "#E7A6C9"),
       border = NA, bty = "n", cex = 0.85)
box(col = "#D7DCE0")
dev.off()

# Figure 3: 16:9 workflow overview for PowerPoint. This base-R rendering mirrors
# the Mermaid flowchart included in README.md; no diagram/rendering package is needed.
draw_box <- function(x, y, w, h, lines, fill = "#F5F7FA", border = "#66727D", cex = 0.76) {
  rect(x, y, x + w, y + h, col = fill, border = border, lwd = 1.35)
  text(x + w / 2, y + h * 0.78, lines[1], font = 2, cex = cex + 0.04, col = "#26323B")
  if (length(lines) > 1) {
    ys <- seq(y + h * 0.57, y + h * 0.25, length.out = length(lines) - 1)
    text(rep(x + w / 2, length(ys)), ys, lines[-1], cex = cex, col = "#34424D")
  }
}

draw_arrow <- function(x0, y0, x1, y1, dashed = FALSE, color = "#52616B") {
  arrows(x0, y0, x1, y1, length = 0.10, angle = 22, code = 2,
         lty = if (dashed) 2 else 1, lwd = 1.25, col = color)
}

fl_text <- c("C. FLORIDA", "RefSeq proteins", "PF00201: 257 hits; 160 HCGWNS", "BLASTP + GFF mapping: 20 loci")
su_text <- c("C. SUECICA", "Genome + GFF3/GTF", "33,025 transcripts; QC passed", "PF00201: 146 hits; 100 HCGWNS", "BLASTP + GFF: 13 loci; 3 flagged")
panel_text <- c("REFERENCE / CONTEXT PANEL", "32 proteins", "24 UGT84 A/B/J", "6 UGT74 Group-L comparators", "2 context sequences (E/H)")
candidate_text <- c("CANDIDATE DATASET", "65 proteins", "32 panel + 20 Florida", "+ 13 Suecica")
tree_text <- c("TREE-ONLY INPUT", "62 tips", "20 Florida + 10 Suecica", "3 flagged Suecica held out")
inference_text <- c("PHYLOGENY", "MAFFT L-INS-i", "IQ-TREE MFP", "1,000 UFBoot + SH-aLRT", "Unrooted topology")

png(file.path(outdir, "pipeline_overview_ppt.png"), width = 2400, height = 1350,
    res = 180, type = "quartz")
par(mar = c(0, 0, 0, 0), family = "sans", xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 16), ylim = c(0, 9), asp = 1)
text(8, 8.82, "Cornus UGT84 candidate-screening and phylogeny workflow",
     cex = 1.35, font = 2, col = "#202A32")

draw_box(0.35, 6.45, 3.55, 2.0, panel_text, fill = "#ECEEF8", border = "#6A3D9A", cex = 0.74)
draw_box(4.65, 6.2, 3.55, 2.05, fl_text, fill = "#E7F1F8", border = "#0072B2", cex = 0.72)
draw_box(4.65, 3.25, 3.55, 2.35, su_text, fill = "#F8EAF2", border = "#CC79A7", cex = 0.70)
draw_box(8.65, 5.25, 2.45, 1.9, candidate_text, fill = "#F1F3F5", cex = 0.74)
draw_box(11.55, 5.25, 2.35, 1.9, tree_text, fill = "#FFF4E5", border = "#C78E2B", cex = 0.70)
draw_box(14.2, 5.25, 1.5, 1.9, inference_text, fill = "#E9F4EC", border = "#4E8B5C", cex = 0.63)
draw_box(11.2, 2.15, 4.5, 1.55,
         c("OUTPUTS", "Unrooted ML tree + detailed figure", "16:9 PowerPoint tree figure"),
         fill = "#EEF2F6", cex = 0.72)
draw_box(7.35, 1.85, 3.25, 1.8,
         c("ROOTING CHECK", "Group E/H context pair lacks", "an exclusive supported split;", "no root applied"),
         fill = "#FFF7E8", border = "#B68A46", cex = 0.68)

# Reference proteins inform both BLASTP searches; candidate branches then combine.
draw_arrow(3.9, 7.45, 4.65, 7.35, dashed = TRUE, color = "#6A3D9A")
draw_arrow(3.55, 6.45, 4.65, 5.0, dashed = TRUE, color = "#6A3D9A")
draw_arrow(8.2, 7.15, 8.65, 6.55, color = "#0072B2")
draw_arrow(8.2, 4.35, 8.65, 5.65, color = "#CC79A7")
draw_arrow(11.1, 6.2, 11.55, 6.2)
draw_arrow(13.9, 6.2, 14.2, 6.2)
draw_arrow(14.95, 5.25, 14.2, 3.7)
draw_arrow(14.2, 5.35, 10.6, 3.05, dashed = TRUE, color = "#B68A46")
text(4.15, 7.7, "24 UGT84 BLASTP queries", cex = 0.62, col = "#5D536F", adj = 0)
text(4.12, 5.55, "reference set", cex = 0.62, col = "#5D536F", adj = 0)
text(8.3, 8.2, "plus reference panel", cex = 0.62, col = "#4D5962", adj = 0)
text(12.4, 4.55, "unrooted tree", cex = 0.62, col = "#4D5962", adj = 0)
dev.off()

cat("Created readme figures in ", outdir, "\n", sep = "")
cat("Reference panel: ", paste(names(ref_counts), ref_counts, sep = "=", collapse = "; "), "\n", sep = "")
cat("Cornus selected/tree/held out: ", paste(species, selected, in_tree, held, sep = ":", collapse = "; "), "\n", sep = "")
