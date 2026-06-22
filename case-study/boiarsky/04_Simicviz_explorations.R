setwd("/home/rshuwaikan/data_b/rshuwaikan/SIMIC/Boiarsky")
library(SimiCviz)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)

# Selected fit: lambda1 = 1e-06, lambda2 = 0.0001, opt R squared on eval 0.7533
plot_dir <- "./output/simicviz_out" # output dir

simic <- load_SimiCPipeline(project_dir = "../Boiarsky_run",
                               run_name = "NBM-SMM-MM",
                               lambda1 = "1e-06",
                               lambda2 = "0.0001")
                               
simic <- setLabelNames(simic, label_names  = c('NBM', 'SMM', 'MM'),
              colors = c("#3B7EA1", "#E66101", "#B2182B"))
simic

# > simic
# An object of class SimiCvizExperiment
#  3 label(s), 100 TF(s), 1000 target(s)
#  Weights: 3 matrices [0: 100 x 1000, 1: 100 x 1000, 2: 100 x 1000]
#  AUC: collected (28550 cells x 100 TFs)
#  Cell labels: 28550 cells across 3 label(s) [0, 1, 2]
#  Label names: 0 = NBM, 1 = SMM, 2 = MM
#  Colors: 0 = #3B7EA1, 1 = #E66101, 2 = #B2182B
#  TFs: MEF2D, FLI1, BACH1, CREM, ESRRA, KLF2, ...
#  Targets: IGKC, IGHA1, IGLC2, IGHG1, IGHG3, IGHG4, ...
#  Meta keys: adjusted_r_squared, run_name, lambda1, lambda2

# Check AUC scores: already calculated with SimiCPipeline
auc_scores <- simic@auc$collected
head(auc_scores[, 1:5])

# Assess quality of fitted models:
# load_SimiCPipeline()` function will already load the `adjusted_r_squared` values in the meta slot.
adjusted_r_squared <- simic@meta$adjusted_r_squared

plot_r2_distribution(adjusted_r_squared, simic, save =TRUE, out_dir = plot_dir, grid=c(1,3), w = 18)

# Select targets by adjusted R2
unselected_targets <- list()
selected_targets <- list()
lab_keys <- names(simic@label_names)
for (lab in lab_keys){
    # Save selected for plotting
    selected_targets[[lab]] <- simic@target_ids[which(adjusted_r_squared[[lab]] >= 0.7)]
    # Save unselected for reporting
    label <- simic@label_names[[lab]]
    unselected_targets[[label]] <- simic@target_ids[which(adjusted_r_squared[[lab]] < 0.7)]
}
print("Number of unselected targets per label:")
print(sapply(unselected_targets, length)) 

# NBM SMM  MM 
# 446  24  14 

# Start by looking at most dissimilar TFs: 
# Global dissimilarity scores: MinMax dissimilarity 
dis_score <- calculate_dissimilarity(simic)
top_tfs <- rownames(dis_score)
top_10 <- top_tfs[1:10]
# lets do top 10:
plot_dissimilarity_heatmap(simic, 
                           top_n = 10, 
                           cmap = "viridis",
                           save = TRUE,
                           out_dir = plot_dir,
                           filename = "dissimilarity_heatmap_top10TFs.pdf")


# Clustered mean AUC heatmap for the top 10: 
# Extracts the mean AUC per (TF, stage) from SimiCviz's plot_auc_heatmap()
# without saving its default un-clustered version (save = FALSE), then
# renders a clustered heatmap with ComplexHeatmap. Rows are clustered by
# Spearman correlation across stages so that TFs with similar activity
# trajectories sit together
# -----------------------------------------------------------------------------
auc_obj <- plot_auc_heatmap(
  simic,
  tf_names = top_10,
  save     = TRUE,
  cmap     = "plasma", 
  out_dir  = plot_dir,
  filename  = "AUC_mean_top10DissScore_unclustered.pdf"
)

mat <- auc_obj$data %>%
  dplyr::select(tf, condition, auc) %>%
  dplyr::mutate(tf = as.character(tf),
                condition = as.character(condition)) %>%
  tidyr::pivot_wider(names_from = condition, values_from = auc) %>%
  tibble::column_to_rownames("tf") %>%
  as.matrix()

col_fun <- colorRamp2(
  seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 100),
  viridis(100, option = "plasma")
)

ht <- Heatmap(
  mat,
  name                     = "Mean AUC",
  col                      = col_fun,
  cluster_rows             = TRUE,
  cluster_columns          = FALSE,
  clustering_distance_rows = "spearman",
  clustering_method_rows   = "complete",
  column_order             = c("NBM", "SMM", "MM"),
  row_names_side           = "left",
  column_names_side        = "bottom",
  column_names_rot         = 0,
  row_names_centered       = FALSE,
  row_title                = "TF",
  column_title             = "Phenotype",
  column_title_side        = "bottom",
  row_names_gp             = gpar(fontsize = 9),
  column_names_gp          = gpar(fontsize = 9),
  row_title_gp             = gpar(fontsize = 9),
  column_title_gp          = gpar(fontsize = 9),
  heatmap_legend_param     = list(title    = "Mean AUC",
                                  title_gp = gpar(fontsize = 9, fontface = "plain")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 9))
  }
)

pdf(file.path(plot_dir, "AUC_mean_top10DissScore_clustered.pdf"), height = 5)
draw(ht,
     column_title    = "Mean Activity Score per TF: clustered",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

# Explore: 
# Plot AUC distributions for the top 10 most dissimilar TFs:
plot_auc_distributions(
  simic,
  tf_names  = top_10,
  fill      = TRUE,
  alpha     = 0.6,
  bw_adjust = 1 / 8,
  rug       = TRUE,
  save      = TRUE,
  out_dir   = plot_dir,
  filename  = "AUC_distributions_top10DissScore.pdf"
)

# Plot cumulative AUC distributions for the top 10 most dissimilar TFs:
plot_auc_cumulative(
  simic,
  tf_names      = top_10,
  rug           = TRUE,
  include_table = TRUE,
  save          = TRUE,
  out_dir       = plot_dir,
  filename      = "AUC_cumulative_top10DissScore.pdf"
)


# -----------------------------------------------------------------------------
# Per-TF case-study panels
#
# For each of the selected TFs, we produce three PDFs that will be
# composed into a single supplementary figure per TF:
#   <TF>_density.pdf                       (panel A: per-cell activity density)
#   <TF>_cumulative.pdf                    (panel B: ECDF of activity)
#   <TF>_network_heatmap_top50targets.pdf  (panel C: top-50 target heatmap)
#

# -----------------------------------------------------------------------------
CASE_STUDY_TFS <- c("MAF", "IRF2")
PLOT_DIR <- plot_dir
R2_THRESHOLD <- 0.7

for (tf in CASE_STUDY_TFS) {

  # Panel A: per-cell activity-score density
  plot_auc_distributions(
    simic,
    tf_names  = tf,
    fill      = TRUE,
    alpha     = 0.6,
    bw_adjust = 1 / 8,
    rug       = TRUE,
    save      = TRUE,
    out_dir   = PLOT_DIR,
    filename  = paste0(tf, "_density.pdf")
  )

  # Panel B: cumulative activity-score distribution (ECDF)
  plot_auc_cumulative(
    simic,
    tf_names      = tf,
    rug           = TRUE,
    include_table = TRUE,
    save          = TRUE,
    out_dir       = PLOT_DIR,
    filename      = paste0(tf, "_cumulative.pdf")
  )

  # Panel C: regulatory-network heatmap of the top 50 targets
  plot_tf_network_heatmap(
    simic,
    tf_name      = tf,
    top_n        = 50,
    r2_threshold = R2_THRESHOLD,
    show_values  = TRUE,
    cmap         = c("purple", "white", "yellow"),
    save         = TRUE,
    out_dir      = PLOT_DIR,
    filename     = paste0(tf, "_network_heatmap_top50targets.pdf")
  )
}


# -----------------------------------------------------------------------------
# 7. Exploratory target heatmaps (not in supplementary)
#
# E2F4 clusters with the progressive group but its top targets in MM
# include hemoglobin chains (HBB, HBA1, HBA2). This is most likely an
# erythroid contamination / ambient-RNA artefact rather than E2F4 biology,
# so the heatmap is generated for inspection but not shown in the
# supplementary.
# -----------------------------------------------------------------------------
for (tf in EXPLORATORY_TFS) {
  plot_tf_network_heatmap(
    simic,
    tf_name      = tf,
    top_n        = 50,
    r2_threshold = R2_THRESHOLD,
    show_values  = TRUE,
    cmap         = c("purple", "white", "yellow"),
    save         = TRUE,
    out_dir      = PLOT_DIR,
    filename     = paste0(tf, "_network_heatmap_top50targets.pdf")
  )
}


# Edit IRF2 target network, making sure columns follow correct order of phenotypes (NBM, SMM, MM): 
irf2_net <- get_tf_network(simic, tf_name = "IRF2",
                           labels = c(0L, 1L, 2L),
                           r2_threshold = R2_THRESHOLD)
irf2_abs <- apply(abs(irf2_net), 1, max, na.rm = TRUE)
irf2_abs[!is.finite(irf2_abs)] <- 0
irf2_top <- names(sort(irf2_abs, decreasing = TRUE))[seq_len(min(50, nrow(irf2_net)))]
irf2_mat  <- as.matrix(irf2_net[irf2_top, c("NBM", "SMM", "MM")])

irf2_lim    <- max(abs(irf2_mat), na.rm = TRUE)
irf2_colfun <- colorRamp2(c(-irf2_lim, 0, irf2_lim),
                           c("purple", "white", "yellow"))

irf2_ht <- Heatmap(
  irf2_mat,
  name             = "Weight",
  col              = irf2_colfun,
  na_col           = "darkgrey",
  cluster_rows     = FALSE,
  cluster_columns  = FALSE,
  column_order     = c("NBM", "SMM", "MM"),
  row_names_side   = "left",
  column_names_side = "bottom",
  column_names_rot = 45,
  row_names_gp     = gpar(fontsize = 9),
  column_names_gp  = gpar(fontsize = 10),
  column_title_gp  = gpar(fontsize = 9),
  heatmap_legend_param = list(title    = "Weight",
                              title_gp = gpar(fontsize = 9, fontface = "plain")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    val <- irf2_mat[i, j]
    if (is.na(val)) {
      grid.text("< R² thr.", x, y, gp = gpar(fontsize = 6, col = "black"))
    } else {
      grid.text(sprintf("%.2f", val), x, y, gp = gpar(fontsize = 7, col = "black"))
    }
  }
)

pdf(file.path(PLOT_DIR, "IRF2_network_heatmap_top50targets.pdf"),
    width = 8, height = max(6, 50 * 0.4 + 2))
draw(irf2_ht,
     column_title    = " IRF2 Regulatory Network\nTop 50 targets across phenotypes",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()
