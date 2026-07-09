required_packages <- c("edgeR", "ggplot2", "dplyr", "tidyr", "patchwork", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("plot_Supplementary_Fig17.R")
base_dir <- dirname(script_path)
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_csv <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

sample_info <- read_csv(file.path(data_dir, "sample_info.csv")) %>%
  mutate(condition_label = factor(condition_label, levels = c("55B", "0B", "55D", "0D")))
counts_long <- read_csv(file.path(data_dir, "counts_long_for_pca.csv"))

strains <- c("Halomonas", "Sphingopyxis", "Aequorivita", "Psychrobacter")
condition_labels <- c(
  "55B" = "Bulk 55 MPa",
  "0B" = "Bulk 0.1 MPa",
  "55D" = "DeepDrop 55 MPa",
  "0D" = "Droplet 0.1 MPa"
)
condition_cols <- c("55B" = "#C9D7EA", "0B" = "#D8D8D8", "55D" = "#D79A9A", "0D" = "#A78BB7")
axis_flip <- data.frame(
  strain = strains,
  flip_PC1 = c(1, -1, 1, -1),
  flip_PC2 = c(-1, 1, 1, -1),
  stringsAsFactors = FALSE
)

normalise_dge <- function(y) {
  if ("normLibSizes" %in% getNamespaceExports("edgeR")) {
    edgeR::normLibSizes(y)
  } else {
    edgeR::calcNormFactors(y)
  }
}

compute_one_strain <- function(st) {
  dd <- counts_long %>% filter(strain == st)
  mat <- dd %>%
    select(orthofinder_orthogroup, sample, raw_count) %>%
    pivot_wider(names_from = sample, values_from = raw_count, values_fill = 0)
  orthogroups <- mat$orthofinder_orthogroup
  mat <- as.data.frame(mat[, -1, drop = FALSE])
  rownames(mat) <- orthogroups

  si <- sample_info %>% filter(strain == st, sample %in% colnames(mat))
  mat <- round(as.matrix(mat[, si$sample, drop = FALSE]))

  y <- DGEList(mat)
  keep <- filterByExpr(y, group = factor(si$condition))
  y <- normalise_dge(y[keep, , keep.lib.sizes = FALSE])
  expr <- log2(cpm(y, log = FALSE) + 1)

  pca <- prcomp(t(expr), center = TRUE, scale. = TRUE)
  ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(
    strain = st,
    sample = rownames(pca$x),
    condition = si$condition[match(rownames(pca$x), si$sample)],
    condition_label = si$condition_label[match(rownames(pca$x), si$sample)],
    PC1_for_plot = pca$x[, 1],
    PC2_for_plot = pca$x[, 2],
    PC1_percent = ve[1],
    PC2_percent = ve[2],
    n_orthogroups_used = nrow(expr),
    stringsAsFactors = FALSE
  ) %>%
    left_join(axis_flip, by = "strain") %>%
    mutate(
      PC1_for_plot = PC1_for_plot * flip_PC1,
      PC2_for_plot = PC2_for_plot * flip_PC2
    ) %>%
    select(-flip_PC1, -flip_PC2)

  axis <- data.frame(
    strain = st,
    n_samples = ncol(expr),
    n_orthogroups_used = nrow(expr),
    PC1_percent = ve[1],
    PC2_percent = ve[2],
    normalization = "filterByExpr, TMM normalization, log2(CPM + 1)",
    centering = "features centered and scaled for PCA",
    stringsAsFactors = FALSE
  )

  list(scores = scores, axis = axis)
}

pca_result <- lapply(strains, compute_one_strain)
scores <- bind_rows(lapply(pca_result, `[[`, "scores")) %>%
  mutate(
    strain = factor(strain, levels = strains),
    condition_label = factor(condition_label, levels = c("55B", "0B", "55D", "0D"))
  )
axis_percent <- bind_rows(lapply(pca_result, `[[`, "axis"))

write.csv(scores, file.path(data_dir, "pca_scores.csv"), row.names = FALSE)
write.csv(axis_percent, file.path(data_dir, "pca_axis_percent.csv"), row.names = FALSE)

theme_set(
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.title = element_blank(),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 7, face = "italic", hjust = 0.5),
      panel.border = element_rect(fill = NA, linewidth = 0.3, colour = "black")
    )
)

make_panel <- function(st) {
  dd <- scores %>% filter(strain == st)
  ax <- axis_percent %>% filter(strain == st)
  ggplot(dd, aes(PC1_for_plot, PC2_for_plot, colour = condition_label)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey86", linetype = "dotted") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey86", linetype = "dotted") +
    geom_point(size = 1.8, alpha = 0.82) +
    scale_colour_manual(values = condition_cols, labels = condition_labels, drop = FALSE) +
    labs(
      title = st,
      x = paste0("PC1 (", sprintf("%.1f", ax$PC1_percent), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", ax$PC2_percent), "%)")
    )
}

p <- (make_panel("Halomonas") + make_panel("Sphingopyxis")) /
  (make_panel("Aequorivita") + make_panel("Psychrobacter")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top")

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 450) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svg_file <- file.path(out_dir, paste0(stem, ".svg"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))
  svglite::svglite(svg_file, width = w, height = h)
  print(plot)
  dev.off()
  ragg::agg_png(png_file, width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  message("Exported: ", file.path("output", basename(svg_file)))
  message("Exported: ", file.path("output", basename(png_file)))
}

save_figure(p, "Supplementary_Fig17_PCA", width_mm = 145, height_mm = 135)
