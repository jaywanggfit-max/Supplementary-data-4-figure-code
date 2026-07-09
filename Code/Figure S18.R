required_packages <- c("ComplexHeatmap", "circlize", "grid", "grDevices")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with: BiocManager::install(c('ComplexHeatmap', 'circlize'))"
  )
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(grDevices)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("plot_heatmap_both.R", winslash = "/", mustWork = FALSE)
}
base_dir <- dirname(script_path)
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

matrix_file <- file.path(data_dir, "heatmap_matrix.csv")
sample_file <- file.path(data_dir, "sample_annotation.csv")
metadata_file <- file.path(data_dir, "row_metadata.csv")

pdf_file <- file.path(out_dir, "heatmap_both.pdf")
png_file <- file.path(out_dir, "heatmap_both_preview.png")
svg_file <- file.path(out_dir, "heatmap_both.svg")
qa_file <- file.path(out_dir, "heatmap_both_QA.csv")

strain_order <- c("Aequorivita", "Psychrobacter", "Sphingopyxis", "Halomonas")
condition_order <- c("D55", "D0", "B55", "B0")
condition_labels <- c(
  D55 = "DeepDrop 55 MPa",
  D0 = "Droplet 0.1 MPa",
  B55 = "Bulk 55 MPa",
  B0 = "Bulk 0.1 MPa"
)
condition_cols <- c(
  D55 = "#D19AA0",
  D0 = "#A986C2",
  B55 = "#C5D4EA",
  B0 = "#BDBDBD"
)

heat_cols <- colorRamp2(c(-4, 0, 4), c("#3B5B9B", "#F7F4F1", "#D84A35"))

matrix_df <- read_csv(matrix_file)
if (!"orthofinder_orthogroup" %in% names(matrix_df)) {
  stop("heatmap_matrix.csv must contain an orthofinder_orthogroup column.")
}
rownames(matrix_df) <- make.unique(matrix_df$orthofinder_orthogroup)
mat_raw <- as.matrix(matrix_df[, setdiff(names(matrix_df), "orthofinder_orthogroup"), drop = FALSE])
storage.mode(mat_raw) <- "numeric"

sample_ann <- read_csv(sample_file)
required_sample_cols <- c("sample", "strain", "condition", "heatmap_column_order")
missing_sample_cols <- setdiff(required_sample_cols, names(sample_ann))
if (length(missing_sample_cols) > 0) {
  stop("sample_annotation.csv is missing: ", paste(missing_sample_cols, collapse = ", "))
}
sample_ann <- sample_ann[match(colnames(mat_raw), sample_ann$sample), , drop = FALSE]
if (any(is.na(sample_ann$sample))) {
  stop("Some heatmap matrix columns are missing from sample_annotation.csv.")
}
sample_ann$strain <- factor(sample_ann$strain, levels = strain_order)
sample_ann$condition <- factor(sample_ann$condition, levels = condition_order)
if (any(is.na(sample_ann$strain)) || any(is.na(sample_ann$condition))) {
  stop("Unexpected strain or condition level in sample_annotation.csv.")
}

row_meta <- read_csv(metadata_file)
if (nrow(row_meta) != nrow(mat_raw)) {
  stop("row_metadata.csv row count does not match heatmap_matrix.csv.")
}

# Target layout: samples are rows, orthogroups are columns.
mat <- t(mat_raw)
rownames(mat) <- sample_ann$sample

row_split <- factor(sample_ann$strain, levels = strain_order)
row_condition <- factor(sample_ann$condition, levels = condition_order)

left_ha <- rowAnnotation(
  Condition = row_condition,
  col = list(Condition = condition_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  width = unit(3.4, "mm"),
  gp = gpar(col = "#595959", lwd = 0.35)
)

draw_condition_legend <- function() {
  x0 <- unit(0.014, "npc")
  y0 <- unit(0.048, "npc")
  grid.text(
    "Condition",
    x = x0, y = unit(0.085, "npc"),
    just = c("left", "bottom"),
    gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica")
  )
  legend_items <- c("D55", "B55", "D0", "B0")
  item_x <- c(0.025, 0.170, 0.025, 0.170)
  item_y <- c(0.058, 0.058, 0.026, 0.026)
  for (i in seq_along(legend_items)) {
    key <- legend_items[i]
    grid.rect(
      x = unit(item_x[i], "npc"), y = unit(item_y[i], "npc"),
      width = unit(0.014, "npc"), height = unit(0.016, "npc"),
      just = c("left", "center"),
      gp = gpar(fill = condition_cols[[key]], col = "#8A6B75", lwd = 0.55)
    )
    grid.text(
      condition_labels[[key]],
      x = unit(item_x[i] + 0.018, "npc"), y = unit(item_y[i], "npc"),
      just = c("left", "center"),
      gp = gpar(fontsize = 10, fontfamily = "Helvetica")
    )
  }
}

draw_z_legend <- function() {
  x0 <- 0.812
  x1 <- 0.988
  y0 <- 0.036
  y1 <- 0.054
  grid.text(
    "Row Z-score",
    x = unit(x0, "npc"), y = unit(0.082, "npc"),
    just = c("left", "bottom"),
    gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica")
  )
  vals <- seq(-4, 4, length.out = 512)
  col_fun <- colorRampPalette(c("#3B5B9B", "#F7F4F1", "#D84A35"))
  grad <- matrix(col_fun(length(vals)), nrow = 1)
  grid.raster(
    as.raster(grad),
    x = unit(x0, "npc"), y = unit(y0, "npc"),
    width = unit(x1 - x0, "npc"), height = unit(y1 - y0, "npc"),
    just = c("left", "bottom"),
    interpolate = TRUE
  )
  grid.rect(
    x = unit(x0, "npc"), y = unit(y0, "npc"),
    width = unit(x1 - x0, "npc"), height = unit(y1 - y0, "npc"),
    just = c("left", "bottom"),
    gp = gpar(fill = NA, col = "#333333", lwd = 0.45)
  )
  ticks <- c(-4, -2, 0, 2, 4)
  tx <- x0 + (ticks + 4) / 8 * (x1 - x0)
  for (i in seq_along(ticks)) {
    grid.lines(
      x = unit(rep(tx[i], 2), "npc"),
      y = unit(c(y0 + 0.001, y1 - 0.001), "npc"),
      gp = gpar(col = "white", lwd = 0.7)
    )
    grid.text(
      as.character(ticks[i]),
      x = unit(tx[i], "npc"), y = unit(y0 - 0.010, "npc"),
      just = c("center", "top"),
      gp = gpar(fontsize = 10, fontfamily = "Helvetica", col = "#555555")
    )
  }
}

draw_full_heatmap <- function() {
  set.seed(1)
  ht <- Heatmap(
    mat,
    name = "Row Z-score",
    col = heat_cols,
    left_annotation = left_ha,
    row_split = row_split,
    cluster_row_slices = FALSE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    row_dend_width = unit(20, "mm"),
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10, fontface = "bold.italic", fontfamily = "Helvetica"),
    column_km = 4,
    cluster_column_slices = FALSE,
    cluster_columns = TRUE,
    show_column_dend = FALSE,
    column_title = NULL,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE,
    border = TRUE,
    border_gp = gpar(col = "black", lwd = 1.0),
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(at = c(-4, -2, 0, 2, 4)),
    use_raster = TRUE,
    raster_quality = 2
  )

  draw(
    ht,
    annotation_legend_list = list(),
    show_annotation_legend = FALSE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    padding = unit(c(42, 4, 4, 4), "mm")
  )
  upViewport(0)
  draw_condition_legend()
  draw_z_legend()
}

export_all <- function() {
  pdf(pdf_file, width = 13.8, height = 10.8, useDingbats = FALSE, family = "Helvetica")
  draw_full_heatmap()
  dev.off()

  png(png_file, width = 2760, height = 2160, res = 200, type = "cairo", bg = "white")
  draw_full_heatmap()
  dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(svg_file, width = 13.8, height = 10.8, bg = "white")
    draw_full_heatmap()
    dev.off()
  }
}

export_all()

qa <- data.frame(
  check = c(
    "matrix_rows_orthogroups",
    "matrix_columns_samples",
    "heatmap_rows_after_transpose_samples",
    "heatmap_columns_after_transpose_orthogroups",
    "sample_annotation_matched",
    "row_metadata_matched",
    "row_split",
    "column_split",
    "condition_annotation_only",
    "strain_colour_annotation"
  ),
  value = c(
    nrow(mat_raw),
    ncol(mat_raw),
    nrow(mat),
    ncol(mat),
    all(rownames(mat) == sample_ann$sample),
    nrow(row_meta) == ncol(mat),
    paste(levels(row_split), collapse = ";"),
    "column_km = 4",
    "TRUE",
    "FALSE"
  )
)
write.csv(qa, qa_file, row.names = FALSE, quote = TRUE)

message("Exported: ", file.path("output", basename(pdf_file)))
message("Exported: ", file.path("output", basename(png_file)))
if (file.exists(svg_file)) {
  message("Exported: ", file.path("output", basename(svg_file)))
}
