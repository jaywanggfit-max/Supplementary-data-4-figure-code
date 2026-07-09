suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

figure_id <- "Supplementary_Fig_24"
figure_title <- "Flagellar structural Orthogroups (55 MPa Droplet vs 55 MPa Bulk)"
output_svg_name <- "Supplementary_Fig_24.svg"
output_png_name <- "Supplementary_Fig_24.png"
expected_contrast <- "D55_vs_B55"
expected_contrast_label <- "55D vs 55B"
width_mm <- 196.2
height_mm <- 112.9
width_px <- 2067
height_px <- 1189

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else ""
script_dir <- if (nzchar(script_file)) {
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
source_file <- if (length(args) >= 1) args[[1]] else file.path(script_dir, "source_data.tsv")
out_dir <- if (length(args) >= 2) args[[2]] else script_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strain_order <- c("Halomonas", "Psychrobacter", "Sphingopyxis", "Aequorivita")
base_family <- "Helvetica"
axis_col <- "#2F2F2F"
grid_col <- "#E8ECEA"
title_fill <- "#F7FAF8"
title_border <- "#CAD8D2"
cross_col <- "#A7AFB4"
fc_low <- "#2F6DA3"
fc_mid <- "#F2F0ED"
fc_high <- "#C34E61"
line_pt <- 0.5 * 0.352778

read_tsv <- function(path) {
  read.delim(
    path,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

num <- function(x) suppressWarnings(as.numeric(x))

sanitize_svg_text_spacing <- function(svg_file) {
  txt <- readLines(svg_file, warn = FALSE, encoding = "UTF-8")
  txt <- gsub("\\s+textLength='[^']*'\\s+lengthAdjust='spacingAndGlyphs'", "", txt)
  writeLines(txt, svg_file, useBytes = TRUE)
}

make_dashed_ellipse <- function(points, x_col = "x", y_col = "y",
                                rx = 0.20, ry = 0.13, n = 144) {
  if (nrow(points) == 0) {
    return(tibble(x = numeric(), y = numeric(), group_id = character()))
  }
  theta <- seq(0, 2 * pi, length.out = n)
  bind_rows(lapply(seq_len(nrow(points)), function(i) {
    tibble(
      x = points[[x_col]][i] + rx * cos(theta),
      y = points[[y_col]][i] + ry * sin(theta),
      group_id = paste0("ns_", i)
    )
  }))
}

make_title_panel <- function(title) {
  ggplot() +
    annotate("rect", xmin = -0.20, xmax = 0.98, ymin = 0.10, ymax = 0.92,
             fill = title_fill, colour = title_border, linewidth = line_pt) +
    annotate("text", x = -0.12, y = 0.51, label = title, family = base_family,
             hjust = 0, size = 3.85, colour = axis_col) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void(base_family = base_family) +
    theme(plot.margin = margin(0, 0, 0, 0))
}

make_legend_panel <- function(fc_limit, size_limit) {
  grad_fun <- scales::col_numeric(c(fc_low, fc_mid, fc_high), domain = c(-fc_limit, 0, fc_limit))
  grad_raster <- matrix(grad_fun(seq(-fc_limit, fc_limit, length.out = 768)), nrow = 1)
  size_values <- c(2, 6, 10)
  if (size_limit > 10) size_values <- c(size_values, size_limit)
  size_values <- unique(size_values[size_values <= size_limit])
  size_df <- tibble(
    x = seq(0.20, 0.80, length.out = length(size_values)),
    y = 0.51,
    value = size_values
  )
  ns_ring <- make_dashed_ellipse(tibble(x = 0.21, y = 0.13), rx = 0.038, ry = 0.052)

  ggplot() +
    annotate("text", x = 0.10, y = 0.98, label = expression(log[2] * "(fold change)"),
             hjust = 0, size = 2.65, family = base_family, colour = axis_col) +
    annotation_raster(
      raster = grad_raster,
      xmin = 0.10, xmax = 0.92, ymin = 0.82, ymax = 0.91,
      interpolate = TRUE
    ) +
    geom_rect(aes(xmin = 0.10, xmax = 0.92, ymin = 0.82, ymax = 0.91),
              fill = NA, colour = axis_col, linewidth = line_pt) +
    annotate("text", x = c(0.10, 0.51, 0.92), y = 0.74,
             label = c(paste0("-", fc_limit), "0", paste0("+", fc_limit)),
             size = 2.15, family = base_family, colour = axis_col) +
    annotate("text", x = 0.10, y = 0.63, label = expression(-log[10] * "(P value)"),
             hjust = 0, size = 2.65, family = base_family, colour = axis_col) +
    geom_point(data = size_df, aes(x = x, y = y, size = value),
               shape = 21, fill = "#D7D7D7", colour = "#777777", stroke = 0.25) +
    geom_text(data = size_df, aes(x = x, y = 0.39, label = value),
              size = 2.05, family = base_family, colour = axis_col) +
    geom_path(data = ns_ring, aes(x = x, y = y, group = group_id),
              colour = cross_col, linewidth = 0.24, linetype = "22", lineend = "round") +
    annotate("text", x = 0.31, y = 0.13, label = "n.s.",
             hjust = 0, vjust = 0.5, size = 2.05, family = base_family, colour = "#646C70") +
    geom_point(aes(x = 0.58, y = 0.13), shape = 4, size = 2.55,
               stroke = 0.60, colour = "#747C80") +
    annotate("text", x = 0.68, y = 0.13, label = "Absent",
             hjust = 0, vjust = 0.5, size = 2.05, family = base_family, colour = axis_col) +
    scale_size_area(max_size = 5.1, limits = c(0, size_limit), guide = "none") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void(base_family = base_family) +
    theme(plot.margin = margin(0, 0, 0, 3))
}

validate_source <- function(dat) {
  required_cols <- c(
    "contrast", "contrast_label", "orthogroup_id", "gene_label", "display_gene_label",
    "x", "strain", "point_status", "member_gene_count_in_strain",
    "log2_fold_change", "P_value", "FDR", "neg_log10_P_value", "direction", "sharing"
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  if (!all(dat$contrast == expected_contrast)) stop("Unexpected contrast values in source_data.tsv.")
  if (!all(dat$contrast_label == expected_contrast_label)) stop("Unexpected contrast_label values in source_data.tsv.")
  valid_status <- c("DEG", "present_not_significant", "absent_no_member_gene")
  bad_status <- setdiff(unique(dat$point_status), valid_status)
  if (length(bad_status) > 0) stop("Unexpected point_status values: ", paste(bad_status, collapse = ", "))
  dup <- dat %>% count(orthogroup_id, strain) %>% filter(n != 1)
  if (nrow(dup) > 0) stop("Non-unique orthogroup-strain cells detected.")
  deg_missing <- dat %>%
    filter(point_status == "DEG") %>%
    filter(is.na(num(log2_fold_change)) | is.na(num(P_value)) | is.na(num(neg_log10_P_value)))
  if (nrow(deg_missing) > 0) stop("DEG rows with missing log2_fold_change, P_value, or neg_log10_P_value.")
  invisible(TRUE)
}

make_matrix_panel <- function(dat, fc_limit, size_limit) {
  og_order <- dat %>%
    distinct(x, orthogroup_id, display_gene_label) %>%
    mutate(x = num(x)) %>%
    arrange(x)
  n_ogs <- nrow(og_order)

  plot_dat <- dat %>%
    mutate(
      x = num(x),
      strain = factor(strain, levels = strain_order),
      strain_y = case_when(
        strain == "Halomonas" ~ 4,
        strain == "Psychrobacter" ~ 3,
        strain == "Sphingopyxis" ~ 2,
        strain == "Aequorivita" ~ 1,
        TRUE ~ NA_real_
      ),
      logFC = num(log2_fold_change),
      neg_log10_p = num(neg_log10_P_value)
    )

  ns_rings <- make_dashed_ellipse(
    filter(plot_dat, point_status == "present_not_significant"),
    x_col = "x",
    y_col = "strain_y",
    rx = 0.20,
    ry = 0.13
  )

  ggplot(plot_dat, aes(x = x, y = strain_y)) +
    geom_vline(xintercept = seq_len(n_ogs), colour = grid_col, linewidth = line_pt) +
    geom_hline(yintercept = 1:4, colour = grid_col, linewidth = line_pt) +
    geom_path(data = ns_rings, aes(x = x, y = y, group = group_id),
              inherit.aes = FALSE, colour = cross_col, linewidth = 0.25,
              linetype = "22", lineend = "round") +
    geom_point(data = filter(plot_dat, point_status == "absent_no_member_gene"),
               shape = 4, size = 2.8, stroke = 0.62, colour = "#747C80") +
    geom_point(data = filter(plot_dat, point_status == "DEG"),
               aes(fill = logFC, size = neg_log10_p),
               shape = 21, colour = "#707A80", stroke = 0.22, alpha = 0.96) +
    annotate("rect", xmin = 0.45, xmax = n_ogs + 0.55, ymin = 0.50, ymax = 4.50,
             fill = NA, colour = axis_col, linewidth = line_pt) +
    scale_x_continuous(
      breaks = og_order$x,
      labels = og_order$display_gene_label,
      expand = expansion(mult = c(0.015, 0.015))
    ) +
    scale_y_continuous(
      breaks = c(4, 3, 2, 1),
      labels = strain_order,
      limits = c(0.45, 4.55),
      expand = c(0, 0)
    ) +
    scale_fill_gradient2(
      low = fc_low, mid = fc_mid, high = fc_high, midpoint = 0,
      limits = c(-fc_limit, fc_limit), oob = squish, guide = "none"
    ) +
    scale_size_area(max_size = 6.2, limits = c(0, size_limit), guide = "none") +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = base_family, base_size = 8) +
    theme(
      text = element_text(family = base_family, colour = axis_col),
      axis.title = element_blank(),
      axis.text.x = element_text(angle = 65, hjust = 1, vjust = 1, size = 9.0, face = "italic"),
      axis.text.y = element_text(size = 10.0, face = "italic", margin = margin(r = 5)),
      panel.grid = element_blank(),
      plot.margin = margin(0, 2, 2, 2)
    )
}

dat <- read_tsv(source_file)
validate_source(dat)

fc_limit <- 5
max_neglog <- max(num(dat$neg_log10_P_value), na.rm = TRUE)
size_limit <- max(10, ceiling(max_neglog / 2) * 2)

top_row <- make_title_panel(figure_title) | plot_spacer()
top_row <- top_row + plot_layout(widths = c(0.76, 0.24))
bottom_row <- make_matrix_panel(dat, fc_limit, size_limit) | make_legend_panel(fc_limit, size_limit)
bottom_row <- bottom_row + plot_layout(widths = c(0.72, 0.28))

fig <- (top_row / bottom_row) +
  plot_layout(heights = c(0.12, 0.88))

svg_file <- file.path(out_dir, output_svg_name)
svglite::svglite(svg_file, width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
print(fig)
dev.off()
sanitize_svg_text_spacing(svg_file)

png_file <- file.path(out_dir, output_png_name)
ragg::agg_png(png_file, width = width_px, height = height_px, units = "px", res = 300, background = "white")
print(fig)
dev.off()

qa <- tibble(
  figure = figure_id,
  output_svg = output_svg_name,
  output_png = output_png_name,
  source_data = basename(source_file),
  contrast = expected_contrast,
  contrast_label = expected_contrast_label,
  structural_orthogroups = n_distinct(dat$orthogroup_id),
  total_matrix_cells = nrow(dat),
  DEG_cells = sum(dat$point_status == "DEG"),
  present_not_significant_cells = sum(dat$point_status == "present_not_significant"),
  absent_no_member_gene_cells = sum(dat$point_status == "absent_no_member_gene"),
  max_neg_log10_P_value = max_neglog,
  log2_fold_change_colour_limit = fc_limit,
  width_mm = width_mm,
  height_mm = height_mm
)

write.table(qa, file.path(out_dir, "plot_QA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
print(qa)
