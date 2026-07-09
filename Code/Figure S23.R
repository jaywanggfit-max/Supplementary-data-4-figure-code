suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

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

figure_id <- "Figure_j_55D_vs_0p1D"
output_svg <- file.path(out_dir, paste0(figure_id, ".svg"))
output_png <- file.path(out_dir, paste0(figure_id, ".png"))
width_mm <- 206.7
height_mm <- 53.9
width_px <- 2067
height_px <- 539

base_family <- "Helvetica"
strain_order <- c("Halomonas", "Psychrobacter", "Sphingopyxis", "Aequorivita")
module_order <- c(
  "membrane maintenance",
  "ROS removal",
  "Flagellar assembly",
  "proteostasis",
  "osmoregulation"
)
module_labels <- c(
  "membrane maintenance" = "Membrane\nmaintenance",
  "ROS removal" = "ROS removal",
  "Flagellar assembly" = "Flagellar",
  "proteostasis" = "Proteostasis",
  "osmoregulation" = "Osmoregulation"
)
module_colours <- c(
  "membrane maintenance" = "#B56F69",
  "ROS removal" = "#D68A3E",
  "Flagellar assembly" = "#8B4B96",
  "proteostasis" = "#4D9A4A",
  "osmoregulation" = "#3D6EA5"
)

fc_low <- "#3F6692"
fc_mid <- "#F2F0ED"
fc_high <- "#B84E58"
axis_col <- "#3C3C3C"
not_deg_col <- "#9B9B9B"
absent_col <- "#737A80"
line_pt <- 0.5 * 0.352778
not_deg_linetype <- "22"

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

format_fc_tick <- function(x) {
  ifelse(x > 0, paste0("+", format(x, trim = TRUE, scientific = FALSE)),
         format(x, trim = TRUE, scientific = FALSE))
}

sanitize_svg_text_spacing <- function(svg_file) {
  txt <- readLines(svg_file, warn = FALSE, encoding = "UTF-8")
  txt <- gsub("\\s+textLength='[^']*'\\s+lengthAdjust='spacingAndGlyphs'", "", txt)
  writeLines(txt, svg_file, useBytes = TRUE)
}

make_dashed_ellipse <- function(points, x_col = "x", y_col = "y",
                                rx = 0.22, ry = 0.13, n = 144) {
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

validate_source <- function(dat) {
  required_cols <- c(
    "contrast", "module", "orthogroup_id", "gene_label", "display_gene_label",
    "x", "strain", "point_class", "status_raw", "member_gene_count_in_strain",
    "logFC", "pvalue", "FDR", "neg_log10_p", "n_DEG_events",
    "n_DEG_strains", "DEG_strains", "directions", "strict_function_evidence"
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  if (!all(dat$contrast == "D55_vs_D0")) stop("source_data.tsv must contain only D55_vs_D0 rows.")
  bad_modules <- setdiff(unique(dat$module), module_order)
  if (length(bad_modules) > 0) stop("Unexpected module(s): ", paste(bad_modules, collapse = ", "))
  bad_status <- setdiff(unique(dat$point_class), c("DEG", "Present_not_DEG", "Absent"))
  if (length(bad_status) > 0) stop("Unexpected point_class value(s): ", paste(bad_status, collapse = ", "))
  dup <- dat %>%
    count(contrast, module, orthogroup_id, strain) %>%
    filter(n != 1)
  if (nrow(dup) > 0) stop("Non-unique contrast-module-OG-strain cells detected.")
  deg_missing <- dat %>%
    filter(point_class == "DEG") %>%
    filter(is.na(num(logFC)) | is.na(num(pvalue)) | is.na(num(FDR)) | is.na(num(neg_log10_p)))
  if (nrow(deg_missing) > 0) stop("DEG rows with missing logFC/P value/FDR detected.")
  invisible(TRUE)
}

make_legend_panel <- function(fc_limit, size_limit) {
  grad_fun <- scales::col_numeric(c(fc_low, fc_mid, fc_high), domain = c(-fc_limit, 0, fc_limit))
  grad_raster <- matrix(grad_fun(seq(-fc_limit, fc_limit, length.out = 768)), nrow = 1)
  size_values <- c(2, 5, 10, 15, 20)
  if (size_limit > 20) size_values <- c(size_values, size_limit)
  size_values <- unique(size_values[size_values <= size_limit])
  size_df <- tibble(
      x = seq(7.55, 9.75, length.out = length(size_values)),
    y = 1.94,
    value = size_values
  )
  ns_ring <- make_dashed_ellipse(tibble(x = 10.28, y = 1.94), rx = 0.14, ry = 0.18)

  ggplot() +
    annotate("text", x = 4.00, y = 1.95, label = expression(Log[2] * "(Fold change)"),
             hjust = 0, size = 2.55, family = base_family, colour = axis_col) +
    annotation_raster(
      raster = grad_raster,
      xmin = 4.00,
      xmax = 5.85,
      ymin = 1.78,
      ymax = 2.12,
      interpolate = TRUE
    ) +
    geom_rect(
      aes(xmin = 4.00, xmax = 5.85, ymin = 1.78, ymax = 2.12),
      fill = NA, colour = axis_col, linewidth = line_pt
    ) +
    annotate("text", x = seq(4.00, 5.85, length.out = 5), y = 1.28,
             label = format_fc_tick(c(-fc_limit, -fc_limit / 2, 0, fc_limit / 2, fc_limit)),
             size = 2.15, family = base_family, colour = axis_col) +
    geom_segment(
      aes(x = seq(4.00, 5.85, length.out = 5), xend = seq(4.00, 5.85, length.out = 5),
          y = 1.78, yend = 1.64),
      inherit.aes = FALSE,
      linewidth = line_pt,
      colour = axis_col
    ) +
    geom_point(
      data = size_df,
      aes(x = x, y = y, size = value),
      shape = 21,
      fill = "#9D9D9D",
      colour = "#9D9D9D",
      stroke = 0.12
    ) +
    geom_text(data = size_df, aes(x = x, y = 1.26, label = value),
              size = 2.15, family = base_family, colour = axis_col) +
    geom_path(data = ns_ring, aes(x = x, y = y, group = group_id),
              colour = not_deg_col, linewidth = 0.23, linetype = not_deg_linetype, lineend = "round") +
    annotate("text", x = 10.28, y = 1.25, label = "n.s.",
             size = 2.15, family = base_family, colour = axis_col) +
    geom_point(aes(x = 10.82, y = 1.94), shape = 4, size = 3.1, stroke = 0.65, colour = absent_col) +
    annotate("text", x = 10.82, y = 1.25, label = "Absent",
             size = 2.15, family = base_family, colour = axis_col) +
    scale_size_area(max_size = 4.0, limits = c(0, size_limit), guide = "none") +
    coord_cartesian(xlim = c(0, 11.45), ylim = c(0.95, 2.5), expand = FALSE, clip = "off") +
    theme_void(base_family = base_family) +
    theme(plot.margin = margin(0, 2, 0, 12))
}

make_matrix <- function(dat, fc_limit, size_limit) {
  long <- dat %>%
    mutate(
      x = num(x),
      module = factor(module, levels = module_order),
      strain = factor(strain, levels = strain_order),
      logFC = num(logFC),
      neg_log10_p = num(neg_log10_p),
      strain_y = case_when(
        strain == "Halomonas" ~ 4,
        strain == "Psychrobacter" ~ 3,
        strain == "Sphingopyxis" ~ 2,
        strain == "Aequorivita" ~ 1,
        TRUE ~ NA_real_
      )
    )

  module_bounds <- long %>%
    distinct(module, orthogroup_id, x) %>%
    group_by(module) %>%
    summarise(
      n = n_distinct(orthogroup_id),
      x_start = min(x),
      x_end = max(x),
      x_center = (x_start + x_end) / 2,
      x_min = x_start - 0.5,
      x_max = x_end + 0.5,
      module_label = module_labels[as.character(first(module))],
      .groups = "drop"
    ) %>%
    arrange(module)

  og_df <- long %>%
    distinct(module, orthogroup_id, display_gene_label, x) %>%
    mutate(
      gene_label_is_og = grepl("^OG[0-9]+$", display_gene_label),
      gene_label_y = ifelse(gene_label_is_og, 0.49, 0.34)
    ) %>%
    arrange(module, x)

  x_limits <- c(min(module_bounds$x_min) - 0.35, max(module_bounds$x_max) + 0.35)
  y_limits <- c(-3.35, 5.55)
  not_deg_rx <- 0.22
  not_deg_ry <- not_deg_rx * diff(y_limits) / diff(x_limits) * 3.85
  ns_rings <- make_dashed_ellipse(
    long %>% filter(point_class == "Present_not_DEG"),
    x_col = "x",
    y_col = "strain_y",
    rx = not_deg_rx,
    ry = not_deg_ry
  )

  separators <- module_bounds %>%
    slice(seq_len(max(0, nrow(module_bounds) - 1))) %>%
    mutate(x_sep = x_max + 1.25 / 2)

  strain_labels <- tibble(
    strain = strain_order,
    strain_y = c(4, 3, 2, 1)
  )

  p_matrix <- ggplot() +
    geom_rect(
      aes(xmin = x_limits[1], xmax = x_limits[2], ymin = 0.5, ymax = 4.42),
      fill = NA, colour = "#4A4A4A", linewidth = line_pt
    ) +
    geom_vline(
      data = separators,
      aes(xintercept = x_sep),
      colour = "#E2E2E2",
      linewidth = line_pt,
      linetype = "dashed"
    ) +
    geom_path(
      data = ns_rings,
      aes(x = x, y = y, group = group_id),
      colour = not_deg_col,
      linewidth = 0.24,
      linetype = not_deg_linetype,
      lineend = "round"
    ) +
    geom_point(
      data = long %>% filter(point_class == "DEG"),
      aes(x = x, y = strain_y, fill = logFC, size = neg_log10_p),
      shape = 21,
      colour = "#858585",
      stroke = 0.16,
      alpha = 0.96
    ) +
    geom_point(
      data = long %>% filter(point_class == "Absent"),
      aes(x = x, y = strain_y),
      shape = 4,
      size = 2.15,
      stroke = 0.58,
      colour = absent_col
    ) +
    geom_text(
      data = strain_labels,
      aes(x = x_limits[1] - 0.58, y = strain_y, label = strain),
      hjust = 1,
      vjust = 0.5,
      size = 2.65,
      fontface = "italic",
      family = base_family,
      colour = "#3A3A3A"
    ) +
    annotate(
      "text",
      x = x_limits[1] - 0.58,
      y = -0.20,
      label = "Orthogroups",
      hjust = 1,
      vjust = 0.5,
      size = 2.65,
      family = base_family,
      colour = "#000000"
    ) +
    geom_text(
      data = module_bounds,
      aes(x = x_center, y = 5.27, label = module_label),
      size = 1.95,
      lineheight = 0.78,
      family = base_family,
      colour = "#202020"
    ) +
    geom_text(
      data = og_df %>% filter(!gene_label_is_og),
      aes(x = x, y = gene_label_y, label = display_gene_label, colour = module),
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 2.30,
      fontface = "italic",
      family = base_family
    ) +
    geom_text(
      data = og_df %>% filter(gene_label_is_og),
      aes(x = x, y = gene_label_y, label = display_gene_label, colour = module),
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 2.05,
      fontface = "italic",
      family = base_family
    ) +
    geom_segment(
      data = module_bounds,
      aes(x = x_start - 0.42, xend = x_end + 0.42, y = -2.35, yend = -2.35,
          colour = module),
      linewidth = 0.42
    ) +
    geom_text(
      data = module_bounds,
      aes(x = x_center, y = -2.90, label = module_label, colour = module),
      size = 1.92,
      lineheight = 0.84,
      family = base_family,
      fontface = "bold"
    ) +
    scale_colour_manual(values = module_colours, guide = "none") +
    scale_fill_gradient2(
      low = fc_low,
      mid = fc_mid,
      high = fc_high,
      midpoint = 0,
      limits = c(-fc_limit, fc_limit),
      oob = scales::squish
    ) +
    scale_size_area(max_size = 4.0, limits = c(0, size_limit), guide = "none") +
    scale_x_continuous(limits = c(x_limits[1] - 4.2, x_limits[2] + 1.1), expand = c(0, 0)) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void(base_family = base_family) +
    theme(
      legend.position = "none",
      plot.margin = margin(1, 4, 0, 22)
    )

  p_matrix / make_legend_panel(fc_limit, size_limit) + plot_layout(heights = c(0.68, 0.32))
}

dat <- read_tsv(source_file)
validate_source(dat)

dat <- dat %>%
  mutate(
    module = factor(module, levels = module_order),
    x = num(x),
    logFC = num(logFC),
    pvalue = num(pvalue),
    FDR = num(FDR),
    neg_log10_p = num(neg_log10_p),
    member_gene_count_in_strain = num(member_gene_count_in_strain)
  )

fc_limit <- 6
size_limit <- max(20, ceiling(max(dat$neg_log10_p, na.rm = TRUE) / 5) * 5)
fig <- make_matrix(dat, fc_limit, size_limit)

svglite::svglite(output_svg, width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
print(fig)
dev.off()
sanitize_svg_text_spacing(output_svg)

ragg::agg_png(output_png, width = width_px, height = height_px, units = "px", res = 300, background = "white")
print(fig)
dev.off()

qa <- dat %>%
  group_by(module) %>%
  summarise(
    shared_orthogroups = n_distinct(orthogroup_id),
    DEG_cells = sum(point_class == "DEG"),
    present_not_DEG_cells = sum(point_class == "Present_not_DEG"),
    absent_cells = sum(point_class == "Absent"),
    .groups = "drop"
  ) %>%
  mutate(
    figure = figure_id,
    contrast = "D55_vs_D0",
    module = as.character(module),
    output_svg = basename(output_svg),
    output_png = basename(output_png),
    source_data = basename(source_file)
  ) %>%
  select(figure, contrast, module, shared_orthogroups, DEG_cells,
         present_not_DEG_cells, absent_cells, output_svg, output_png, source_data)

write.table(qa, file.path(out_dir, "plot_QA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

package_table <- tibble(
  package = c("dplyr", "ggplot2", "patchwork", "scales", "grid", "svglite", "ragg"),
  role = c(
    "data validation and summaries",
    "matrix plotting",
    "layout composition",
    "colour scale squishing",
    "unit/margin support",
    "SVG export",
    "PNG export"
  ),
  version = vapply(package, function(pkg) {
    if (pkg == "grid") {
      as.character(getRversion())
    } else if (requireNamespace(pkg, quietly = TRUE)) {
      as.character(utils::packageVersion(pkg))
    } else {
      NA_character_
    }
  }, character(1))
)
write.table(package_table, file.path(out_dir, "R_packages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat("FIGURE_J_55D_VS_0P1D_OUTPUT\n")
cat(output_svg, "\n")
cat(output_png, "\n")
cat("\nQA\n")
print(qa)
