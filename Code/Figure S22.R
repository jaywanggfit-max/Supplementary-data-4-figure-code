required_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork", "scales", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- if (length(file_arg) > 0) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("plot.R")
base_dir <- dirname(script_path)
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_csv <- function(file, show_col_types = FALSE) {
  read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
}

theme_set(
  theme_classic(base_size = 7, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      strip.text = element_text(size = 6.5, face = "italic"),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      panel.border = element_rect(fill = NA, linewidth = 0.3, colour = "black"),
      panel.grid = element_blank()
    )
)

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

points <- read_csv(file.path(data_dir, "volcano_points.csv"), show_col_types = FALSE) %>%
  mutate(
    strain = factor(strain, levels = c("Halomonas", "Psychrobacter", "Aequorivita", "Sphingopyxis")),
    panel_label = paste0(panel, "  ", comparison_label),
    module_for_plot = factor(module_for_plot, levels = c("Other", "Proteostasis", "Membrane maintenance", "ROS removal", "Flagellar / chemotaxis", "Osmoregulation"))
  )
labels <- read_csv(file.path(data_dir, "selected_labels.csv"), show_col_types = FALSE) %>%
  mutate(
    strain = factor(strain, levels = c("Halomonas", "Psychrobacter", "Aequorivita", "Sphingopyxis")),
    panel_label = paste0(panel, "  ", comparison_label),
    module_for_plot = factor(module_for_plot, levels = levels(points$module_for_plot))
  )

module_cols <- c(
  "Other" = "grey78",
  "Proteostasis" = "#77A35D",
  "Membrane maintenance" = "#C99A45",
  "ROS removal" = "#54B6B0",
  "Flagellar / chemotaxis" = "#9A6BB7",
  "Osmoregulation" = "#6AA2D8"
)

p <- ggplot(points, aes(logFC, minus_log10_FDR)) +
  geom_point(data = points %>% filter(module_for_plot == "Other"), colour = "grey78", size = 0.35, alpha = 0.55) +
  geom_point(data = points %>% filter(module_for_plot != "Other"), aes(colour = module_for_plot), size = 0.7, alpha = 0.86) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.22, colour = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.22, colour = "grey55") +
  geom_text(data = labels, aes(label = label_for_plot, colour = module_for_plot), size = 1.7, check_overlap = TRUE, show.legend = FALSE) +
  facet_grid(panel_label ~ strain, scales = "free") +
  scale_colour_manual(values = module_cols, breaks = names(module_cols)[-1]) +
  labs(x = expression(log[2] * "FC"), y = expression(-log[10] * "(FDR)"), colour = NULL) +
  theme(legend.position = "top", strip.background = element_blank())

save_figure(p, "Supplementary_Fig22_Volcano", width_mm = 183, height_mm = 150)
