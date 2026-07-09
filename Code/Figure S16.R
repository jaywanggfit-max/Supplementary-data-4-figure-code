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

ma <- read_csv(file.path(data_dir, "ma_points.csv"), show_col_types = FALSE) %>%
  mutate(
    strain = factor(strain, levels = c("Halomonas", "Psychrobacter", "Aequorivita", "Sphingopyxis")),
    panel_label = paste0(panel, "  ", comparison_label),
    ma_class_from_code = factor(ma_class_from_code, levels = c("up", "down", "mid"))
  )

p <- ggplot(ma, aes(logCPM, logFC, colour = ma_class_from_code)) +
  geom_point(size = 0.35, alpha = 0.68) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey35") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.25, colour = "grey65") +
  facet_grid(panel_label ~ strain, scales = "free_x") +
  scale_colour_manual(
    values = c(up = "#D95F5F", down = "#4C78A8", mid = "grey78"),
    labels = c(up = expression(log[2] * "FC > 1"), down = expression(log[2] * "FC < -1"), mid = expression(-1 <= log[2] * "FC <= 1"))
  ) +
  labs(x = "Average expression (logCPM)", y = expression(log[2] * "FC"), colour = NULL) +
  theme(legend.position = "top", strip.background = element_blank())

save_figure(p, "Supplementary_Fig16_MA", width_mm = 183, height_mm = 118)
