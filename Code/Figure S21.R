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

contrast_levels <- c("D55_vs_D0", "D55_vs_B55")
df <- read_csv(file.path(data_dir, "go_dotplot.csv"), show_col_types = FALSE) %>%
  mutate(contrast = factor(contrast, levels = contrast_levels)) %>%
  arrange(contrast, plot_rank) %>%
  mutate(
    comparison_label = factor(comparison_label, levels = unique(comparison_label)),
    term_label = if_else(is.na(term_wrapped) | term_wrapped == "", term_name, term_wrapped),
    term_plot = paste(term_label, contrast, sep = "___"),
    term_plot = factor(term_plot, levels = rev(unique(term_plot)))
  )
lab_fun <- function(x) sub("___.*$", "", x)

p <- ggplot(df, aes(GeneRatio, term_plot)) +
  geom_point(aes(size = Count, colour = neg_log10_p_adjust), alpha = 0.92) +
  facet_wrap(~ comparison_label, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = lab_fun) +
  scale_colour_gradient(low = "#2C4A9A", high = "#E7292A", name = expression(-log[10] * "(P.adjust)")) +
  scale_size_continuous(range = c(2.0, 6.0)) +
  labs(x = "GeneRatio", y = "Biological process", size = "Count") +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 5.2, face = "bold", margin = margin(1, 1, 1, 1))
  )

save_figure(p, "Supplementary_Fig21_GO", width_mm = 200, height_mm = 118)
