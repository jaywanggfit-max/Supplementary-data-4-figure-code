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

totals <- read_csv(file.path(data_dir, "strain_totals.csv"), show_col_types = FALSE)
inter <- read_csv(file.path(data_dir, "total_intersections.csv"), show_col_types = FALSE) %>%
  arrange(desc(exclusive_total_intersection_size)) %>%
  mutate(intersection_id = row_number())

strain_levels <- c("Halomonas", "Psychrobacter", "Aequorivita", "Sphingopyxis")
strain_cols <- c(Halomonas = "#F2B179", Psychrobacter = "#86B7D7", Aequorivita = "#D6A0C8", Sphingopyxis = "#66B89C")
direction_cols <- c(up = "#B96C75", down = "#9BB0C5")

bar_df <- totals %>%
  select(strain, up, down) %>%
  pivot_longer(c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(count_plot = if_else(direction == "down", -count, count), strain = factor(strain, levels = rev(strain_levels)))

p_left <- ggplot(bar_df, aes(count_plot, strain, fill = direction)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linewidth = 0.25) +
  geom_text(aes(label = count, hjust = if_else(direction == "down", 1.08, -0.08)), size = 2.0) +
  scale_fill_manual(values = direction_cols) +
  scale_x_continuous(labels = abs, expand = expansion(mult = 0.12)) +
  labs(x = "DE orthogroups per strain", y = NULL, fill = NULL) +
  theme(legend.position = "top")

p_top <- ggplot(inter, aes(intersection_id, exclusive_total_intersection_size)) +
  geom_col(width = 0.72, fill = "grey78", colour = "grey35", linewidth = 0.15) +
  geom_text(aes(label = exclusive_total_intersection_size), vjust = -0.35, size = 1.9) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = NULL, y = "Shared DE orthogroups") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

matrix_df <- inter %>%
  select(intersection_id, ends_with("_included")) %>%
  pivot_longer(ends_with("_included"), names_to = "strain", values_to = "included") %>%
  mutate(
    included = included %in% c(TRUE, "TRUE", "True", "true", "1", 1),
    strain = sub("_included$", "", strain),
    strain = factor(strain, levels = rev(strain_levels)),
    strain_num = as.numeric(strain)
  )
segments <- matrix_df %>%
  filter(included) %>%
  group_by(intersection_id) %>%
  summarise(ymin = min(strain_num), ymax = max(strain_num), .groups = "drop")

p_matrix <- ggplot(matrix_df, aes(intersection_id, strain_num)) +
  geom_segment(data = segments, aes(x = intersection_id, xend = intersection_id, y = ymin, yend = ymax), inherit.aes = FALSE, linewidth = 0.35, colour = "grey30") +
  geom_point(aes(fill = strain, alpha = included), shape = 21, size = 2.0, colour = "grey30", stroke = 0.15) +
  scale_fill_manual(values = strain_cols) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.18), guide = "none") +
  scale_y_continuous(breaks = seq_along(rev(strain_levels)), labels = rev(strain_levels)) +
  labs(x = "Intersection", y = NULL) +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())

circle_df <- function(name, x0, y0, r = 1, n = 180) {
  theta <- seq(0, 2*pi, length.out = n)
  tibble(strain = name, x = x0 + r*cos(theta), y = y0 + r*sin(theta))
}
venn_circles <- bind_rows(
  circle_df("Halomonas", -0.6, 0.15),
  circle_df("Psychrobacter", 0.0, 0.45),
  circle_df("Aequorivita", 0.55, 0.15),
  circle_df("Sphingopyxis", 0.0, -0.35)
)
venn_labels <- inter %>%
  mutate(
    x = case_when(
      intersection_pattern == "Halomonas" ~ -1.15,
      intersection_pattern == "Psychrobacter" ~ 0.0,
      intersection_pattern == "Aequorivita" ~ 1.10,
      intersection_pattern == "Sphingopyxis" ~ 0.0,
      intersection_pattern == "Halomonas&Psychrobacter" ~ -0.45,
      intersection_pattern == "Halomonas&Aequorivita" ~ -0.05,
      intersection_pattern == "Halomonas&Sphingopyxis" ~ -0.45,
      intersection_pattern == "Psychrobacter&Aequorivita" ~ 0.45,
      intersection_pattern == "Psychrobacter&Sphingopyxis" ~ 0.30,
      intersection_pattern == "Aequorivita&Sphingopyxis" ~ 0.55,
      intersection_pattern == "Halomonas&Psychrobacter&Aequorivita" ~ 0.05,
      intersection_pattern == "Halomonas&Psychrobacter&Sphingopyxis" ~ -0.15,
      intersection_pattern == "Halomonas&Aequorivita&Sphingopyxis" ~ 0.12,
      intersection_pattern == "Psychrobacter&Aequorivita&Sphingopyxis" ~ 0.28,
      TRUE ~ 0.05
    ),
    y = case_when(
      intersection_pattern == "Halomonas" ~ 0.10,
      intersection_pattern == "Psychrobacter" ~ 1.18,
      intersection_pattern == "Aequorivita" ~ 0.10,
      intersection_pattern == "Sphingopyxis" ~ -1.10,
      intersection_pattern == "Halomonas&Psychrobacter" ~ 0.62,
      intersection_pattern == "Halomonas&Aequorivita" ~ 0.25,
      intersection_pattern == "Halomonas&Sphingopyxis" ~ -0.45,
      intersection_pattern == "Psychrobacter&Aequorivita" ~ 0.62,
      intersection_pattern == "Psychrobacter&Sphingopyxis" ~ -0.18,
      intersection_pattern == "Aequorivita&Sphingopyxis" ~ -0.42,
      intersection_pattern == "Halomonas&Psychrobacter&Aequorivita" ~ 0.58,
      intersection_pattern == "Halomonas&Psychrobacter&Sphingopyxis" ~ -0.02,
      intersection_pattern == "Halomonas&Aequorivita&Sphingopyxis" ~ -0.52,
      intersection_pattern == "Psychrobacter&Aequorivita&Sphingopyxis" ~ -0.05,
      TRUE ~ 0.02
    )
  )
p_venn <- ggplot() +
  geom_polygon(data = venn_circles, aes(x, y, group = strain, fill = strain), alpha = 0.33, colour = "grey35", linewidth = 0.25) +
  geom_text(data = venn_labels, aes(x, y, label = exclusive_total_intersection_size), size = 1.9) +
  scale_fill_manual(values = strain_cols) +
  coord_equal(xlim = c(-1.7, 1.7), ylim = c(-1.45, 1.55), expand = FALSE) +
  theme_void(base_size = 7) +
  theme(legend.position = "none")

fig <- (p_venn | (p_top / p_matrix + plot_layout(heights = c(1.15, 0.85)))) / p_left +
  plot_layout(heights = c(2, 0.9))
save_figure(fig, "Supplementary_Fig19_Shared_DEG", width_mm = 170, height_mm = 130)
