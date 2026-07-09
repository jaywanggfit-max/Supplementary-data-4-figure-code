required_packages <- c("ggplot2", "dplyr", "patchwork", "scales", "svglite", "ragg")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing) > 0) {
  stop("Missing R packages: ", paste(still_missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    script_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(script_path)) {
      return(dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

base_dir <- get_script_dir()
data_dir <- file.path(base_dir, "data")
out_dir <- file.path(base_dir, "output")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw_file <- file.path(data_dir, "kegg_dotplot_latest_raw.tsv")
plot_file <- file.path(data_dir, "kegg_dotplot.csv")
qa_file <- file.path(data_dir, "kegg_dotplot_update_QA.tsv")

if (!file.exists(raw_file)) {
  stop("Cannot find raw data file: ", raw_file)
}

lines <- readLines(raw_file, warn = FALSE, encoding = "UTF-8")
lines <- lines[nzchar(gsub("\t", "", trimws(lines)))]
lines <- lines[!grepl("^contrast\\tterm_id\\tterm_name\\tGeneRatio\\tCount\\t", lines)]

header <- "contrast\tterm_id\tterm_name\tGeneRatio\tCount\tpvalue\tp.adjust\tneg_log10_p_adjust\tgeneID"
clean_text <- paste(c(header, lines), collapse = "\n")
raw <- read.delim(
  text = clean_text,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

label_map <- c(
  "D55_vs_B55" = "KEGG enrichment: 55D vs 55B",
  "D55_vs_D0" = "KEGG enrichment: 55D vs 0D"
)
contrast_order <- c("D55_vs_B55", "D55_vs_D0")

plot_data <- raw %>%
  mutate(
    input_order = row_number(),
    GeneRatio = as.numeric(GeneRatio),
    Count = as.integer(Count),
    pvalue = as.numeric(pvalue),
    p.adjust = as.numeric(p.adjust),
    neg_log10_p_adjust = as.numeric(neg_log10_p_adjust),
    comparison_label = unname(label_map[contrast])
  ) %>%
  filter(contrast %in% contrast_order, !is.na(comparison_label)) %>%
  group_by(contrast) %>%
  arrange(desc(GeneRatio), input_order, .by_group = TRUE) %>%
  mutate(plot_rank = row_number()) %>%
  ungroup() %>%
  mutate(contrast = factor(contrast, levels = contrast_order)) %>%
  arrange(contrast, plot_rank) %>%
  mutate(contrast = as.character(contrast)) %>%
  select(
    contrast, comparison_label, plot_rank, term_id, term_name,
    GeneRatio, Count, pvalue, p.adjust, neg_log10_p_adjust, geneID
  )

if (nrow(plot_data) != 30) {
  stop("Expected 30 KEGG rows after cleaning; observed ", nrow(plot_data), ".")
}
if (any(is.na(plot_data$GeneRatio)) || any(is.na(plot_data$p.adjust))) {
  stop("Numeric parsing failed for GeneRatio or p.adjust.")
}
if (any(abs(-log10(plot_data$p.adjust) - plot_data$neg_log10_p_adjust) > 1e-4)) {
  stop("Provided -log10(p.adjust) values do not match p.adjust within tolerance.")
}

write.csv(plot_data, plot_file, row.names = FALSE, quote = TRUE)

qa <- plot_data %>%
  group_by(contrast, comparison_label) %>%
  summarise(
    rows = n(),
    min_GeneRatio = min(GeneRatio),
    max_GeneRatio = max(GeneRatio),
    min_p_adjust = min(p.adjust),
    max_p_adjust = max(p.adjust),
    .groups = "drop"
  )
write.table(qa, qa_file, sep = "\t", quote = FALSE, row.names = FALSE)

theme_set(
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.title = element_text(size = 8, face = "bold"),
      axis.text = element_text(size = 7, colour = "black"),
      legend.title = element_text(size = 5.8, face = "bold"),
      legend.text = element_text(size = 6.0, colour = "black"),
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
      panel.border = element_rect(fill = NA, linewidth = 0.35, colour = "black"),
      panel.grid = element_blank()
    )
)

x_label <- function(x) ifelse(abs(x) < 1e-9, "0", sprintf("%.2f", x))

make_panel <- function(panel_df, title, fill_breaks) {
  panel_df <- panel_df %>%
    arrange(desc(GeneRatio), plot_rank) %>%
    mutate(term_name = factor(term_name, levels = rev(term_name)))

  ggplot(panel_df, aes(GeneRatio, term_name)) +
    geom_point(
      aes(size = Count, fill = neg_log10_p_adjust),
      shape = 21,
      colour = "black",
      stroke = 0.32,
      alpha = 0.95
    ) +
    scale_x_continuous(
      limits = c(0, 0.10),
      breaks = seq(0, 0.10, by = 0.02),
      labels = x_label,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_fill_gradient(
      low = "#2C4A9A",
      high = "#E7292A",
      limits = range(c(0, fill_breaks)),
      breaks = fill_breaks,
      labels = sprintf("%.1f", fill_breaks),
      oob = scales::squish,
      name = expression("-log"[10] * "(P.adjust)")
    ) +
    scale_size_continuous(
      range = c(1.4, 4.0),
      limits = c(0, 50),
      breaks = c(10, 20, 30, 40),
      name = "Count"
    ) +
    guides(
      fill = guide_colourbar(
        order = 1,
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(0.13, "in"),
        barheight = unit(0.42, "in"),
        ticks = TRUE,
        frame.colour = "black"
      ),
      size = guide_legend(
        order = 2,
        title.position = "top",
        keywidth = unit(0.38, "cm"),
        keyheight = unit(0.28, "cm"),
        override.aes = list(fill = "#B8B8B8", colour = "#666666", alpha = 1)
      )
    ) +
    labs(x = "GeneRatio", y = "KEGG pathway", title = title) +
    theme(
      legend.position = c(0.56, 0.05),
      legend.justification = c(0, 0),
      legend.box = "vertical",
      legend.box.just = "left",
      legend.spacing.y = unit(-0.03, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.key.width = unit(0.38, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.background = element_blank(),
      axis.text.y = element_text(size = 7.5),
      plot.margin = margin(6, 8, 6, 8)
    )
}

p_a <- make_panel(
  plot_data %>% filter(contrast == "D55_vs_B55"),
  "KEGG enrichment: 55D vs 55B",
  1:4
)

p_b <- make_panel(
  plot_data %>% filter(contrast == "D55_vs_D0"),
  "KEGG enrichment: 55D vs 0D",
  1:5
)

p <- (p_a | p_b) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 17, face = "plain", colour = "black"),
    plot.tag.position = c(0.01, 1.02)
  )

save_figure <- function(plot, stem, width_mm, height_mm, dpi = 450) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svg_file <- file.path(out_dir, paste0(stem, ".svg"))
  pdf_file <- file.path(out_dir, paste0(stem, ".pdf"))
  png_file <- file.path(out_dir, paste0(stem, ".png"))

  svglite::svglite(svg_file, width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(pdf_file, width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_png(png_file, width = w, height = h, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()

  message("Saved processed data: ", plot_file)
  message("Saved QA summary: ", qa_file)
  message("Exported figure: ", png_file)
  message("Exported figure: ", svg_file)
  message("Exported figure: ", pdf_file)
}

save_figure(p, "S20_KEGG", width_mm = 190, height_mm = 90)
