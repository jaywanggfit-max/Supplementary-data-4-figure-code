suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(stringr)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else ""
script_dir <- if (nzchar(script_file)) {
  dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
raw_file <- if (length(args) >= 1) args[[1]] else file.path(script_dir, "DEG_OG_edgeR.tsv")
out_dir <- if (length(args) >= 2) args[[2]] else script_dir
if (!file.exists(raw_file)) {
  stop("Cannot find input DEG_OG_edgeR.tsv: ", raw_file)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

contrast_titles <- c(
  "D55_vs_B55" = "55D vs 55B",
  "D55_vs_D0" = "55D vs 0.1D"
)

module_order <- c(
  "membrane maintenance",
  "ROS removal",
  "Flagellar assembly",
  "proteostasis",
  "osmoregulation",
  "Aerobic respiration",
  "Translation / ribosome"
)

module_labels <- c(
  "membrane maintenance" = "Membrane\nmaintenance",
  "ROS removal" = "ROS\nremoval",
  "Flagellar assembly" = "Flagellar\nassembly",
  "proteostasis" = "Proteostasis",
  "osmoregulation" = "Osmoregulation",
  "Aerobic respiration" = "Aerobic\nrespiration",
  "Translation / ribosome" = "Translation /\nribosome"
)

module_colours <- c(
  "membrane maintenance" = "#B36F72",
  "ROS removal" = "#B36F72",
  "Flagellar assembly" = "#315B86",
  "proteostasis" = "#B36F72",
  "osmoregulation" = "#B36F72",
  "Aerobic respiration" = "#315B86",
  "Translation / ribosome" = "#315B86"
)

base_family <- "Helvetica"
axis_col <- "#2B2B2B"
grid_col <- "#B5B5B5"
zero_col <- "#6B6B6B"
line_pt <- 0.5 * 0.352778

fill_values <- c(
  "Down | shared" = "#8EA9BD",
  "Down | single" = "#C9D8E4",
  "Up | shared" = "#B77980",
  "Up | single" = "#E0B4B9"
)

legend_labels <- c(
  "Down | shared" = "Down, OG in >=2 strains",
  "Down | single" = "Down, OG in 1 strain",
  "Up | shared" = "Up, OG in >=2 strains",
  "Up | single" = "Up, OG in 1 strain"
)

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

clean_chr <- function(x) {
  ifelse(is.na(x), "", as.character(x))
}

clean_gene_label <- function(gene_label, orthogroup_id) {
  ifelse(
    is.na(gene_label) | gene_label == "" | grepl("^gene[0-9]+$", gene_label),
    sub("^OG0*", "OG", orthogroup_id),
    gene_label
  )
}

sanitize_svg_text_spacing <- function(svg_file) {
  txt <- readLines(svg_file, warn = FALSE, encoding = "UTF-8")
  txt <- gsub("\\s+textLength='[^']*'\\s+lengthAdjust='spacingAndGlyphs'", "", txt)
  writeLines(txt, svg_file, useBytes = TRUE)
}

module_y_table <- tibble(
  module = module_order,
  module_label = module_labels[module_order],
  y = rev(seq_along(module_order))
)

theme_event <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(family = base_family, colour = axis_col),
      axis.line = element_line(linewidth = line_pt, colour = axis_col),
      axis.ticks = element_line(linewidth = line_pt, colour = axis_col),
      axis.title = element_text(size = base_size + 0.25, face = "bold"),
      axis.text = element_text(size = base_size - 0.2),
      plot.title = element_text(size = base_size + 1.25, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.55),
      legend.key.size = unit(2.8, "mm"),
      legend.spacing.x = unit(1.2, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(1.5, 2, 2, 2)
    )
}

raw <- read_tsv(raw_file) %>%
  mutate(
    is_deg = tolower(is_DEG) == "true",
    gene = tolower(OG_gene_symbol),
    logFC = suppressWarnings(as.numeric(logFC)),
    pvalue = suppressWarnings(as.numeric(PValue)),
    FDR = suppressWarnings(as.numeric(FDR)),
    direction = case_when(
      DEG_direction == "up" ~ "Up",
      DEG_direction == "down" ~ "Down",
      logFC > 0 ~ "Up",
      logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    ),
    ko_s = paste0(";", clean_chr(ko_terms), ";"),
    pathway_s = paste0(";", clean_chr(pathway_ids), ";"),
    annotation_l = tolower(paste(
      clean_chr(OG_gene_symbol),
      clean_chr(representative_name),
      clean_chr(representative_product),
      clean_chr(representative_Preferred_name),
      clean_chr(ko_terms),
      clean_chr(pathway_ids),
      clean_chr(pathway_names),
      clean_chr(go_terms),
      sep = " | "
    ))
  ) %>%
  filter(
    is_deg,
    contrast %in% names(contrast_titles),
    !is.na(direction),
    !is.na(orthofinder_orthogroup),
    orthofinder_orthogroup != "",
    !is.na(strain),
    strain != ""
  )

membrane_symbols <- c(
  "bepa", "bama", "bamb", "bamc", "bamd", "bame",
  "pal", "tola", "tolb", "tolq", "tolr",
  "lola", "lolb", "lolc", "lold", "lole", "lolc_e",
  "cfa", "pgpa", "pgpb", "pgpc", "cls", "clsa", "clsb",
  "cdsa", "pssa", "psd", "plsc", "plsb", "plsy", "lgt", "lspa", "lnt"
)
membrane_ko <- ";(K01423|K03634|K03640|K03641|K01095|K00574|K17713);"
membrane_hit <- raw$gene %in% membrane_symbols |
  str_detect(raw$ko_s, membrane_ko) |
  str_detect(
    raw$annotation_l,
    paste(
      "tol-pal",
      "peptidoglycan-associated lipoprotein",
      "outer membrane lipoprotein chaperone",
      "lipoprotein-releasing system",
      "beta-barrel assembly",
      "\\bbam[a-e]\\b",
      "phosphatidylglycerophosphatase",
      "cyclopropane-fatty-acyl-phospholipid synthase",
      sep = "|"
    )
  )
membrane_exclude <- str_detect(
  raw$annotation_l,
  paste(
    "porin",
    "tonb-dependent",
    "antiporter",
    "efflux",
    "phosphate-selective",
    "outer membrane receptor",
    "nutrient uptake",
    sep = "|"
  )
)
membrane_hit <- membrane_hit & !membrane_exclude

ros_symbols <- c(
  "kate", "katg", "soda", "sodb", "sodc",
  "ahpc", "ahpf", "prx", "bcp", "btue", "gpo",
  "osmc", "ohra", "ohr", "tpx"
)
ros_hit <- raw$gene %in% ros_symbols |
  str_detect(
    raw$annotation_l,
    paste(
      "catalase",
      "superoxide dismutase",
      "peroxiredoxin",
      "alkyl hydroperoxide reductase",
      "organic hydroperoxide resistance",
      "glutathione peroxidase",
      "thiol peroxidase",
      sep = "|"
    )
  )
ros_exclude <- str_detect(
  raw$annotation_l,
  paste(
    "transcriptional regulator",
    "\\bregulator\\b",
    "\\bdps\\b",
    "ferritin",
    "methionine sulfoxide",
    "\\bmsra\\b",
    "\\bmsrb\\b",
    "rubrerythrin",
    "glutaredoxin",
    "thioredoxin reductase",
    "cytochrome c-type biogenesis",
    "\\bccmg\\b",
    "\\btrna\\b",
    "threonylcarbamoyladenosine",
    sep = "|"
  )
)
ros_hit <- ros_hit & !ros_exclude

flagellar_hit <- str_detect(raw$gene, "^(flg[b-l]|fli[cdefghijklmnopqrst]|flh[ab]|mot[ab])$") |
  str_detect(raw$pathway_s, ";map02040;")
flagellar_hit <- flagellar_hit & !str_detect(raw$gene, "^(che|aer)")

proteostasis_symbols <- c(
  "groes", "groel", "dnak", "dnaj", "grpe", "htpg",
  "clpa", "clpb", "clpp", "clpx", "lon", "ftsh",
  "hslu", "hslv", "hslo", "htpx", "degp", "mucd",
  "slpa", "fkp", "fkpa", "fkpb", "dsba", "dsbc",
  "ybbn", "hflk", "hflc", "ppia", "ppib", "ppic", "ppid"
)
proteostasis_hit <- raw$gene %in% proteostasis_symbols |
  str_detect(
    raw$annotation_l,
    paste(
      "co-chaperone groes",
      "chaperonin groel",
      "molecular chaperone",
      "atp-dependent.*protease",
      "lon protease",
      "clp protease",
      "hsluv",
      "heat shock.*protease",
      "serine protease do",
      "peptidyl-prolyl.*isomerase",
      "thiol:disulfide interchange",
      "protein disulfide isomerase",
      sep = "|"
    )
  )
proteostasis_exclude <- str_detect(
  raw$annotation_l,
  "dj-1/pfpi|\\bhcha\\b|generic peptidase|metalloprotease$"
)
proteostasis_hit <- proteostasis_hit & !proteostasis_exclude

osmo_symbols <- c(
  "ecta", "ectb", "ectc", "ectd", "beta", "betb", "bett", "betl",
  "proi", "prop", "prov", "prow", "prox", "opua", "opub", "opuc", "opud",
  "mscl", "msck", "mscs", "yggb"
)
osmo_hit <- raw$gene %in% osmo_symbols |
  str_detect(
    raw$annotation_l,
    paste(
      "ectoine",
      "hydroxyectoine",
      "glycine betaine",
      "choline dehydrogenase",
      "betaine aldehyde dehydrogenase",
      "osmoprotectant",
      "compatible solute",
      "mechanosensitive channel",
      sep = "|"
    )
  )
osmo_exclude <- str_detect(
  raw$annotation_l,
  "transcriptional regulator|\\bbeti\\b|antiporter|\\bnha\\b|\\bpha[a-z]\\b|\\bmnh[a-z]\\b|\\bcvra\\b"
)
osmo_hit <- osmo_hit & !osmo_exclude

resp_ko <- paste0(
  ";(K0033[0-9]|K00340|K0034[67]|K00239|K0024[0-2]|",
  "K0040[4-7]|K0041[1-3]|K0042[56]|K0227[4-6]|K0229[7-9]|",
  "K02300|K0210[8-9]|K0211[0-5]);"
)
respiration_hit <- str_detect(
  raw$gene,
  "^(cco(no|[nopq])|cox(ac|[abc])|cyd[ab]|cyo[ad]|pet[abc]|sdh[abcd]|nuo[hijk]|nqr[ab]|atp[a-h])$"
) | str_detect(raw$ko_s, resp_ko)

translation_hit <- (
  str_detect(raw$gene, "^(rps|rpl|rpm)[a-z0-9]+$") |
    raw$gene %in% c(
      "infa", "infb", "infc", "fusa", "tuf", "tufa", "tufb", "tsf",
      "efp", "prfa", "prfb", "prfc", "frr", "lepa", "fmt"
    ) |
    str_detect(
      raw$annotation_l,
      paste(
        "map03010",
        "ribosomal protein",
        "translation initiation factor",
        "translation elongation factor",
        "elongation factor p",
        "ribosome recycling factor",
        "peptide chain release factor",
        sep = "|"
      )
    )
) & !str_detect(
  raw$annotation_l,
  "transcription elongation|grea/greb|release factor.*methyltransferase"
)

hit_table <- tibble(
  rowid = seq_len(nrow(raw)),
  `membrane maintenance` = membrane_hit,
  `ROS removal` = ros_hit,
  `Flagellar assembly` = flagellar_hit,
  proteostasis = proteostasis_hit,
  osmoregulation = osmo_hit,
  `Aerobic respiration` = respiration_hit,
  `Translation / ribosome` = translation_hit
) %>%
  pivot_longer(-rowid, names_to = "module", values_to = "hit") %>%
  filter(hit)

priority <- c(
  "membrane maintenance" = 1,
  "ROS removal" = 2,
  "Flagellar assembly" = 3,
  "proteostasis" = 4,
  "osmoregulation" = 5,
  "Aerobic respiration" = 6,
  "Translation / ribosome" = 7
)

source_data <- hit_table %>%
  mutate(priority = priority[module]) %>%
  group_by(rowid) %>%
  arrange(priority, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  bind_cols(raw[.$rowid, ]) %>%
  mutate(
    module = factor(module, levels = module_order),
    orthogroup_id = orthofinder_orthogroup,
    gene_label = clean_gene_label(OG_gene_symbol, orthogroup_id),
    strict_function_evidence = case_when(
      module == "membrane maintenance" ~
        "strict membrane maintenance: Tol-Pal/envelope integrity, lipid synthesis/modification, or OM protein/lipoprotein assembly",
      module == "ROS removal" ~
        "strict ROS removal: catalase/SOD/Ahp/peroxiredoxin/OsmC-Ohr/Bcp/BtuE/Tpx direct detoxification enzyme",
      module == "Flagellar assembly" ~
        "strict flagellar assembly: flagellar structural, motor, basal-body, hook, or export component",
      module == "proteostasis" ~
        "strict proteostasis: chaperone/co-chaperone, folding/isomerization factor, or stress/ATP-dependent protease",
      module == "osmoregulation" ~
        "strict osmoregulation: compatible-solute synthesis/transport or mechanosensitive channel",
      module == "Aerobic respiration" ~
        "strict aerobic respiration: respiratory-chain or oxidative-phosphorylation complex subunit",
      module == "Translation / ribosome" ~
        "strict translation/ribosome: ribosomal structural protein or core translation factor",
      TRUE ~ ""
    ),
    source_table = "DEG_OG_edgeR.tsv strict full-module filter"
  ) %>%
  select(
    contrast, strain, module, orthogroup_id, gene_label, direction, logFC, pvalue, FDR,
    n_DEG_strains_raw = n_member_strains,
    strict_function_evidence, source_table,
    ko_terms, pathway_ids, pathway_names,
    representative_product, representative_Preferred_name,
    member_gene_symbols, member_products
  ) %>%
  arrange(contrast, module, orthogroup_id, strain)

og_scope <- source_data %>%
  distinct(contrast, module, orthogroup_id, strain) %>%
  count(contrast, module, orthogroup_id, name = "n_DEG_strains") %>%
  mutate(
    sharing = if_else(n_DEG_strains >= 2, "shared", "single"),
    sharing_label = if_else(n_DEG_strains >= 2, "OG in >=2 strains", "OG in 1 strain")
  )

event_data <- source_data %>%
  left_join(og_scope, by = c("contrast", "module", "orthogroup_id")) %>%
  mutate(
    module = factor(module, levels = module_order),
    direction = factor(direction, levels = c("Down", "Up")),
    sharing = factor(sharing, levels = c("shared", "single"))
  ) %>%
  arrange(contrast, module, orthogroup_id, strain)

event_counts <- event_data %>%
  count(contrast, module, direction, sharing, name = "n") %>%
  complete(
    contrast = names(contrast_titles),
    module = factor(module_order, levels = module_order),
    direction = factor(c("Down", "Up"), levels = c("Down", "Up")),
    sharing = factor(c("shared", "single"), levels = c("shared", "single")),
    fill = list(n = 0)
  ) %>%
  left_join(module_y_table, by = c("module" = "module")) %>%
  group_by(contrast, module, direction) %>%
  arrange(sharing, .by_group = TRUE) %>%
  mutate(
    cum_n = cumsum(n),
    prev_n = lag(cum_n, default = 0),
    xmin = if_else(direction == "Up", prev_n, -cum_n),
    xmax = if_else(direction == "Up", cum_n, -prev_n),
    ymin = y - 0.33,
    ymax = y + 0.33,
    fill_key = paste(as.character(direction), "|", as.character(sharing))
  ) %>%
  ungroup()

direction_totals <- event_counts %>%
  group_by(contrast, module, module_label, direction, y) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    signed_n = if_else(direction == "Down", -n, n),
    label_x = signed_n + if_else(direction == "Down", -2.2, 2.2),
    label_hjust = if_else(direction == "Down", 1, 0)
  )

summary_counts <- event_counts %>%
  group_by(contrast, module, direction, sharing) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    module = as.character(module),
    direction = as.character(direction),
    sharing = as.character(sharing),
    count_unit = "OG-strain DEG event"
  ) %>%
  arrange(contrast, factor(module, levels = module_order), direction, sharing)

wide_summary <- event_data %>%
  count(contrast, module, direction, name = "n") %>%
  mutate(direction = as.character(direction)) %>%
  complete(
    contrast = names(contrast_titles),
    module = factor(module_order, levels = module_order),
    direction = c("Down", "Up"),
    fill = list(n = 0)
  ) %>%
  pivot_wider(names_from = direction, values_from = n, names_prefix = "event_") %>%
  left_join(
    event_data %>%
      distinct(contrast, module, orthogroup_id, n_DEG_strains, sharing) %>%
      count(contrast, module, sharing, name = "unique_OGs") %>%
      mutate(sharing = as.character(sharing)) %>%
      complete(
        contrast = names(contrast_titles),
        module = factor(module_order, levels = module_order),
        sharing = c("shared", "single"),
        fill = list(unique_OGs = 0)
      ) %>%
      pivot_wider(names_from = sharing, values_from = unique_OGs, names_prefix = "unique_OG_"),
    by = c("contrast", "module")
  ) %>%
  mutate(
    event_total = event_Down + event_Up,
    unique_OG_total = unique_OG_shared + unique_OG_single
  ) %>%
  arrange(contrast, factor(module, levels = module_order))

review_table <- event_data %>%
  group_by(contrast, module, orthogroup_id, gene_label) %>%
  summarise(
    n_DEG_events = n(),
    n_DEG_strains = first(n_DEG_strains),
    sharing = first(as.character(sharing)),
    DEG_strains = paste(sort(unique(strain)), collapse = ";"),
    directions = paste(sort(unique(as.character(direction))), collapse = ";"),
    strict_function_evidence = first(strict_function_evidence),
    representative_product = first(representative_product),
    ko_terms = first(ko_terms),
    pathway_ids = first(pathway_ids),
    member_gene_symbols = first(member_gene_symbols),
    .groups = "drop"
  ) %>%
  arrange(contrast, factor(module, levels = module_order), orthogroup_id)

write.table(
  event_data %>% mutate(module = as.character(module), direction = as.character(direction), sharing = as.character(sharing)),
  file.path(out_dir, "strict_extended_module_strain_event_source_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  summary_counts,
  file.path(out_dir, "strict_extended_module_strain_event_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  wide_summary %>% mutate(module = as.character(module)),
  file.path(out_dir, "strict_extended_module_strain_event_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

write.table(
  review_table %>% mutate(module = as.character(module)),
  file.path(out_dir, "strict_extended_module_strain_event_retained_og_review.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

qa <- tibble(
  check = c(
    "raw_DEG_rows_in_two_contrasts",
    "retained_strict_function_events",
    "counted_strain_events",
    "unique_contrast_module_OGs",
    "shared_unique_OGs",
    "single_strain_unique_OGs",
    "shared_strain_events",
    "single_strain_events",
    "D55_vs_B55_ROS_events",
    "D55_vs_B55_ROS_single_events"
  ),
  value = c(
    nrow(raw),
    nrow(event_data),
    sum(summary_counts$n),
    nrow(og_scope),
    sum(og_scope$sharing == "shared"),
    sum(og_scope$sharing == "single"),
    sum(summary_counts$n[summary_counts$sharing == "shared"]),
    sum(summary_counts$n[summary_counts$sharing == "single"]),
    sum(summary_counts$n[summary_counts$contrast == "D55_vs_B55" & summary_counts$module == "ROS removal"]),
    sum(summary_counts$n[summary_counts$contrast == "D55_vs_B55" & summary_counts$module == "ROS removal" & summary_counts$sharing == "single"])
  )
)

write.table(
  qa,
  file.path(out_dir, "strict_extended_module_strain_event_counts_QA.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

target_contrast <- "D55_vs_B55"
plot_summary <- wide_summary %>%
  filter(contrast == target_contrast) %>%
  mutate(
    module = factor(module, levels = module_order),
    y = rev(seq_along(module_order))
  )

bar_df <- bind_rows(
  plot_summary %>%
    transmute(module, y, direction = "Down", value = event_Down,
              xmin = -event_Down, xmax = 0),
  plot_summary %>%
    transmute(module, y, direction = "Up", value = event_Up,
              xmin = 0, xmax = event_Up)
) %>%
  filter(value > 0)

label_df <- bar_df %>%
  mutate(
    label = value,
    x = if_else(direction == "Down", xmin - 2.2, xmax + 2.2),
    hjust = if_else(direction == "Down", 1, 0)
  )

x_lim <- 110
x_breaks <- c(-100, -50, 0, 50, 100)

make_label_panel_final <- function() {
  ggplot(module_y_table, aes(y = y)) +
    annotate(
      "text",
      x = 0.03,
      y = length(module_order) + 0.46,
      label = "i",
      hjust = 0,
      family = base_family,
      fontface = "bold",
      size = 9.2,
      colour = "#222222"
    ) +
    geom_text(
      aes(x = 0.97, label = module_label, colour = module),
      hjust = 1,
      lineheight = 0.82,
      family = base_family,
      size = 3.95,
      show.legend = FALSE
    ) +
    geom_segment(
      aes(x = 1.005, xend = 1.06, yend = y),
      linewidth = line_pt,
      colour = axis_col
    ) +
    scale_colour_manual(values = module_colours, guide = "none") +
    scale_x_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, length(module_order) + 0.75), expand = c(0, 0)) +
    theme_void(base_family = base_family) +
    theme(plot.margin = margin(0, 0, 12, 0))
}

make_bar_panel_final <- function() {
  ggplot() +
    geom_vline(
      xintercept = x_breaks,
      linewidth = line_pt,
      linetype = "longdash",
      colour = grid_col
    ) +
    geom_vline(xintercept = 0, linewidth = line_pt, colour = zero_col) +
    geom_rect(
      data = bar_df,
      aes(xmin = xmin, xmax = xmax, ymin = y - 0.32, ymax = y + 0.32, fill = direction),
      colour = axis_col,
      linewidth = line_pt
    ) +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label, hjust = hjust, colour = direction),
      family = base_family,
      fontface = "bold",
      size = 3.15,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = -72,
      y = length(module_order) + 0.55,
      label = "Down",
      family = base_family,
      fontface = "bold",
      size = 4.05,
      colour = axis_col
    ) +
    annotate(
      "text",
      x = 76,
      y = length(module_order) + 0.55,
      label = "Up",
      family = base_family,
      fontface = "bold",
      size = 4.05,
      colour = axis_col
    ) +
    scale_fill_manual(values = c("Down" = "#8EA9BD", "Up" = "#B77980"), guide = "none") +
    scale_colour_manual(values = c("Down" = "#315B86", "Up" = "#B36F72"), guide = "none") +
    scale_x_continuous(
      limits = c(-x_lim, x_lim),
      breaks = x_breaks,
      labels = abs(x_breaks),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0.5, length(module_order) + 0.75),
      breaks = module_y_table$y,
      labels = rep("", length(module_order)),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(x = "DE orthogroups", y = NULL) +
    theme_classic(base_size = 7.4, base_family = base_family) +
    theme(
      text = element_text(family = base_family, colour = axis_col),
      axis.line = element_line(linewidth = line_pt, colour = axis_col),
      axis.ticks = element_line(linewidth = line_pt, colour = axis_col),
      axis.title.x = element_text(size = 9.2, face = "bold", margin = margin(t = 6)),
      axis.text.x = element_text(size = 7.7),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 5, 12, 0)
    )
}

make_figure_i <- function() {
  make_label_panel_final() | make_bar_panel_final() +
    plot_layout(widths = c(0.31, 0.69))
}

save_svg <- function(plot, stem, width_mm, height_mm) {
  svg_file <- file.path(out_dir, paste0(stem, ".svg"))

  svglite::svglite(svg_file, width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
  print(plot)
  dev.off()
  sanitize_svg_text_spacing(svg_file)

  tibble(svg = svg_file)
}

fig_files <- save_svg(
  make_figure_i(),
  "Figure_i_DE_orthogroups_55D_vs_55B",
  150,
  88
) %>%
  mutate(figure = "Figure_i_55D_vs_55B")

cat("INPUT_RAW_FILE\n")
cat(raw_file, "\n\n")
cat("FIGURE_FILES\n")
print(fig_files)
cat("\nSTRAIN_EVENT_SUMMARY\n")
print(wide_summary, n = Inf)
cat("\nQA\n")
print(qa, n = Inf)
