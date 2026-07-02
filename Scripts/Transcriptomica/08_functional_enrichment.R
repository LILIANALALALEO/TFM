if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "tximport",
  "DESeq2",
  "clusterProfiler",
  "org.At.tair.db",
  "AnnotationDbi",
  "enrichplot"
))

install.packages(c("dplyr", "readr", "stringr", "tidyr", "ggplot2"))

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tidyr)
  library(stringr)
  library(scales)
})

# =====================================================
# CONFIGURACIÓN
# =====================================================

deseq_dir <- "Data/Processed/04.DESeq2/cold"
curated_file <- "Results/cold/Estres frio_TFM.xlsx"
outdir <- "Results/cold_figures_publication"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05
lfc_cutoff <- -1
key_contrast <- "24h_vs_4h"

files <- list(
  "4h_vs_control" = file.path(deseq_dir, "cold_4h_vs_control.csv"),
  "24h_vs_control" = file.path(deseq_dir, "cold_24h_vs_control.csv"),
  "24h_vs_4h" = file.path(deseq_dir, "cold_24h_vs_cold_4h.csv")
)

# =====================================================
# VALIDAR ARCHIVOS
# =====================================================

missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop("Faltan archivos DESeq2:\n", paste(missing_files, collapse = "\n"))
}

if (!file.exists(curated_file)) {
  stop("No existe el archivo curado: ", curated_file)
}

# =====================================================
# LEER TABLA CURADA
# =====================================================

sheets <- excel_sheets(curated_file)
curated <- read_excel(curated_file, sheet = sheets[1])
names(curated) <- trimws(names(curated))

find_col <- function(df, patterns) {
  nms <- names(df)
  hit <- nms[str_detect(tolower(nms), paste(patterns, collapse = "|"))]
  if (length(hit) == 0) return(NA)
  hit[1]
}

gene_col <- find_col(curated, c("gene", "gene id", "geneid", "id"))
tipo_col <- find_col(curated, c("tipo", "categoria", "relacion", "relación"))
bp_col <- find_col(curated, c("bp", "biological", "proceso"))
mf_col <- find_col(curated, c("mf", "molecular", "funcion", "función"))
kegg_col <- find_col(curated, c("kegg", "ruta"))
prioridad_col <- find_col(curated, c("prioridad", "priority"))
rol_col <- find_col(curated, c("rol", "role"))

if (is.na(gene_col)) {
  stop("No se detectó columna de Gene ID en el archivo curado.")
}

curated2 <- curated %>%
  mutate(
    GeneID = trimws(as.character(.data[[gene_col]])),
    Tipo = if (!is.na(tipo_col)) as.character(.data[[tipo_col]]) else NA_character_,
    BP = if (!is.na(bp_col)) as.character(.data[[bp_col]]) else NA_character_,
    MF = if (!is.na(mf_col)) as.character(.data[[mf_col]]) else NA_character_,
    KEGG = if (!is.na(kegg_col)) as.character(.data[[kegg_col]]) else NA_character_,
    Prioridad = if (!is.na(prioridad_col)) as.character(.data[[prioridad_col]]) else NA_character_,
    Rol = if (!is.na(rol_col)) as.character(.data[[rol_col]]) else NA_character_
  ) %>%
  separate_rows(GeneID, sep = ",|/|;") %>%
  mutate(
    GeneID = trimws(GeneID),
    GeneID = sub("\\..*$", "", GeneID)
  ) %>%
  filter(!is.na(GeneID), GeneID != "") %>%
  distinct(GeneID, .keep_all = TRUE)

curated2 <- curated2 %>%
  mutate(
    Texto_funcional = paste(BP, MF, KEGG, sep = " "),
    Tipo = case_when(
      !is.na(Tipo) & Tipo != "" ~ Tipo,
      str_detect(
        tolower(Texto_funcional),
        "phenylpropanoid|flavonoid|anthocyanin|ath00940|ath00941"
      ) ~ "Directa",
      str_detect(
        tolower(Texto_funcional),
        "secondary|jasmonic|salicylic|glutathione|suberin|cell wall|hormone|mapk|ath04075|ath00480|ath00073"
      ) ~ "Indirecta fuerte",
      TRUE ~ "Indirecta"
    ),
    Tipo = str_replace_all(Tipo, "🔴|🟠|🟡", ""),
    Tipo = trimws(Tipo),
    Tipo = case_when(
      str_detect(tolower(Tipo), "directa") & !str_detect(tolower(Tipo), "indirecta") ~ "Directa",
      str_detect(tolower(Tipo), "indirecta fuerte") ~ "Indirecta fuerte",
      str_detect(tolower(Tipo), "indirecta") ~ "Indirecta",
      TRUE ~ Tipo
    )
  )

cat("Genes curados detectados:", nrow(curated2), "\n")

# =====================================================
# LEER DESEQ2
# =====================================================

read_deseq <- function(file, contrast) {
  df <- read.csv(file, row.names = 1, check.names = FALSE)
  df$GeneID <- sub("\\..*$", "", rownames(df))
  df$contrast <- contrast
  
  df %>%
    mutate(
      padj_plot = ifelse(is.na(padj), 1, padj),
      padj_plot = ifelse(padj_plot == 0, .Machine$double.xmin, padj_plot),
      neglog10padj = -log10(padj_plot),
      regulation = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange < lfc_cutoff ~ "Downregulated",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange > abs(lfc_cutoff) ~ "Upregulated",
        TRUE ~ "Not significant"
      )
    )
}

all_data <- bind_rows(
  read_deseq(files[[1]], names(files)[1]),
  read_deseq(files[[2]], names(files)[2]),
  read_deseq(files[[3]], names(files)[3])
)

matched_genes <- intersect(curated2$GeneID, all_data$GeneID)

cat("Genes curados encontrados en DESeq2:", length(matched_genes), "\n")

if (length(matched_genes) == 0) {
  stop("ERROR: Ningún gen curado coincide con los GeneID de DESeq2.")
}

plot_data <- all_data %>%
  left_join(curated2, by = "GeneID") %>%
  mutate(
    Fenilpropanoide = ifelse(GeneID %in% curated2$GeneID, "Sí", "No"),
    Tipo_plot = ifelse(is.na(Tipo), "No fenilpropanoide", Tipo),
    Tipo_plot = factor(
      Tipo_plot,
      levels = c("Directa", "Indirecta fuerte", "Indirecta", "No fenilpropanoide")
    )
  )

# =====================================================
# GENES CANDIDATOS
# =====================================================

candidate_genes <- plot_data %>%
  filter(
    contrast == key_contrast,
    Fenilpropanoide == "Sí",
    !is.na(padj),
    padj < padj_cutoff,
    log2FoldChange < lfc_cutoff
  ) %>%
  arrange(log2FoldChange) %>%
  pull(GeneID) %>%
  unique()

cat("Genes candidatos downregulated en", key_contrast, ":", length(candidate_genes), "\n")
print(candidate_genes)

if (length(candidate_genes) == 0) {
  warning("No hay genes curados downregulated en 24h_vs_4h. Se usarán todos los genes curados encontrados.")
  candidate_genes <- matched_genes
}

# =====================================================
# PALETAS TIPO PUBLICACIÓN
# =====================================================

col_reg <- c(
  "Downregulated" = "#2C7BB6",
  "Upregulated" = "#D7191C",
  "Not significant" = "grey82"
)

col_type <- c(
  "Directa" = "#D7191C",
  "Indirecta fuerte" = "#FDAE61",
  "Indirecta" = "#636363",
  "No fenilpropanoide" = "grey70"
)

shape_type <- c(
  "Directa" = 21,
  "Indirecta fuerte" = 24,
  "Indirecta" = 22,
  "No fenilpropanoide" = 16
)

contrast_labels <- c(
  "4h_vs_control" = "4h vs control",
  "24h_vs_control" = "24h vs control",
  "24h_vs_4h" = "24h vs 4h"
)

# =====================================================
# 1. VOLCANO TRIPLE
# =====================================================

make_volcano <- function(df, title) {
  
  xlim_val <- max(abs(plot_data$log2FoldChange), na.rm = TRUE)
  xlim_val <- min(max(xlim_val, 4), 10)
  
  ggplot(df, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(
      aes(color = regulation),
      alpha = 0.45,
      size = 1.1
    ) +
    geom_point(
      data = df %>% filter(GeneID %in% candidate_genes),
      aes(fill = Tipo_plot, shape = Tipo_plot),
      color = "black",
      size = 3.2,
      stroke = 0.7
    ) +
    geom_text_repel(
      data = df %>%
        filter(GeneID %in% candidate_genes) %>%
        filter(
          contrast == key_contrast |
            (!is.na(padj) & padj < 0.01 & abs(log2FoldChange) > 1.5)
        ),
      aes(label = GeneID),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.size = 0.25
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.4) +
    scale_color_manual(values = col_reg) +
    scale_fill_manual(values = col_type, drop = FALSE) +
    scale_shape_manual(values = shape_type, drop = FALSE) +
    coord_cartesian(xlim = c(-xlim_val, xlim_val)) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = "Regulation",
      fill = "Phenylpropanoid relation",
      shape = "Phenylpropanoid relation"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

p1 <- make_volcano(plot_data %>% filter(contrast == "4h_vs_control"), "A. 4h vs control")
p2 <- make_volcano(plot_data %>% filter(contrast == "24h_vs_control"), "B. 24h vs control")
p3 <- make_volcano(plot_data %>% filter(contrast == "24h_vs_4h"), "C. 24h vs 4h")

volcano_triple <- (p1 + p2 + p3) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(outdir, "Figure1_volcano_triple_publication.pdf"),
  volcano_triple,
  width = 16,
  height = 5.8
)

ggsave(
  file.path(outdir, "Figure1_volcano_triple_publication.png"),
  volcano_triple,
  width = 16,
  height = 5.8,
  dpi = 300
)

# =====================================================
# 2. HEATMAP
# =====================================================

heatmap_df <- plot_data %>%
  filter(GeneID %in% candidate_genes) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("4h_vs_control", "24h_vs_control", "24h_vs_4h")
    ),
    contrast_label = factor(
      contrast_labels[as.character(contrast)],
      levels = c("4h vs control", "24h vs control", "24h vs 4h")
    ),
    Tipo_plot = factor(Tipo_plot, levels = c("Directa", "Indirecta fuerte", "Indirecta"))
  )

gene_order <- heatmap_df %>%
  filter(contrast == key_contrast) %>%
  arrange(Tipo_plot, log2FoldChange) %>%
  pull(GeneID) %>%
  unique()

heatmap_df <- heatmap_df %>%
  mutate(GeneID = factor(GeneID, levels = rev(gene_order)))

heatmap_plot <- ggplot(
  heatmap_df,
  aes(x = contrast_label, y = GeneID, fill = log2FoldChange)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(Tipo_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "log2FC"
  ) +
  labs(
    title = "Candidate downregulated phenylpropanoid-related genes",
    subtitle = "Genes selected by padj < 0.05 and log2FC < -1 in 24h vs 4h",
    x = "",
    y = "Gene ID"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text.y.left = element_text(face = "bold", angle = 90),
    axis.title = element_text(face = "bold"),
    plot.margin = margin(10, 30, 10, 30)
  )

ggsave(
  file.path(outdir, "Figure2_heatmap_candidates_publication.pdf"),
  heatmap_plot,
  width = 8.8,
  height = 6.8
)

ggsave(
  file.path(outdir, "Figure2_heatmap_candidates_publication.png"),
  heatmap_plot,
  width = 8.8,
  height = 6.8,
  dpi = 300
)

# =====================================================
# 3. INTEGRATIVE PLOT
# =====================================================

integrative_df <- plot_data %>%
  filter(GeneID %in% candidate_genes) %>%
  mutate(
    contrast = factor(
      contrast,
      levels = c("4h_vs_control", "24h_vs_control", "24h_vs_4h")
    ),
    contrast_label = factor(
      contrast_labels[as.character(contrast)],
      levels = c("4h vs control", "24h vs control", "24h vs 4h")
    ),
    Tipo_plot = factor(Tipo_plot, levels = c("Directa", "Indirecta fuerte", "Indirecta")),
    GeneID = factor(GeneID, levels = rev(gene_order))
  )

integrative_plot <- ggplot(
  integrative_df,
  aes(x = contrast_label, y = GeneID)
) +
  geom_point(
    aes(
      size = abs(log2FoldChange),
      fill = log2FoldChange,
      shape = Tipo_plot
    ),
    color = "black",
    stroke = 0.7
  ) +
  facet_grid(Tipo_plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(
    low = "#2C7BB6",
    mid = "white",
    high = "#D7191C",
    midpoint = 0,
    limits = c(-4, 4),
    oob = squish,
    name = "log2FC"
  ) +
  scale_shape_manual(
    values = shape_type,
    drop = FALSE,
    name = "Functional relation"
  ) +
  scale_size_continuous(
    name = "|log2FC|",
    range = c(2.5, 8)
  ) +
  labs(
    title = "Integrative expression profile of phenylpropanoid-related candidates",
    subtitle = "Bubble size represents |log2FC|; color represents direction and magnitude",
    x = "",
    y = "Gene ID"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text.y.left = element_text(face = "bold", angle = 90),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(10, 60, 10, 30)
  )

ggsave(
  file.path(outdir, "Figure3_integrative_candidates_publication.pdf"),
  integrative_plot,
  width = 10.5,
  height = 6.8
)

ggsave(
  file.path(outdir, "Figure3_integrative_candidates_publication.png"),
  integrative_plot,
  width = 10.5,
  height = 6.8,
  dpi = 300
)

# =====================================================
# 4. TABLAS DE SOPORTE
# =====================================================

candidate_table <- plot_data %>%
  filter(GeneID %in% candidate_genes) %>%
  select(
    GeneID,
    contrast,
    log2FoldChange,
    padj,
    regulation,
    Tipo_plot,
    BP,
    MF,
    KEGG,
    Rol,
    Prioridad
  ) %>%
  arrange(GeneID, contrast)

write.csv(
  candidate_table,
  file.path(outdir, "candidate_genes_for_publication_figures.csv"),
  row.names = FALSE
)

summary_table <- plot_data %>%
  group_by(contrast, regulation, Tipo_plot) %>%
  summarise(n = n(), .groups = "drop")

write.csv(
  summary_table,
  file.path(outdir, "summary_regulation_by_category.csv"),
  row.names = FALSE
)

write.csv(
  curated2,
  file.path(outdir, "curated_genes_expanded.csv"),
  row.names = FALSE
)

cat("\n========================================\n")
cat("FIGURAS OPTIMIZADAS GENERADAS\n")
cat("========================================\n")
cat("Salida:", outdir, "\n")
cat("Genes candidatos usados:", length(candidate_genes), "\n")
cat("- Figure1_volcano_triple_publication.pdf/png\n")
cat("- Figure2_heatmap_candidates_publication.pdf/png\n")
cat("- Figure3_integrative_candidates_publication.pdf/png\n")
cat("- candidate_genes_for_publication_figures.csv\n")
cat("- summary_regulation_by_category.csv\n")
