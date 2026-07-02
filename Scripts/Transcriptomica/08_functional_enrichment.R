#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(AnnotationDbi)
})

# -----------------------------------
# Directorios
# -----------------------------------
indir <- "Data/Processed/04.DESeq2"
outdir <- "Data/Processed/05.Functional_Enrichment"
results_dir <- "Results"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------
# Buscar automáticamente resultados génicos DESeq2
# -----------------------------------
result_files <- list.files(
  path = indir,
  pattern = "^deseq2_gene_level_.*\\.csv$",
  full.names = TRUE
)

if (length(result_files) == 0) {
  stop("ERROR: No se encontraron archivos 'deseq2_gene_level_*.csv' en ", indir)
}

cat("Archivos detectados para enriquecimiento:\n")
print(basename(result_files))
cat("\n")

# -----------------------------------
# Función para procesar cada contraste
# -----------------------------------
process_contrast <- function(file_path) {

  contrast_name <- basename(file_path)
  contrast_name <- sub("^deseq2_gene_level_", "", contrast_name)
  contrast_name <- sub("\\.csv$", "", contrast_name)

  contrast_outdir <- file.path(outdir, contrast_name)
  dir.create(contrast_outdir, showWarnings = FALSE, recursive = TRUE)

  cat("========================================\n")
  cat("Procesando contraste:", contrast_name, "\n")
  cat("========================================\n")

  # -----------------------------
  # Leer resultados DESeq2
  # -----------------------------
  res_gene <- read.csv(
    file_path,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Convertir rownames en columna gene_id
  res_gene$gene_id <- rownames(res_gene)

  required_cols <- c("gene_id", "log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(res_gene))) {
    warning("Archivo omitido por columnas faltantes: ", basename(file_path))
    return(NULL)
  }

  res_gene <- res_gene %>%
    filter(!is.na(gene_id), !is.na(log2FoldChange), !is.na(padj))

  if (nrow(res_gene) == 0) {
    warning("Archivo sin filas útiles: ", basename(file_path))
    return(NULL)
  }

  # -----------------------------
  # Definir conjuntos de genes
  # -----------------------------
  down_genes_df <- res_gene %>%
    filter(padj < 0.05, log2FoldChange < -1) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(padj)

  up_genes_df <- res_gene %>%
    filter(padj < 0.05, log2FoldChange > 1) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    arrange(padj)

  down_genes <- unique(down_genes_df$gene_id)
  up_genes <- unique(up_genes_df$gene_id)
  all_genes <- unique(res_gene$gene_id)

  # -----------------------------
  # Guardar listas
  # -----------------------------
  write.table(
    down_genes,
    file.path(contrast_outdir, "downregulated_genes.txt"),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  write.table(
    up_genes,
    file.path(contrast_outdir, "upregulated_genes.txt"),
    row.names = FALSE, col.names = FALSE, quote = FALSE
  )

  write.csv(
    down_genes_df,
    file.path(contrast_outdir, "downregulated_genes_table.csv"),
    row.names = FALSE
  )

  write.csv(
    up_genes_df,
    file.path(contrast_outdir, "upregulated_genes_table.csv"),
    row.names = FALSE
  )

  cat("Genes downregulated:", length(down_genes), "\n")
  cat("Genes upregulated:", length(up_genes), "\n")
  cat("Genes del universo:", length(all_genes), "\n\n")

  if (length(down_genes) == 0) {
    writeLines(
      c(
        paste("Contraste:", contrast_name),
        "No se detectaron genes downregulated con padj < 0.05 y log2FC < -1."
      ),
      con = file.path(contrast_outdir, "enrichment_summary.txt")
    )
    return(NULL)
  }

  # -----------------------------
  # GO enrichment
  # -----------------------------
  ego_bp <- enrichGO(
    gene = down_genes,
    universe = all_genes,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5,
    maxGSSize = 500
  )

  ego_mf <- enrichGO(
    gene = down_genes,
    universe = all_genes,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5,
    maxGSSize = 500
  )

  ego_cc <- enrichGO(
    gene = down_genes,
    universe = all_genes,
    OrgDb = org.At.tair.db,
    keyType = "TAIR",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 5,
    maxGSSize = 500
  )

  go_bp_df <- as.data.frame(ego_bp)
  go_mf_df <- as.data.frame(ego_mf)
  go_cc_df <- as.data.frame(ego_cc)

  write.csv(go_bp_df, file.path(contrast_outdir, "downregulated_GO_BP.csv"), row.names = FALSE)
  write.csv(go_mf_df, file.path(contrast_outdir, "downregulated_GO_MF.csv"), row.names = FALSE)
  write.csv(go_cc_df, file.path(contrast_outdir, "downregulated_GO_CC.csv"), row.names = FALSE)

  # -----------------------------
  # KEGG enrichment
  # -----------------------------
  mapping <- AnnotationDbi::select(
    org.At.tair.db,
    keys = all_genes,
    keytype = "TAIR",
    columns = c("TAIR", "ENTREZID", "SYMBOL", "GENENAME")
  ) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct()

  down_entrez <- unique(mapping$ENTREZID[mapping$TAIR %in% down_genes])
  all_entrez <- unique(mapping$ENTREZID[mapping$TAIR %in% all_genes])

  if (length(down_entrez) > 0 && length(all_entrez) > 0) {
    ekegg <- enrichKEGG(
      gene = down_entrez,
      universe = all_entrez,
      organism = "ath",
      keyType = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      minGSSize = 5,
      maxGSSize = 500
    )
    kegg_df <- as.data.frame(ekegg)
  } else {
    kegg_df <- data.frame()
  }

  write.csv(
    kegg_df,
    file.path(contrast_outdir, "downregulated_KEGG.csv"),
    row.names = FALSE
  )

  # -----------------------------
  # Resumen
  # -----------------------------
  summary_lines <- c(
    paste("Contraste:", contrast_name),
    paste("Genes downregulated:", length(down_genes)),
    paste("Genes upregulated:", length(up_genes)),
    paste("Genes del universo:", length(all_genes)),
    paste("GO BP enriquecidos:", nrow(go_bp_df)),
    paste("GO MF enriquecidos:", nrow(go_mf_df)),
    paste("GO CC enriquecidos:", nrow(go_cc_df)),
    paste("KEGG enriquecidos:", nrow(kegg_df))
  )

  writeLines(
    summary_lines,
    con = file.path(contrast_outdir, "enrichment_summary.txt")
  )

  writeLines(
    summary_lines,
    con = file.path(results_dir, paste0("summary_", contrast_name, ".txt"))
  )

  cat("Enriquecimiento completado para:", contrast_name, "\n\n")
}

# -----------------------------------
# Ejecutar en todos los contrastes
# -----------------------------------
for (f in result_files) {
  process_contrast(f)
}

cat("Procesamiento de GO/KEGG completado correctamente.\n")
