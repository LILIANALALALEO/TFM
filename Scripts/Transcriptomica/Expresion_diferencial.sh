#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(readr)
})

# -----------------------------
# Configuración
# -----------------------------
outdir <- "Data/Processed/04.DESeq2"
metadata_file <- file.path(outdir, "sample_metadata.csv")
salmon_dir <- "Data/Processed/03.Salmon"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Leer metadata
# -----------------------------
samples <- read.csv(metadata_file, stringsAsFactors = FALSE)

required_cols <- c("sample", "condition")
if (!all(required_cols %in% colnames(samples))) {
  stop("ERROR: sample_metadata.csv debe contener las columnas: sample y condition")
}

if (nrow(samples) == 0) {
  stop("ERROR: sample_metadata.csv está vacío")
}

samples$sample <- trimws(samples$sample)
samples$condition <- trimws(samples$condition)

# Mantener el orden del metadata
samples$condition <- factor(samples$condition, levels = unique(samples$condition))

cat("Metadata cargado correctamente:\n")
print(samples)
cat("\n")

# -----------------------------
# Rutas a quant.sf
# -----------------------------
files <- file.path(salmon_dir, samples$sample, "quant.sf")
names(files) <- samples$sample

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop(
    "ERROR: No se encontraron los siguientes archivos quant.sf:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("Todos los archivos quant.sf fueron encontrados.\n\n")

# -----------------------------
# Construir tx2gene
# -----------------------------
# Se usa el primer quant.sf para extraer los transcript IDs
first_quant <- files[1]

tx_ids <- read_tsv(
  first_quant,
  show_col_types = FALSE,
  progress = FALSE
)$Name

if (length(tx_ids) == 0) {
  stop("ERROR: No se pudieron leer transcript IDs desde ", first_quant)
}

tx2gene <- data.frame(
  transcript_id = tx_ids,
  gene_id = sub("\\..*$", "", tx_ids),
  stringsAsFactors = FALSE
)

# Validación mínima
if (any(tx2gene$gene_id == "" | is.na(tx2gene$gene_id))) {
  stop("ERROR: No se pudieron derivar gene_id correctamente desde los transcript IDs")
}

write.csv(
  tx2gene,
  file.path(outdir, "tx2gene_generated.csv"),
  row.names = FALSE
)

cat("Tabla tx2gene generada correctamente.\n")
cat("Ejemplo de tx2gene:\n")
print(head(tx2gene))
cat("\n")

# -----------------------------
# Importar Salmon a nivel génico
# countsFromAbundance='lengthScaledTPM' suele ser una buena opción
# para datos de Salmon + DESeq2
# -----------------------------
txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = FALSE
)

# -----------------------------
# Tabla de diseño
# -----------------------------
sampleTable <- data.frame(
  row.names = samples$sample,
  condition = samples$condition
)

cat("Número de muestras por condición:\n")
print(table(sampleTable$condition))
cat("\n")

# -----------------------------
# Crear objeto DESeq2
# -----------------------------
dds <- DESeqDataSetFromTximport(
  txi = txi,
  colData = sampleTable,
  design = ~ condition
)

# Filtrado básico
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

cat("Número de genes retenidos tras filtrado:", nrow(dds), "\n\n")

# -----------------------------
# Correr DESeq2
# -----------------------------
dds <- DESeq(dds, fitType = "local")

cat("Coeficientes disponibles:\n")
print(resultsNames(dds))
cat("\n")

# -----------------------------
# Guardar objeto principal y conteos normalizados
# -----------------------------
write.csv(
  as.data.frame(counts(dds, normalized = TRUE)),
  file.path(outdir, "normalized_counts_gene_level.csv")
)

saveRDS(dds, file.path(outdir, "dds_gene_level_object.rds"))

# -----------------------------
# Generar todos los contrastes posibles
# -----------------------------
conditions <- levels(samples$condition)

if (length(conditions) < 2) {
  stop("ERROR: Se requieren al menos dos condiciones para hacer contrastes")
}

contrast_pairs <- combn(conditions, 2, simplify = FALSE)

summary_file <- file.path(outdir, "contrast_summary_gene_level.txt")
cat("Resumen de contrastes DESeq2 a nivel de gen\n", file = summary_file)
cat("==========================================\n\n", file = summary_file, append = TRUE)

for (pair in contrast_pairs) {
  group1 <- pair[2]
  group2 <- pair[1]

  contrast_name <- paste0(group1, "_vs_", group2)
  cat("Procesando contraste:", contrast_name, "\n")

  res <- results(dds, contrast = c("condition", group1, group2))
  res <- res[order(res$pvalue), ]

  csv_file <- file.path(outdir, paste0("deseq2_gene_level_", contrast_name, ".csv"))
  rds_file <- file.path(outdir, paste0("res_gene_level_", contrast_name, ".rds"))

  write.csv(as.data.frame(res), csv_file)
  saveRDS(res, rds_file)

  sink(summary_file, append = TRUE)
  cat("Contraste:", contrast_name, "\n")
  print(summary(res))
  cat("\nPrimeras filas:\n")
  print(head(as.data.frame(res)))
  cat("\n----------------------------------------\n\n")
  sink()

  cat("  Guardado en:", csv_file, "\n")
}

cat("\nAnálisis DESeq2 a nivel de gen completado correctamente.\n")
cat("Resultados guardados en:", outdir, "\n")
