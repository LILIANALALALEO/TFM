#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
})

cat("========================================\n")
cat(" methylKit filtros flexibles - SALT\n")
cat("========================================\n\n")

indir <- "Data/Processed/11.methylKit/salt"
diff_file <- file.path(indir, "diff_meth_CpG_salt.rds")

if (!file.exists(diff_file)) {
  stop("No existe diff_meth_CpG_salt.rds")
}

cat("Cargando diff_meth...\n")
diff_meth <- readRDS(diff_file)

cat("Generando filtros flexibles...\n")

# Filtro 1: diferencia >= 10%, qvalue <= 0.1
dmc_diff10_q01 <- getMethylDiff(
  diff_meth,
  difference = 10,
  qvalue = 0.1
)

write.table(
  getData(dmc_diff10_q01),
  file = file.path(indir, "DMCs_CpG_salt_q01_diff10.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Hipermetilados con diferencia >= 10%, qvalue <= 0.1
hyper_diff10_q01 <- getMethylDiff(
  diff_meth,
  difference = 10,
  qvalue = 0.1,
  type = "hyper"
)

write.table(
  getData(hyper_diff10_q01),
  file = file.path(indir, "DMCs_CpG_salt_hyper_q01_diff10.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Hipometilados con diferencia >= 10%, qvalue <= 0.1
hypo_diff10_q01 <- getMethylDiff(
  diff_meth,
  difference = 10,
  qvalue = 0.1,
  type = "hypo"
)

write.table(
  getData(hypo_diff10_q01),
  file = file.path(indir, "DMCs_CpG_salt_hypo_q01_diff10.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Filtro exploratorio adicional: diferencia >= 5%, qvalue <= 0.1
dmc_diff5_q01 <- getMethylDiff(
  diff_meth,
  difference = 5,
  qvalue = 0.1
)

write.table(
  getData(dmc_diff5_q01),
  file = file.path(indir, "DMCs_CpG_salt_q01_diff5.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat(" Filtros flexibles finalizados\n")
cat("========================================\n")
cat("Resultados guardados en:\n")
cat(indir, "\n")
