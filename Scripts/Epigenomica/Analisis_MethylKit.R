#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methylKit)
})

cat("========================================\n")
cat(" methylKit analysis - SALT\n")
cat("========================================\n\n")

outdir <- "Data/Processed/11.methylKit/salt"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

control_file <- "Data/Processed/10.Methylation_Calls/salt/control/SRR26113672_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
nacl_file    <- "Data/Processed/10.Methylation_Calls/salt/NaCl/SRR26113669_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"

if (!file.exists(control_file)) stop("No existe archivo control")
if (!file.exists(nacl_file)) stop("No existe archivo NaCl")

cat("Leyendo control...\n")
meth_control <- methRead(
  location = control_file,
  sample.id = "control",
  assembly = "TAIR10",
  treatment = 0,
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  mincov = 5
)

cat("Leyendo NaCl...\n")
meth_nacl <- methRead(
  location = nacl_file,
  sample.id = "NaCl",
  assembly = "TAIR10",
  treatment = 1,
  context = "CpG",
  pipeline = "bismarkCoverage",
  header = FALSE,
  mincov = 5
)

cat("Construyendo methylRawList...\n")
meth <- new(
  "methylRawList",
  list(meth_control, meth_nacl),
  treatment = c(0, 1)
)

cat("Filtrando cobertura...\n")
meth_filtered <- filterByCoverage(
  meth,
  lo.count = 5,
  hi.perc = 99.9
)

cat("Uniendo muestras...\n")
meth_united <- unite(
  meth_filtered,
  destrand = FALSE
)

saveRDS(meth_united, file.path(outdir, "meth_united_CpG_salt.rds"))

write.table(
  getData(meth_united),
  file = file.path(outdir, "meth_united_CpG_salt.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Calculando metilación diferencial exploratoria...\n")
diff_meth <- calculateDiffMeth(
  meth_united,
  overdispersion = "none"
)

saveRDS(diff_meth, file.path(outdir, "diff_meth_CpG_salt.rds"))

write.table(
  getData(diff_meth),
  file = file.path(outdir, "diff_meth_CpG_salt_all.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Filtrando DMCs exploratorios...\n")
sig_diff <- getMethylDiff(
  diff_meth,
  difference = 25,
  qvalue = 0.01
)

write.table(
  getData(sig_diff),
  file = file.path(outdir, "DMCs_CpG_salt_q001_diff25.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

hyper <- getMethylDiff(
  diff_meth,
  difference = 25,
  qvalue = 0.01,
  type = "hyper"
)

hypo <- getMethylDiff(
  diff_meth,
  difference = 25,
  qvalue = 0.01,
  type = "hypo"
)

write.table(
  getData(hyper),
  file = file.path(outdir, "DMCs_CpG_salt_hyper_q001_diff25.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  getData(hypo),
  file = file.path(outdir, "DMCs_CpG_salt_hypo_q001_diff25.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat(" methylKit FINALIZADO\n")
cat("========================================\n")
cat("Resultados en:\n")
cat(outdir, "\n")