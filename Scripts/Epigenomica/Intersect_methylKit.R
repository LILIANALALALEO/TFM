#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

cat("========================================\n")
cat(" 12. Intersección DMCs CpG con genes\n")
cat("========================================\n\n")

gff_file <- "Data/Annotation/Arabidopsis_thaliana.TAIR10.53.gff3"
dmcs_file <- "Data/Processed/11.methylKit/salt/DMCs_CpG_salt_q01_diff10.tsv"
outdir <- "Data/Processed/12.Epigenetic_Intersection/salt"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Leyendo GFF3...\n")
gff <- import(gff_file)

genes <- gff[gff$type == "gene"]

gene_id <- mcols(genes)$ID
gene_name <- if ("Name" %in% colnames(mcols(genes))) mcols(genes)$Name else gene_id

mcols(genes)$gene_id <- gene_id
mcols(genes)$gene_name <- gene_name

cat("Genes detectados:", length(genes), "\n")

cat("Leyendo DMCs...\n")
dmcs <- read.table(
  dmcs_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("DMCs leídos:", nrow(dmcs), "\n")

if (nrow(dmcs) == 0) {
  stop("No hay DMCs en el archivo seleccionado.")
}

dmcs_gr <- GRanges(
  seqnames = dmcs$chr,
  ranges = IRanges(start = dmcs$start, end = dmcs$end),
  strand = "*"
)

mcols(dmcs_gr) <- dmcs

cat("Construyendo regiones génicas...\n")

gene_body <- genes

promoters_2kb <- promoters(
  genes,
  upstream = 2000,
  downstream = 200
)

upstream_2kb <- promoters(
  genes,
  upstream = 2000,
  downstream = 1
)

mcols(gene_body)$region <- "gene_body"
mcols(promoters_2kb)$region <- "promoter_2kb"
mcols(upstream_2kb)$region <- "upstream_2kb"

intersect_region <- function(dmcs_gr, regions_gr, region_name) {
  hits <- findOverlaps(dmcs_gr, regions_gr, ignore.strand = TRUE)

  if (length(hits) == 0) {
    return(data.frame())
  }

  d <- as.data.frame(mcols(dmcs_gr[queryHits(hits)]))
  g <- as.data.frame(mcols(regions_gr[subjectHits(hits)]))

  out <- cbind(
    d,
    gene_id = g$gene_id,
    gene_name = g$gene_name,
    region = region_name
  )

  return(out)
}

cat("Intersectando DMCs con promotores, upstream y gene body...\n")

res_promoter <- intersect_region(dmcs_gr, promoters_2kb, "promoter_2kb")
res_upstream <- intersect_region(dmcs_gr, upstream_2kb, "upstream_2kb")
res_gene_body <- intersect_region(dmcs_gr, gene_body, "gene_body")

res_all <- rbind(
  res_promoter,
  res_upstream,
  res_gene_body
)

write.table(
  res_all,
  file = file.path(outdir, "DMCs_CpG_salt_intersect_genomic_regions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  unique(res_all[, c("gene_id", "gene_name", "region")]),
  file = file.path(outdir, "genes_with_DMCs_CpG_salt.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n========================================\n")
cat(" Intersección finalizada\n")
cat("========================================\n")
cat("Resultados en:\n")
cat(outdir, "\n\n")
cat("Archivos:\n")
cat("- DMCs_CpG_salt_intersect_genomic_regions.tsv\n")
cat("- genes_with_DMCs_CpG_salt.tsv\n")