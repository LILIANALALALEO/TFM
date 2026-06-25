#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

cat("========================================\n")
cat(" DMCs por región genómica\n")
cat("========================================\n\n")

infile <- "Data/Processed/12.Epigenetic_Intersection/salt/DMCs_CpG_salt_intersect_genomic_regions.tsv"

outdir <- "Results/Epigenomics/salt"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cat("Leyendo datos...\n")

df <- read.table(
  infile,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

cat("Número total de intersecciones:\n")
print(nrow(df))

region_counts <- df %>%
  count(region)

print(region_counts)

p <- ggplot(region_counts,
            aes(x = region,
                y = n)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Distribución de DMCs por región genómica",
    x = "Región genómica",
    y = "Número de DMCs"
  ) +
  theme_bw(base_size = 14)

ggsave(
  filename = file.path(outdir, "DMCs_by_genomic_region.png"),
  plot = p,
  width = 7,
  height = 5,
  dpi = 300
)

cat("\nFigura guardada en:\n")
cat(file.path(outdir, "DMCs_by_genomic_region.png"), "\n")