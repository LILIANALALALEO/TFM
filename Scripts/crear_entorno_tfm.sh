#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="at_epi_tfm"

echo "=== 1. Desactivar entorno actual si existe ==="
conda deactivate 2>/dev/null || true

echo "=== 2. Eliminar entorno previo con el mismo nombre si existe ==="
conda env remove -n "${ENV_NAME}" -y 2>/dev/null || true

echo "=== 3. Configurar canales de conda ==="
conda config --add channels defaults || true
conda config --add channels bioconda || true
conda config --add channels conda-forge || true
conda config --set channel_priority strict

echo "=== 4. Crear entorno limpio ==="
conda create -n "${ENV_NAME}" -y \
  -c conda-forge -c bioconda -c defaults \
  python=3.10 \
  r-base=4.3 \
  r-essentials \
  r-tidyverse \
  r-optparse \
  r-data.table \
  r-pheatmap \
  r-ggplot2 \
  r-readxl \
  r-writexl \
  bioconductor-deseq2 \
  bioconductor-tximport \
  bioconductor-biomart \
  bioconductor-clusterprofiler \
  bioconductor-org.at.tair.db \
  bioconductor-complexheatmap \
  bioconductor-genomicranges \
  bioconductor-rtracklayer \
  bioconductor-annotationdbi \
  bioconductor-biocparallel \
  bioconductor-methylkit \
  bioconductor-bsseq \
  salmon \
  sra-tools \
  fastqc \
  multiqc \
  trim-galore \
  cutadapt \
  samtools \
  bedtools \
  bismark \
  bowtie2 \
  pigz \
  wget \
  curl \
  unzip

echo "=== 5. Activar entorno ==="
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

echo "=== 6. Verificaciones rápidas ==="
echo "-- Conda env activo:"
echo "$CONDA_DEFAULT_ENV"

echo "-- Versiones clave:"
R --version | head -n 1
python --version
salmon --version | head -n 1
fastqc --version
samtools --version | head -n 1
bedtools --version
bismark --version | head -n 1 || true
prefetch --version || true
fasterq-dump --version || true

echo "=== 7. Comprobación en R ==="
Rscript -e 'library(DESeq2); library(tximport); library(clusterProfiler); library(org.At.tair.db); library(GenomicRanges); library(methylKit); cat("R packages OK\n")'

echo "=== 8. Entorno listo ==="
echo "Actívalo con:"
echo "conda activate ${ENV_NAME}"





