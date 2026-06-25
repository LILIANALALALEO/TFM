mkdir -p Code

mkdir -p Data/Raw/WGBS
mkdir -p Data/Raw/WGBS/salt/control
mkdir -p Data/Raw/WGBS/salt/NaCl
mkdir -p Data/Raw/WGBS/cold/control
mkdir -p Data/Raw/WGBS/cold/cold

mkdir -p Data/Processed/06.WGBS_QC
mkdir -p Data/Processed/07.WGBS_Trimming
mkdir -p Data/Processed/08.Bismark_Alignment
mkdir -p Data/Processed/09.Bismark_Deduplication
mkdir -p Data/Processed/10.Methylation_Calls
mkdir -p Data/Processed/11.methylKit
mkdir -p Data/Processed/12.Epigenetic_Intersection
mkdir -p Data/Reference_genome
mkdir -p Data/Annotation
mkdir -p Results/Epigenomics



#!/bin/bash
set -euo pipefail

echo "========================================"
echo " DESCARGA WGBS ORGANIZADA"
echo "========================================"

THREADS=4
TMPDIR="Data/Raw/WGBS/tmp_fasterq"
mkdir -p "$TMPDIR"

download_sample () {
  RUN=$1
  OUTDIR=$2

  mkdir -p "$OUTDIR"

  echo "----------------------------------------"
  echo "Descargando $RUN"
  echo "Salida: $OUTDIR"
  echo "----------------------------------------"

  fasterq-dump "$RUN" \
    --split-files \
    --threads "$THREADS" \
    --temp "$TMPDIR" \
    --outdir "$OUTDIR"

  echo "Comprimiendo $RUN..."
  gzip -f "$OUTDIR/${RUN}"*.fastq

  echo "Limpiando temporales..."
  rm -rf "$TMPDIR"/*
}

download_sample SRR26113672 Data/Raw/WGBS/salt/control
download_sample SRR26113673 Data/Raw/WGBS/salt/NaCl
download_sample SRR13764470 Data/Raw/WGBS/cold/control
download_sample SRR13764471 Data/Raw/WGBS/cold/control
download_sample SRR13764472 Data/Raw/WGBS/cold/cold

echo "========================================"
echo " DESCARGA COMPLETADA"
echo "========================================"