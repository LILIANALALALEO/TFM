#!/bin/bash
set -euo pipefail

echo "========================================"
echo "        QC WGBS RAW POR ESTRÉS"
echo "========================================"

read -p "Estrés a procesar: salt o cold [salt]: " TARGET_STRESS
TARGET_STRESS=${TARGET_STRESS:-salt}

if [[ "$TARGET_STRESS" != "salt" && "$TARGET_STRESS" != "cold" ]]; then
  echo "Error: el estrés debe ser 'salt' o 'cold'"
  exit 1
fi

THREADS=4

BASE_IN="Data/Raw/WGBS"
BASE_OUT="Data/Processed/06.WGBS_QC"

INDIR="${BASE_IN}/${TARGET_STRESS}"
OUTDIR="${BASE_OUT}/${TARGET_STRESS}"

FASTQC_OUT="${OUTDIR}/fastqc_raw"
MULTIQC_OUT="${OUTDIR}/multiqc_raw"

mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

echo
echo "Parámetros:"
echo "  Estrés      = $TARGET_STRESS"
echo "  Entrada     = $INDIR"
echo "  FastQC out  = $FASTQC_OUT"
echo "  MultiQC out = $MULTIQC_OUT"
echo

command -v fastqc >/dev/null || { echo "Error: fastqc no encontrado"; exit 1; }
command -v multiqc >/dev/null || { echo "Error: multiqc no encontrado"; exit 1; }

NFILES=$(find "$INDIR" -name "*.fastq.gz" -type f | wc -l)

if [[ "$NFILES" -eq 0 ]]; then
  echo "Error: no se encontraron archivos .fastq.gz en:"
  echo "$INDIR"
  exit 1
fi

echo "Archivos FASTQ encontrados: $NFILES"
echo

read -p "¿Desea continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  echo "Ejecución cancelada."
  exit 0
fi

echo
echo "Ejecutando FastQC..."
echo

find "$INDIR" -name "*.fastq.gz" -type f | while read file
do
  echo "Procesando: $file"

  fastqc "$file" \
    --threads "$THREADS" \
    --outdir "$FASTQC_OUT"
done

echo
echo "Generando MultiQC..."
multiqc "$FASTQC_OUT" \
  -o "$MULTIQC_OUT" \
  --force

echo
echo "========================================"
echo " QC WGBS RAW COMPLETADO"
echo " Estrés procesado: $TARGET_STRESS"
echo "========================================"
echo "Reporte MultiQC:"
echo "${MULTIQC_OUT}/multiqc_report.html"
