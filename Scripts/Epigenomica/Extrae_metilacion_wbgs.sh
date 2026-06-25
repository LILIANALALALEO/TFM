#!/bin/bash
set -euo pipefail

echo "========================================"
echo " 10. LIGHT METHYLATION EXTRACTION"
echo "========================================"

ulimit -n 4096 || true

read -p "Estrés a procesar [salt]: " STRESS
STRESS=${STRESS:-salt}

INDIR="Data/Processed/09.Bismark_Deduplication/${STRESS}"
OUTBASE="Data/Processed/10.Methylation_Calls/${STRESS}"

mkdir -p "$OUTBASE"

echo
echo "Input : $INDIR"
echo "Output: $OUTBASE"
echo

NBAM=$(find "$INDIR" -name "*deduplicated.bam" | wc -l)

if [[ "$NBAM" -eq 0 ]]; then
  echo "No se encontraron BAM deduplicados"
  exit 1
fi

echo "BAM encontrados: $NBAM"
echo

read -p "¿Desea continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  exit 0
fi

find "$INDIR" -name "*deduplicated.bam" | sort | while read BAM
do

  CONDITION=$(basename "$(dirname "$BAM")")

  OUTDIR="${OUTBASE}/${CONDITION}"

  mkdir -p "$OUTDIR"

  echo "----------------------------------------"
  echo "Procesando:"
  echo "$BAM"
  echo "----------------------------------------"

  df -h .

  bismark_methylation_extractor \
    --paired-end \
    --no_overlap \
    --gzip \
    --bedGraph \
    --counts \
    --buffer_size 2G \
    --output "$OUTDIR" \
    "$BAM"

  echo
  echo "Extracción completada"
  echo

done

echo
echo "========================================"
echo " METHYLATION EXTRACTION FINALIZADO"
echo "========================================"