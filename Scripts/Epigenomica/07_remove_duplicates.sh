#!/bin/bash
set -euo pipefail

echo "========================================"
echo " 09. DEDUPLICACIÓN BISMARK WGBS"
echo "========================================"

read -p "Estrés a procesar [salt]: " STRESS
STRESS=${STRESS:-salt}

INDIR="Data/Processed/08.Bismark_Alignment/${STRESS}"
OUTBASE="Data/Processed/09.Bismark_Deduplication/${STRESS}"

mkdir -p "$OUTBASE"

echo
echo "Parámetros:"
echo "  STRESS  = $STRESS"
echo "  INDIR   = $INDIR"
echo "  OUTBASE = $OUTBASE"
echo

if [[ ! -d "$INDIR" ]]; then
  echo "Error: no existe la carpeta de alineamiento:"
  echo "$INDIR"
  exit 1
fi

NBAM=$(find "$INDIR" -name "*.bam" -type f | wc -l)

if [[ "$NBAM" -eq 0 ]]; then
  echo "Error: no se encontraron archivos BAM en:"
  echo "$INDIR"
  exit 1
fi

echo "BAM encontrados: $NBAM"
echo

read -p "¿Desea continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  echo "Ejecución cancelada."
  exit 0
fi

find "$INDIR" -name "*.bam" -type f | sort | while read BAM
do
  CONDITION=$(basename "$(dirname "$BAM")")
  OUTDIR="${OUTBASE}/${CONDITION}"

  mkdir -p "$OUTDIR"

  echo "----------------------------------------"
  echo "Procesando BAM:"
  echo "$BAM"
  echo "Condición: $CONDITION"
  echo "Salida: $OUTDIR"
  echo "----------------------------------------"

  df -h .

  deduplicate_bismark \
    --paired \
    --bam \
    --output_dir "$OUTDIR" \
    "$BAM"

done

echo
echo "========================================"
echo " DEDUPLICACIÓN FINALIZADA"
echo "========================================"
echo "Resultados en:"
echo "$OUTBASE"
