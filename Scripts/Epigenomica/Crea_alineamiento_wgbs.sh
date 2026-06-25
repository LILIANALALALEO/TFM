#!/bin/bash
set -euo pipefail

echo "========================================"
echo "       08. BISMARK ALIGNMENT WGBS"
echo "========================================"

read -p "Estrés a procesar [salt]: " STRESS
STRESS=${STRESS:-salt}

THREADS=4

GENOME="$(pwd)/Data/Reference_genome/Bismark_TAIR10_genome"
INDIR="Data/Processed/07.WGBS_Trimming/data/${STRESS}"
OUTBASE="Data/Processed/08.Bismark_Alignment/${STRESS}"

mkdir -p "$OUTBASE"

echo
echo "Parámetros:"
echo "  STRESS = $STRESS"
echo "  GENOME = $GENOME"
echo "  INDIR  = $INDIR"
echo "  OUTDIR = $OUTBASE"
echo

if [[ ! -d "$GENOME/Bisulfite_Genome" ]]; then
  echo "Error: no existe el índice Bismark en:"
  echo "$GENOME/Bisulfite_Genome"
  exit 1
fi

NFILES=$(find "$INDIR" -name "*_1_val_1.fq.gz" -type f | wc -l)

if [[ "$NFILES" -eq 0 ]]; then
  echo "Error: no se encontraron archivos *_1_val_1.fq.gz en:"
  echo "$INDIR"
  exit 1
fi

echo "Muestras paired-end encontradas: $NFILES"
echo

read -p "¿Desea continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  echo "Ejecución cancelada."
  exit 0
fi

find "$INDIR" -name "*_1_val_1.fq.gz" -type f | sort | while read R1
do
  R2="${R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"
  SAMPLE=$(basename "$R1" _1_val_1.fq.gz)
  CONDITION=$(basename "$(dirname "$R1")")
  OUTDIR="${OUTBASE}/${CONDITION}"

  mkdir -p "$OUTDIR"

  echo "----------------------------------------"
  echo "Procesando muestra: $SAMPLE"
  echo "Condición: $CONDITION"
  echo "R1: $R1"
  echo "R2: $R2"
  echo "Salida: $OUTDIR"
  echo "----------------------------------------"

  if [[ ! -f "$R2" ]]; then
    echo "Error: no se encontró R2 para $SAMPLE"
    exit 1
  fi

  df -h .

  bismark \
    --genome "$GENOME" \
    -1 "$R1" \
    -2 "$R2" \
    --parallel 2 \
    -o "$OUTDIR"

done

echo
echo "========================================"
echo " ALINEAMIENTO BISMARK FINALIZADO"
echo "========================================"
echo "Resultados en:"
echo "$OUTBASE"