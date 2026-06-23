#!/bin/bash
set -euo pipefail

echo "========================================"
echo "        TRIM GALORE AUTOMATIZADO        "
echo "========================================"

read -p "Archivo metadata [Code/sra_metadata_cold.csv]: " METADATA
METADATA=${METADATA:-Code/sra_metadata_cold.csv}

if [[ ! -f "$METADATA" ]]; then
  echo "Error: no existe el archivo $METADATA"
  exit 1
fi

read -p "Número de hilos [6]: " THREADS
THREADS=${THREADS:-6}

read -p "Corte en 5' para R1 y R2 [15]: " CLIP5
CLIP5=${CLIP5:-15}

read -p "Corte en 3' para R1 y R2 [5]: " CLIP3
CLIP3=${CLIP3:-5}

read -p "Calidad mínima Phred [25]: " QUALITY
QUALITY=${QUALITY:-25}

read -p "Longitud mínima de lectura [129]: " MINLEN
MINLEN=${MINLEN:-129}

read -p "Carpeta de salida [Data/Processed/02.Trimming]: " OUTDIR
OUTDIR=${OUTDIR:-Data/Processed/02.Trimming}

FASTQC_OUT="${OUTDIR}/fastqc_trimmed"
mkdir -p "$OUTDIR" "$FASTQC_OUT"

echo
echo "Parámetros seleccionados:"
echo "  METADATA  = $METADATA"
echo "  THREADS   = $THREADS"
echo "  CLIP5     = $CLIP5"
echo "  CLIP3     = $CLIP3"
echo "  QUALITY   = $QUALITY"
echo "  MINLEN    = $MINLEN"
echo "  OUTDIR    = $OUTDIR"
echo "  FASTQC    = $FASTQC_OUT"
echo

read -p "¿Desea continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  echo "Ejecución cancelada."
  exit 0
fi

echo
echo "Iniciando trimming..."
echo

tail -n +2 "$METADATA" | while IFS=',' read -r SRA CONDITION STRESS
do
  r1="Data/Raw/${STRESS}/${CONDITION}/${SRA}_1.fastq.gz"
  r2="Data/Raw/${STRESS}/${CONDITION}/${SRA}_2.fastq.gz"

  echo "----------------------------------------"
  echo "Procesando muestra: $SRA"
  echo "Condición: $CONDITION"
  echo "Estrés: $STRESS"
  echo "R1: $r1"
  echo "R2: $r2"

  if [[ ! -f "$r1" ]]; then
    echo "No se encontró R1 para $SRA. Se omite."
    continue
  fi

  if [[ ! -f "$r2" ]]; then
    echo "No se encontró R2 para $SRA. Se omite."
    continue
  fi

  trim_galore --paired \
    --clip_R1 "$CLIP5" \
    --clip_R2 "$CLIP5" \
    --three_prime_clip_R1 "$CLIP3" \
    --three_prime_clip_R2 "$CLIP3" \
    --quality "$QUALITY" \
    --length "$MINLEN" \
    --fastqc \
    --fastqc_args "--outdir $FASTQC_OUT" \
    --cores "$THREADS" \
    -o "$OUTDIR" \
    "$r1" "$r2"

done

echo
echo "Generando reporte MultiQC..."
multiqc "$OUTDIR" -o "$OUTDIR/multiqc"

echo
echo "Trimming finalizado correctamente."
echo "Resultados en: $OUTDIR"
