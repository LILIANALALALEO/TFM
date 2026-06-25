#!/bin/bash
set -euo pipefail

echo "========================================"
echo "     TRIM GALORE WGBS POR ESTRÉS"
echo "========================================"

read -p "Archivo metadata [Code/epigenomica/wgbs_runs_final.csv]: " METADATA
METADATA=${METADATA:-Code/epigenomica/wgbs_runs_final.csv}

if [[ ! -f "$METADATA" ]]; then
  echo "Error: no existe el archivo $METADATA"
  exit 1
fi

read -p "Estrés a procesar: salt o cold [cold]: " TARGET_STRESS
TARGET_STRESS=${TARGET_STRESS:-cold}

if [[ "$TARGET_STRESS" != "salt" && "$TARGET_STRESS" != "cold" ]]; then
  echo "Error: el estrés debe ser 'salt' o 'cold'"
  exit 1
fi

read -p "Número de hilos [4]: " THREADS
THREADS=${THREADS:-4}

read -p "Corte en 5' para R1 y R2 [10]: " CLIP5
CLIP5=${CLIP5:-10}

read -p "Corte en 3' para R1 y R2 [5]: " CLIP3
CLIP3=${CLIP3:-5}

read -p "Calidad mínima Phred [20]: " QUALITY
QUALITY=${QUALITY:-20}

read -p "Longitud mínima de lectura [30]: " MINLEN
MINLEN=${MINLEN:-30}

RAWBASE="Data/Raw/WGBS"
OUTBASE="Data/Processed/07.WGBS_Trimming"

DATA_OUT="${OUTBASE}/data/${TARGET_STRESS}"
FASTQC_OUT="${OUTBASE}/qc/fastqc_trimmed/${TARGET_STRESS}"
REPORT_TRIM="${OUTBASE}/reports/${TARGET_STRESS}/trimming"
REPORT_MULTIQC="${OUTBASE}/reports/${TARGET_STRESS}/multiqc"

mkdir -p "$DATA_OUT" "$FASTQC_OUT" "$REPORT_TRIM" "$REPORT_MULTIQC"

echo
echo "Parámetros seleccionados:"
echo "  METADATA      = $METADATA"
echo "  TARGET_STRESS = $TARGET_STRESS"
echo "  THREADS       = $THREADS"
echo "  CLIP5         = $CLIP5"
echo "  CLIP3         = $CLIP3"
echo "  QUALITY       = $QUALITY"
echo "  MINLEN        = $MINLEN"
echo "  RAWBASE       = $RAWBASE"
echo "  DATA_OUT      = $DATA_OUT"
echo "  FASTQC_OUT    = $FASTQC_OUT"
echo "  REPORT_TRIM   = $REPORT_TRIM"
echo "  REPORT_MULTIQC= $REPORT_MULTIQC"
echo

read -p "¿Desea continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  echo "Ejecución cancelada."
  exit 0
fi

echo
echo "Iniciando trimming WGBS para estrés: $TARGET_STRESS"
echo

tail -n +2 "$METADATA" | while IFS=',' read -r RUN STRESS CONDITION
do
  if [[ "$STRESS" != "$TARGET_STRESS" ]]; then
    continue
  fi

  r1="${RAWBASE}/${STRESS}/${CONDITION}/${RUN}_1.fastq.gz"
  r2="${RAWBASE}/${STRESS}/${CONDITION}/${RUN}_2.fastq.gz"

  OUTDIR="${DATA_OUT}/${CONDITION}"
  FASTQC_COND_OUT="${FASTQC_OUT}/${CONDITION}"

  mkdir -p "$OUTDIR" "$FASTQC_COND_OUT"

  echo "----------------------------------------"
  echo "Procesando muestra: $RUN"
  echo "Estrés: $STRESS"
  echo "Condición: $CONDITION"
  echo "R1: $r1"
  echo "R2: $r2"
  echo "Salida: $OUTDIR"

  if [[ ! -f "$r1" ]]; then
    echo "No se encontró R1 para $RUN. Se omite."
    continue
  fi

  if [[ ! -f "$r2" ]]; then
    echo "No se encontró R2 para $RUN. Se omite."
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
    --fastqc_args "--outdir $FASTQC_COND_OUT" \
    --cores "$THREADS" \
    -o "$OUTDIR" \
    "$r1" "$r2"

done

echo
echo "Moviendo reportes de trimming..."
find "$DATA_OUT" -name "*_trimming_report.txt" -exec mv {} "$REPORT_TRIM"/ \;

echo
echo "Generando MultiQC para $TARGET_STRESS..."
multiqc "$DATA_OUT" "$FASTQC_OUT" "$REPORT_TRIM" \
  -o "$REPORT_MULTIQC" \
  --force

echo
echo "========================================"
echo " TRIMMING WGBS FINALIZADO"
echo " Estrés procesado: $TARGET_STRESS"
echo "========================================"
echo "Datos trimmed:"
echo "$DATA_OUT"
echo
echo "Reporte MultiQC:"
echo "${REPORT_MULTIQC}/multiqc_report.html"