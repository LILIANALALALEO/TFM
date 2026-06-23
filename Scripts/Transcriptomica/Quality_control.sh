#!/bin/bash

set -e
set -o pipefail

THREADS=6

echo "========================================"
echo " CONTROL DE CALIDAD DE LECTURAS CRUDAS "
echo "========================================"

# Preguntar metadata
read -p "Ingrese el archivo metadata (ej: Code/sra_metadata_cold.csv): " METADATA

# Validar existencia
if [[ ! -f "$METADATA" ]]; then
  echo "Error: el archivo $METADATA no existe"
  exit 1
fi

OUT_BASE="Data/Processed/01.Quality_control"

# Obtener los estreses únicos del metadata
STRESSES=$(tail -n +2 "$METADATA" | cut -d',' -f3 | sort | uniq)

for STRESS in $STRESSES
do
  OUT_FASTQC="${OUT_BASE}/${STRESS}/fastqc_raw"
  OUT_MULTI="${OUT_BASE}/${STRESS}/multiqc_raw"

  mkdir -p "$OUT_FASTQC"
  mkdir -p "$OUT_MULTI"

  echo "----------------------------------------"
  echo "Procesando estrés: $STRESS"
  echo "Buscando archivos FASTQ/FASTQ.GZ en Data/Raw/${STRESS}"

  find "Data/Raw/${STRESS}" \( -name "*.fastq" -o -name "*.fastq.gz" \) > Code/tmp_fastq_list.txt

  if [ ! -s Code/tmp_fastq_list.txt ]; then
    echo "No se encontraron FASTQ para el estrés: $STRESS"
    continue
  fi

  fastqc $(cat Code/tmp_fastq_list.txt) \
    -o "$OUT_FASTQC" \
    -t "$THREADS"

  multiqc "$OUT_FASTQC" -o "$OUT_MULTI"

  echo "QC completado correctamente para $STRESS"
done

rm -f Code/tmp_fastq_list.txt

echo "========================================"
echo "Proceso de QC finalizado"
