#!/bin/bash
set -euo pipefail

THREADS=6
INDEX="Data/Reference_genome/salmon_index"
OUTDIR="Data/Processed/03.Salmon"

mkdir -p "${OUTDIR}"

for r1 in Data/Processed/02.Trimming/*_1_val_1.fq; do
  sample=$(basename "${r1}" _1_val_1.fq)
  r2="Data/Processed/02.Trimming/${sample}_2_val_2.fq"

  if [[ ! -f "${r2}" ]]; then
    echo "Falta el par trimmeado de ${sample}; se omite."
    continue
  fi

  salmon quant \
    -i "${INDEX}" \
    -l A \
    -1 "${r1}" \
    -2 "${r2}" \
    -p "${THREADS}" \
    --validateMappings \
    -o "${OUTDIR}/${sample}"
done

echo "Cuantificación Salmon completada."
