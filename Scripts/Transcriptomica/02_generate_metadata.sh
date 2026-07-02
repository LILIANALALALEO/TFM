#!/bin/bash
set -euo pipefail

INPUT="Code/sra_metadata.csv"
OUT="Data/Processed/04.DESeq2/sample_metadata.csv"

mkdir -p Data/Processed/04.DESeq2

echo "sample,condition" > "$OUT"

tail -n +2 "$INPUT" | while IFS=',' read -r SRA CONDITION STRESS
do
  echo "${SRA},${CONDITION}" >> "$OUT"
done

echo "Metadata generado automáticamente"
