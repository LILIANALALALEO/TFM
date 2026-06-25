#!/bin/bash
set -euo pipefail

echo "========================================"
echo " PREPARACIÓN GENOMA BISMARK"
echo "========================================"

GENOME_DIR="Data/Reference_genome"

echo
echo "Usando genoma:"
echo "$GENOME_DIR"
echo

if [[ ! -f "$GENOME_DIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa" ]]; then
  echo "Error: no se encontró el genoma DNA"
  exit 1
fi

echo "Cabeceras detectadas:"
grep ">" "$GENOME_DIR/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa" | head

echo
read -p "¿Continuar? [s/n]: " CONFIRM
CONFIRM=${CONFIRM,,}

if [[ "$CONFIRM" != "s" && "$CONFIRM" != "si" ]]; then
  echo "Cancelado."
  exit 0
fi

echo
echo "Construyendo índice Bismark..."
echo

bismark_genome_preparation \
  --bowtie2 \
  "$GENOME_DIR"

echo
echo "========================================"
echo " ÍNDICE BISMARK COMPLETADO"
echo "========================================"