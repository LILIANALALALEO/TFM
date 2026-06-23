cd ~/TFM

mkdir -p Code
mkdir -p Data/Raw
mkdir -p Data/Processed/01.Quality_control
mkdir -p Data/Processed/02.Trimming
mkdir -p Data/Processed/03.Salmon
mkdir -p Data/Processed/04.DESeq2
mkdir -p Data/Processed/05.Functional_Enrichment
mkdir -p Data/Reference_genome
mkdir -p Data/Annotation
mkdir -p Results


#!/bin/bash
set -euo pipefail

THREADS=6
RAW_DIR="Data/Raw"

echo "========================================"
echo "   DESCARGA DE DATOS SRA (GENÉRICO)"
echo "========================================"

# 🔹 Preguntar metadata al usuario
read -p "Ingrese el archivo metadata (ej: Code/sra_metadata_cold.csv): " METADATA

# 🔹 Validar que exista
if [[ ! -f "$METADATA" ]]; then
    echo "❌ Error: El archivo $METADATA no existe"
    exit 1
fi

echo "📄 Usando metadata: $METADATA"
echo ""

mkdir -p "$RAW_DIR"

echo "📥 Iniciando descarga..."

tail -n +2 "$METADATA" | while IFS=',' read -r SRA CONDITION STRESS
do
    OUTDIR="${RAW_DIR}/${STRESS}/${CONDITION}"
    mkdir -p "$OUTDIR"

    echo "----------------------------------------"
    echo "🔹 Muestra: $SRA"
    echo "   Condición: $CONDITION"
    echo "   Estrés: $STRESS"
    echo "   Output: $OUTDIR"

    # Descargar SRA
    prefetch "$SRA"

    # Convertir a FASTQ
    fasterq-dump "$SRA" \
        --split-files \
        --threads "$THREADS" \
        -O "$OUTDIR"

    # Comprimir
    gzip -f "$OUTDIR/${SRA}_1.fastq"
    gzip -f "$OUTDIR/${SRA}_2.fastq"

    echo "✅ $SRA completado"
done

echo "----------------------------------------"
echo "🎉 Descarga finalizada"
