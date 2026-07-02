#!/bin/bash
set -euo pipefail

echo "========================================"
echo "     CREAR METADATA WGBS"
echo "========================================"

OUTMETA="Code/epigenomica/wgbs_runs_final.csv"

mkdir -p Code/epigenomica

cat > "$OUTMETA" <<EOF
run,stress,condition
SRR26113672,salt,control
SRR26113669,salt,NaCl
SRR13764470,cold,control
SRR13764471,cold,control
SRR13764472,cold,cold
SRR13764473,cold,cold
EOF

echo "Metadata creado en:"
echo "$OUTMETA"
echo
cat "$OUTMETA"
