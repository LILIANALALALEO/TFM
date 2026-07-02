# TFM – Regulación epigenética del metabolismo fenilpropanoide en *Arabidopsis thaliana* bajo estrés abiótico

## Descripción

Este repositorio contiene el código desarrollado para el Trabajo Fin de Máster (TFM) del Máster Universitario en Bioinformática de la Universidad Internacional de Valencia (VIU).

El objetivo del estudio fue identificar genes del metabolismo fenilpropanoide potencialmente regulados por mecanismos epigenéticos en *Arabidopsis thaliana* sometida a estrés salino y estrés por frío, mediante la integración de datos transcriptómicos (RNA-seq) y epigenómicos (WGBS).

---

# Objetivos

- Analizar datos públicos de RNA-seq bajo estrés salino y estrés por frío.
- Identificar genes diferencialmente expresados mediante DESeq2.
- Detectar citosinas diferencialmente metiladas utilizando Bismark y methylKit.
- Integrar transcriptómica y epigenómica.
- Priorizar genes candidatos relacionados con el metabolismo fenilpropanoide.

---

# Organización del repositorio

```
TFM/
│
├── Bibliografia/
│
├── Datos/
│
├── RawData/
│
└── Scripts/
    ├── Transcriptomica/
    └── Epigenomica/
```

---

# Pipeline transcriptómico

Los scripts se ejecutan en el siguiente orden:

| Paso | Script | Descripción |
|------|--------|-------------|
|1|01_download_rnaseq.sh|Descarga de datos RNA-seq desde SRA|
|2|02_generate_metadata.sh|Generación de metadatos|
|3|03_quality_control.sh|Control de calidad mediante FastQC y MultiQC|
|4|04_trim_reads.sh|Eliminación de adaptadores y bases de baja calidad|
|5|05_build_salmon_index.sh|Construcción del índice transcriptómico|
|6|06_salmon_quantification.sh|Cuantificación mediante Salmon|
|7|07_differential_expression.R|Análisis diferencial con DESeq2|
|8|08_functional_enrichment.R|Enriquecimiento funcional GO y KEGG|
|9|09_generate_report.R|Generación de tablas y resultados|

---

# Pipeline epigenómico

Los scripts se ejecutan en el siguiente orden:

| Paso | Script | Descripción |
|------|--------|-------------|
|1|01_download_wgbs.sh|Descarga de datos WGBS desde SRA|
|2|02_generate_metadata.sh|Generación de metadatos|
|3|03_quality_control.sh|Control de calidad|
|4|04_trim_reads.sh|Recorte de adaptadores|
|5|05_prepare_bismark_genome.sh|Preparación del genoma para Bismark|
|6|06_bismark_alignment.sh|Alineamiento con Bismark|
|7|07_remove_duplicates.sh|Eliminación de duplicados|
|8|08_extract_methylation.sh|Extracción de metilación|
|9|09_methylkit_filtering.R|Filtrado inicial de DMCs|
|10|10_methylkit_analysis.R|Análisis diferencial de metilación|
|11|11_intersect_genomic_features.R|Asignación de DMCs a regiones genómicas|
|12|12_generate_figures.R|Generación de figuras y resultados|

---

# Software utilizado

## Transcriptómica

- FastQC
- MultiQC
- Trim Galore
- Salmon
- tximport
- DESeq2
- clusterProfiler

## Epigenómica

- FastQC
- MultiQC
- Trim Galore
- Bismark
- Bowtie2
- methylKit
- GenomicRanges
- rtracklayer

---

# Datos

Los datos utilizados proceden de experimentos públicos disponibles en el repositorio Sequence Read Archive (SRA) del National Center for Biotechnology Information (NCBI).

---

# Reproducibilidad

Los scripts fueron desarrollados para ejecutarse en Linux utilizando Bash y R.

Las versiones de software empleadas se describen en la memoria del Trabajo Fin de Máster.

---

# Autor

Liliana Paulina Lalaleo Córdova

Máster Universitario en Bioinformática

Universidad Internacional de Valencia (VIU)

2026

---

# Licencia

Este repositorio tiene fines exclusivamente académicos y de investigación.
