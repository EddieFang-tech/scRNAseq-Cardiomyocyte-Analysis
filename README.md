# scRNAseq-Cardiomyocyte-Analysis

## Overview

This repository contains R scripts for single-cell RNA sequencing (scRNA-seq) analysis of human embryonic stem cell-derived cardiomyocytes (hESC-CMs) and hypoimmunogenic human induced pluripotent stem cell-derived cardiomyocytes (hypo-hiPSC-CMs).

## Samples

Four scRNA-seq samples were included in the analysis:

| Sample ID | Cell type |
|---|---|
| TN331-11 | hESC-CMs |
| TN331-12 | hESC-CMs |
| TN331-13 | Hypo-hiPSC-CMs |
| TN331-14 | Hypo-hiPSC-CMs |

## Analysis

Single-cell RNA sequencing data were analyzed using the Seurat R package (version 5).

The analysis includes:

- Quality control and cell filtering
- Data normalization
- Principal component analysis (PCA)
- Uniform Manifold Approximation and Projection (UMAP)
- Cell clustering
- Gene expression analysis
- Visualization of cardiac electrophysiology-related genes

## Code

The analysis script is provided in:

`scRNA_electrophysiology_analysis.R`

## Data availability

The accession numbers for the raw sequencing data will be provided following deposition in a public repository.
