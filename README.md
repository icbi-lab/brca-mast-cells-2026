# brca-mast-cells-2026

Analysis code accompanying
> [Kirchmair et al., Genes & Immunity (2026)](https://www.nature.com/gene/)

This repository contains the code to reproduce the analysis of transcriptional signatures of mast cells in public breast cancer cohorts, as described below.

### METABRIC
Data source: [cBioPortal](https://www.cbioportal.org/study/summary?id=brca_metabric)  
Analysis script: `analysis/METABRIC.Rmd`

### TCGA-BRCA
Data source: [GDC](https://portal.gdc.cancer.gov/projects/TCGA-BRCA)  
Analysis script: `analysis/TCGA-BRCA.Rmd`

### Single-cell atlas (Chen et al.)
Data source: [CELLxGENE](https://cellxgene.cziscience.com/collections/9432ae97-4803-4b9f-8f64-2b41e42ad3cb)  
Analysis script: `analysis/scRNAseq.Rmd`

Paper figures can be reproduced from the results by running `analysis/Figures.Rmd`.  
Analyses were performed in R 4.5.1.
