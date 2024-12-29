### README for Differential Expression and Survival Analysis Pipeline

---

## Overview
This project contains an R-based pipeline for processing and analyzing clinical and multi-omics data (RNA-Seq, CNA) to study differential gene expression, gene ontology, pathway enrichment, and survival analysis for breast cancer data.

The primary dataset used is from TCGA's BRCA cohort (`brca_tcga_pan_can_atlas_2018`).

---

## Workflow Description
### 1. **Setup and File Management**
- Set the working directory.
- Decompress the dataset archive (`brca_tcga_pan_can_atlas_2018.tar.gz`).
- Load clinical (`data_clinical_patient.txt`), RNA-Seq (`data_mrna_seq_v2_rsem.txt`), and CNA (`data_cna.txt`) data into R.

### 2. **Data Preprocessing**
- **Clinical Data**: Filter and clean unnecessary rows, retaining relevant patient identifiers.
- **CNA Data**: Standardize patient IDs for consistency with clinical and RNA-Seq data.
- **RNA-Seq Data**: Prepare the matrix for downstream differential expression analysis.

### 3. **HER2 Status and ERBB2 Analysis**
- Identify HER2 status using ERBB2 counts from the CNA dataset.
- Match RNA-Seq patient IDs with CNA and clinical datasets.

### 4. **Differential Expression Analysis**
- Perform normalization using the `DESeq2` package.
- Filter genes with low counts and apply variance stabilization.
- Conduct differential expression analysis, identifying genes with significant adjusted p-values (`padj < 0.05`).
- Visualize top differentially expressed genes using PCA and heatmaps.

### 5. **Functional Enrichment Analysis**
- **Gene Ontology (GO)**: Analyze over- and under-expressed genes for biological processes.
- **Pathway Enrichment**: Perform KEGG and Reactome pathway enrichment for significant genes.
- Visualize results with dot plots and tree plots.

### 6. **Survival Analysis**
- Map clinical and RNA-Seq data to identify matching patient IDs.
- Prepare survival data (`time` and `status`) and perform Cox proportional hazard regression using elastic net regularization.
- Compute risk scores and categorize patients into "High Risk" and "Low Risk" groups.
- Visualize survival probabilities using Kaplan-Meier plots.

---

## Prerequisites
### Required Software
- R (v4.0 or later)

### Required R Packages
The following R packages must be installed:
- **Core Packages**: `DESeq2`, `clusterProfiler`, `ReactomePA`, `enrichplot`, `pathview`
- **Survival Analysis**: `survival`, `glmnet`, `survminer`
- **Data Manipulation and Visualization**: `pheatmap`, `ggplot2`, `xlsx`

Use the following commands to install missing packages:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "ReactomePA", "org.Hs.eg.db"))
install.packages(c("survminer", "glmnet", "pheatmap", "ggplot2", "xlsx"))
```

---

## Usage
1. Clone or download this repository.
2. Set the working directory to the folder containing the dataset:
   ```R
   setwd("C:/Users/KIIT/Desktop/bio-project/")
   ```
3. Run the script step-by-step in an R environment to process the data, analyze differential expression, and perform survival analysis.

---

## Outputs
- **Differential Expression Results**: Gene names, log2 fold changes, and adjusted p-values (`DEG.xlsx`).
- **Enrichment Analysis Plots**: Dot plots for GO, KEGG, and Reactome pathways.
- **Survival Analysis**: Kaplan-Meier survival curves for risk groups.
- **Heatmaps**: Top 20 differentially expressed genes.

