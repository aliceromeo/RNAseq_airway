# Airway RNA-Seq pipeline: automated transcriptomics workflow

This repository features an end-to-end bioinformatics pipeline for analyzing the ```airway``` RNA-Seq transcriptomics data, available through Bioconductor: [10.18129/B9.bioc.airway](https://doi.org/10.18129/B9.bioc.airway). The workflow is fully automated using **Snakemake** [(https://doi.org/10.12688/f1000research.29032.3)](https://doi.org/10.12688/f1000research.29032.3), ensuring a reproducible path from raw SRA data to functional enrichment analysis.

## 🧬 Workflow Overview
The pipeline automates the following biological and computational steps:
1. **Preprocessing**: Automatic download of raw FASTQ files (SRA) and reference genomic data (GENCODE v44).
2. **Quality Control**: Read quality assessment via **FastQC** and aggregated reporting with **MultiQC**.
3. **Quantification**: Transcriptome indexing and expression quantification using **Salmon**.
4. **Differential Expression Analysis (R/DESeq2)**: Data normalization, Exploratory Data Analysis and identification of DEGs.
5. **Functional Enrichment**: Biological interpretation through **Gene Ontology (GO)** and **Gene Set Enrichment Analysis (GSEA)**.

## 🛠️ Tech Stack
* **Workflow Manager**: Snakemake 
* **Languages**: R 
* **Bioinformatics tools**: Salmon, FastQC, MultiQC, parallel-fastq-dump 
* **Environment management**: Conda (two isolated environments for preprocessing and R analysis) 

## 📂 Project Structure
```text
├── Snakefile               # Core workflow logic
├── config/
│   └── config.yaml         # Execution parameters (thresholds, threads, paths)
├── envs/                   # Conda environment definitions (YAML)
│   ├── preproc_env.yaml
│   └── r_env.yaml
├── scripts/                # R scripts for DEA and enrichment
├── metadata/
│   └── metadata.csv        # Mapping of Sample Run IDs (SRA)
└── results/                # Generated tables and plots
```

## 🚀 Getting Started

### 1. Installation and Setup
Clone this repository and navigate into the project directory:
```bash
git clone https://github.com/aliceromeo/RNAseq_airway.git
cd RNAseq_airway
```

### 2. Execution
To run the full pipeline, automatically creating the necessary Conda environments and using (for example) 8 CPU cores:
```bash
snakemake --use-conda --cores 8
```

## 📈 Interactive analysis report
One of the key features of this workflow is the automatic generation of a comprehensive Snakemake Report. This HTML file encapsulates the entire analysis, providing:

1. Execution metrics: Time, date, and tools used.
2. A visual representation of the workflow's rule dependencies.
3. Embedded results: High-resolution plots (PCA, GSEA, Heatmaps) and data tables directly accessible within the browser.

To generate the report after a successful run, use:
```bash
snakemake --report report.html
```

You can then open ```report.html``` (already included in the repository as an example) in any web browser to explore the results.
