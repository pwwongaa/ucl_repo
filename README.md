# 🧠 UCL Interview – Transcriptomics Analysis of Alzheimer's Disease

## 🧪 Project Overview

A personal practice draft project on bulk/pseudobulk transcriptomics analysis, focusing on differential gene expression in **astrocytes** from the human **entorhinal cortex** using public **snRNA-seq** data from **Alzheimer's disease (AD)** and control samples.

---

## 📂 Dataset

**Source:** [NCBI GEO – GSE138852](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852)  
- Human entorhinal cortex  
- Single-nucleus RNA-seq (snRNA-seq)  
- Raw counts aggregated into pseudobulk format for AD vs. control groups in the astrocyte cell type

---

## 🧬 Analysis Steps

1. Load expression and metadata tables  
2. Preprocess and harmonise data and structures
3. Perform **pseudobulk differential expression (DE) analysis**  
   - Compare AD vs. Control in astrocytes  
   - Focus on AD-related genes (e.g. `APOE`)  
4. Generate:
   - Volcano plot  
   - Heatmap of selected genes  
5. Validate results against published findings

---

## 🔍 Primary Finding

- `APOE` is found to be **downregulated** in Alzheimer's disease astrocyte subpopulations selected for analysis.

---

## ⚙️ Setup
d
### Clone the repository:
```bash
git clone https://github.com/your-username/ucl_repo.git
cd ucl_repo
```

```bash
#Create and activate the Conda environment:
conda env create -n rnaseq_analysis -f environment.yaml

conda activate rnasea_analysis
```

### Default input files (already included in data/):
```bash
GSE138852_pseudobulk_counts.csv
GSE138852_pseudobulk_metadata.csv
```
- If using your own data, update the file paths in main.py or set the environment variables: EXPR_FILE and META_FILE.

---

## 🚀 Run the Analysis

```bash
python main.py
```

---

## 🔧 Further Development

- Expand analysis to include additional cell types
- Add new visualisation and filtering functions
- Adapt the pipeline to start from raw transcriptomics FASTQ files on HPC systems

---

## 📚 Reference

- Grubman et al., Nature Neuroscience (2019), PMID: 31768052
- Dataset source: GSE138852 (NCBI GEO)
