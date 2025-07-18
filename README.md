🧠 UCL Interview – Transcriptomics Analysis of Alzheimer's Disease
🧪 Project Overview
A personal practice draft project on bulk/pseudobulk transcriptomics analysis, focusing on differential gene expression in astrocytes from the human entorhinal cortex using public snRNA-seq data from Alzheimer's disease (AD) and control samples.

📂 Dataset
Source: NCBI GEO – GSE138852
Human entorhinal cortex
Single-nucleus RNA-seq (snRNA-seq)
Raw counts aggregated into pseudobulk format for AD vs. control groups in the astrocyte cell type

🧬 Analysis Steps
Load expression and metadata tables
Preprocess and harmonise sample IDs
Perform pseudobulk differential expression (DE) analysis (AD vs. Control in astrocytes), with a focus on AD-related genes (e.g. APOE)
Generate:
Volcano plot
Heatmap of selected genes
Validate against published findings for hands-on learning

Primary finding
APOE is DOWN regulated in Alzheimer's disease astrocyte subpopulations we selected

⚙️ Setup
#download the github repo
git clone https://github.com/your-username/ucl_repo.git
cd ucl_repo

#create the conda env (name: rnasea_analysis) for this pipeline and activate it
conda env create -f environment.yml
conda activate rnasea_analysis

default data files are already in data/:
If using your own data, update the file paths in main.py or set environment variables EXPR_FILE and META_FILE.

🚀 Run the Analysis
python main.py

🔧 Further Development

Expand analysis to include additional cell types
Add new visualisation and filtering functions
Adapt the pipeline to start from raw FASTQ files using high-performance computing resources

📚 Reference

Grubman et al., Nature Neuroscience (2019), PMID: 31768052
