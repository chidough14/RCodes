🧬 Improving Pathway Interpretation in Colorectal Cancer

This is an R/Bioconductor project for the master’s thesis:
“Improving Pathway Interpretation in Colorectal Cancer: A Comparison of Classical and Correlation-Based Enrichment Methods.”
The analysis compares Over-Representation Analysis (ORA), Gene Set Enrichment Analysis (GSEA), and a correlation-based Gscore approach using TCGA-COAD RNA-seq data.

🧠 Get started
1. Install dependencies

Open R or RStudio and install required packages:

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "DESeq2", "org.Hs.eg.db", "biomaRt",
  "clusterProfiler", "fgsea", "msigdbr", "pheatmap"
))

install.packages(c("ggplot2", "dplyr", "data.table", "knitr"))

2. Run the analysis

Run each script in order from the scripts folder:

Step	Script	Description
1️⃣	data_acquisition.R	Download and preprocess TCGA-COAD RNA-seq data
2️⃣	differential_expression.R	Perform DESeq2 differential expression analysis
3️⃣	enrichment_analysis.R	Conduct ORA using GO and KEGG pathways
4️⃣	gsea.R	Run GSEA on Hallmark gene sets using fgsea
5️⃣	gscore_fgsea.R	Perform correlation-based enrichment (Gscore)
6️⃣	comparison_of_methods.R	Compare ORA, GSEA, and Gscore (Venn, Jaccard, Spearman)
7️⃣	visualization.R	Generate and export all figures and summary tables
3. View results

Results (tables, figures, and .RData files) are automatically saved inside the results directory.
You can view plots directly in RStudio or open the saved .png files.

🧩 Project structure
├── data/                  # Raw and processed TCGA data
├── scripts/               # All R analysis scripts
├── results/               # Tables, plots, and saved RData files
├── report/                # Final thesis report and figures
└── README.md              # Project documentation

🧰 Computational environment

All analyses were performed using R (v4.4.1) and Bioconductor (v3.19)
on Windows 11 Pro (64-bit) with 32 GB RAM and an Intel Core i7 processor.

Key R Packages:
TCGAbiolinks, DESeq2, org.Hs.eg.db, biomaRt,
clusterProfiler, fgsea, msigdbr, pheatmap, ggplot2

Reproducibility was ensured through versioned scripts and fixed random seeds.

🔗 Access the code

All scripts and documentation are available in this GitHub repository:
👉 https://github.com/yourusername/colorectal-cancer-enrichment

🧬 Learn more

If you’d like to reproduce the data retrieval step, refer to the data_acquisition.R script for the
GDCquery() command used to access TCGA-COAD via the Genomic Data Commons (GDC) API.

For background reading:

TCGAbiolinks documentation

clusterProfiler documentation

fgsea documentation

🧠 Citation

If you use this repository or code, please cite:

Ohiaeri, C. (2025). Improving Pathway Interpretation in Colorectal Cancer: A Comparison of Classical and Correlation-Based Enrichment Methods. Master’s Thesis, [University Name], Sweden.