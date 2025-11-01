# Improving Pathway Interpretation in Colorectal Cancer

This is an R/Bioconductor project for the master’s thesis:
“Improving Pathway Interpretation in Colorectal Cancer: A Comparison of Classical and Correlation-Based Enrichment Methods.”
The analysis compares Over-Representation Analysis (ORA), Gene Set Enrichment Analysis (GSEA), and a correlation-based Gscore approach using TCGA-COAD RNA-seq data.

## Get started
1. Clone the repository

```bash
git clone https://github.com/chidough14/RCodes.git
cd RCodes
```

2. Install dependencies

Open R or RStudio and install required packages:

```bash
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "DESeq2", "org.Hs.eg.db", "biomaRt",
  "clusterProfiler", "fgsea", "msigdbr", "pheatmap"
))

install.packages(c("ggplot2", "dplyr", "data.table", "knitr"))
```

3. Run the analysis

Run each script in order from the scripts folder:

1. **data_acquisition.R:**	Download and preprocess TCGA-COAD RNA-seq data
2. **differential_expression.R:**	Perform DESeq2 differential expression analysis
3. **enrichment_analysis.R:**	Generates a ranked gene table that is used by the ORA, GSEA, and Gscore
4. **ora.R:**	Conduct ORA using GO and KEGG pathways
5. **gsea.R:**	Run GSEA on Hallmark gene sets using fgsea
6. **gscore_fgsea.R:**	Perform correlation-based enrichment (Gscore)
7. **comparison_of_methods.R:**	Compare ORA, GSEA, and Gscore (Venn, Jaccard, Spearman)

## View results

Results (tables, figures) are automatically saved inside the tables/figures directory.
You can view plots directly in RStudio or open the saved .png files.

## Computational environment

All analyses were performed using R (v4.4.1) and Bioconductor (v3.19)
on Windows 11 Pro (64-bit) with 32 GB RAM and an Intel Core i7 processor.

Key R Packages:
TCGAbiolinks, DESeq2, org.Hs.eg.db, biomaRt,
clusterProfiler, fgsea, msigdbr, pheatmap, ggplot2