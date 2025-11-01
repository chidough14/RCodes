# Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
cran_pkgs <- c("dplyr", "tibble", "ggplot2", "data.table", "pheatmap")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# Bioconductor packages
bioc_pkgs <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi",
               "fgsea", "msigdbr")  # msigdbr fetches gene sets from MSigDB; internet required
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

# load libraries
library(DESeq2)
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(data.table)
library(pheatmap)
library(magrittr)  

# ------------- PREP if needed ----------------
load("TCGA_COAD_data_acquisition_Processed.RData") 
if (!exists("dds")) stop("dds object not found. Load your DESeqDataSet named 'dds' before running enrichment.")

# run DESeq if not yet run
if (!exists("res")) {
  if(!"DESeq" %in% loadedNamespaces()) library(DESeq2)
  dds <- DESeq(dds)
  res <- results(dds, name = "shortLetterCode_TP_vs_NT")
  
  # Shrink log2 fold changes
  res <- lfcShrink(dds, coef = "shortLetterCode_TP_vs_NT", type = "apeglm")
}

# VST for expression-based methods
if (!exists("vst_mat")) {
  vst_obj <- vst(dds, blind = FALSE)
  vst_mat <- assay(vst_obj) # genes x samples
}

# Ensure rownames are Ensembl IDs (or SYMBOL). If Ensembl have versions (ENSG... .1) strip versions for mapping.
ensembl_ids <- rownames(res)
ensembl_ids_simple <- sub("\\..*$", "", ensembl_ids)

# Map Ensembl -> SYMBOL
map <- mapIds(org.Hs.eg.db,
              keys = ensembl_ids_simple,
              column = "SYMBOL",
              keytype = "ENSEMBL",
              multiVals = "first")

# Build gene table
geneTable <- as.data.frame(res) %>%
  rownames_to_column(var = "ensembl") %>%
  mutate(ensembl_simple = sub("\\..*$", "", ensembl),
         SYMBOL = map[ensembl_simple]) %>%
  arrange(padj)

# Remove genes without symbol if you want
geneTable <- geneTable %>% filter(!is.na(SYMBOL))
save(dds, vsd, data_filtered, res, sig_res, geneTable, file = "TCGA_COAD_data_acquisition_Processed.RData")
