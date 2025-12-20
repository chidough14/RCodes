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
library(RColorBrewer)

if (!exists("dds")) {
  if (file.exists("TCGA_COAD_data_acquisition_Processed.RData")) {
    load("TCGA_COAD_data_acquisition_Processed.RData")
    cat("Loaded 'dds' from TCGA_COAD_data_acquisition_Processed.RData\n")
  } else {
    stop("Error: 'dds' object not found. Please run Data_Acquisition.R first.")
  }
} else {
  cat("'dds' already exists in the environment. Proceeding with analysis...\n")
}

# ===================== Filter Low-Count Genes =====================
# dds <- dds[rowSums(counts(dds)) > 10, ]

# ===================== Run DESeq =====================

# Run DESeq only if not already done
# if (!"DESeqDataSet" %in% class(dds) || !"DESeqResults" %in% ls()) {
#   dds <- DESeq(dds)
#   cat("DESeq model fitted successfully.\n")
# } else {
#   cat("DESeq analysis already completed. Skipping re-run.\n")
# }

# Check the coefficients
resultsNames(dds)

# Extract results using the correct coefficient name
res <- results(dds, name = "shortLetterCode_TP_vs_NT")

# Shrink log2 fold changes
res <- lfcShrink(dds, coef = "shortLetterCode_TP_vs_NT", type = "apeglm")

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

# Summary of results
summary(res)

# Appendix_Figure_A2_CooksDistance
png("Appendix_Figure_A2_CooksDistance.png", width = 1800, height = 1600, res = 220)
# Faster Cook's distance plot by limiting to top 2000 most variable genes
topVar <- head(order(rowVars(counts(dds)), decreasing = TRUE), 2000)
plot(log10(assays(dds)[["cooks"]][topVar, ] + 1),
     main = "Cook's Distance (Top 2,000 Most Variable Genes)",
     xlab = "Sample",
     ylab = "log10(Cook's Distance + 1)")

dev.off()


# Filter significant DEGs
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)
nrow(sig_res)

plotMA(res, ylim = c(-5, 5), main = "MA Plot: Tumor vs Normal")

library(EnhancedVolcano)

# Attach gene symbols from geneTable to DESeq2 results
res$ensembl <- rownames(res)

# Merge res with geneTable to add SYMBOL column
res_merged <- merge(as.data.frame(res), 
                    geneTable[, c("ensembl", "SYMBOL")],
                    by = "ensembl",
                    all.x = TRUE)

# Replace rownames for convenience
rownames(res_merged) <- res_merged$ensembl

#Volcano plot
EnhancedVolcano(res,
                lab = res_merged$SYMBOL,
                x = "log2FoldChange",
                y = "padj",
                xlim = c(-6, 6),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 3,
                title = "Volcano Plot: Tumor vs Normal",
                subtitle = "TCGA-COAD Differential Expression"
)

#Appendix Table B3. Complete DESeq2 differential expression results for tumor vs normal samples.
write.csv(as.data.frame(res_merged),
          "Appendix_Table_B3_DESeq2_All_Genes.csv",
          row.names = FALSE)

#TableB3 preview
selected_metadata <- res_merged %>%
  dplyr::select(
    ensembl,
    log2FoldChange,
    pvalue,
    padj,
    SYMBOL,
  )
preview_B3 <- head(selected_metadata, 10)
write.csv(as.data.frame(preview_B3),
          "Appendix_Table_B3_DESeq2_Preview.csv",
          row.names = FALSE)

#Appendix Table B4. Full list of significantly differentially expressed genes.
write.csv(as.data.frame(sig_res),
          "Appendix_Table_B4_Significant_DEGs.csv",
          row.names = TRUE)
#Table b4 preview
preview_B4 <- head(sig_res, 10)
write.csv(as.data.frame(preview_B4),
          "Appendix_Table_B4_Significant_DEGs_Preview.csv",
          row.names = TRUE)


library(pheatmap)
library(RColorBrewer)

#Heatmap
top_genes <- head(order(res$padj), 50)

# ---- Extract VST expression matrix ----
vst_mat <- assay(vst(dds))[top_genes, ]

# ----Map rownames (ENSEMBL) to gene symbols ----
ensembl_ids <- rownames(vst_mat)
ensembl_simple <- sub("\\..*$", "", ensembl_ids)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = ensembl_simple,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Replace NA symbols with ENSEMBL IDs (to avoid empty labels)
gene_symbols[is.na(gene_symbols)] <- ensembl_simple[is.na(gene_symbols)]

# Assign mapped symbols as rownames
rownames(vst_mat) <- gene_symbols

# ---- Z-score transform rows ----
mat <- t(scale(t(vst_mat)))

# --- Annotation for columns ----
annotation_col <- data.frame(SampleType = colData(dds)$sample_type)
rownames(annotation_col) <- colnames(mat)

# ---- Plot heatmap ----
pheatmap(
  mat,
  annotation_col = annotation_col,
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 10,
  main = "Top 50 Differentially Expressed Genes (Gene Symbols)"
)

# Appendix Figure A1 — Dispersion Estimates
png("Appendix_Figure_A1_DispersionEstimates.png", 
    width = 1800, height = 1600, res = 220)
plotDispEsts(dds, main = "Dispersion Estimates (DESeq2)")
dev.off()

save(dds, vsd, data_filtered, res, sig_res, res_merged, geneTable, file = "TCGA_COAD_data_acquisition_Processed.RData")

