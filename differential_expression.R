library(DESeq2)

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

# Summary of results
summary(res)

# Filter significant DEGs
sig_res <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)
nrow(sig_res)

plotMA(res, ylim = c(-5, 5), main = "MA Plot: Tumor vs Normal")

library(EnhancedVolcano)

#Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
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

library(pheatmap)
library(RColorBrewer)

#Heatmap
top_genes <- head(order(res$padj), 50)
mat <- assay(vst(dds))[top_genes, ]
mat <- t(scale(t(mat)))

annotation_col <- data.frame(SampleType = colData(dds)$sample_type)
rownames(annotation_col) <- colnames(mat)

pheatmap(mat,
         annotation_col = annotation_col,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize = 10,
         main = "Top 50 Differentially Expressed Genes"
)

save(dds, vsd, data_filtered, res, sig_res, file = "TCGA_COAD_data_acquisition_Processed.RData")

