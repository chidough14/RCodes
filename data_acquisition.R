
# ===================== 1. Load Required Packages =====================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "TCGAbiolinks", "SummarizedExperiment", "DESeq2",
  "biomaRt", "org.Hs.eg.db", "pheatmap", "RColorBrewer",
  "matrixStats"
))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)

# ===================== 2. Data Download =====================
# Define query for TCGA-COAD (HTSeq-count workflow)
query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download (may take a while — large dataset!)
GDCdownload(query, method = "api", files.per.chunk = 20)

# Prepare data
data_coad <- GDCprepare(query)
save(data_coad, file = "TCGA_COAD_data.RData")

# ===================== 3. Dataset Overview =====================
# Check sample types
sample_types <- data_coad$shortLetterCode
table(sample_types)

# Count tumor and normal samples
sample_counts <- table(sample_types)

# Number of genes
total_genes <- nrow(assay(data_coad))

# Mean read depth (mean of total counts per sample)
mean_read_depth <- mean(colSums(assay(data_coad)))

dataset_summary <- data.frame(
  Metric = c("Number of Tumor Samples", "Number of Normal Samples",
             "Total Genes Detected", "Mean Read Depth per Sample"),
  Value = c(sample_counts["TP"], sample_counts["NT"],
            total_genes, round(mean_read_depth, 2))
)

write.csv(dataset_summary, "Table1_Dataset_Summary.csv", row.names = FALSE)
print(dataset_summary)

# ===================== 4. Data Filtering =====================

# ---- Filter low-count genes ----
# Keep genes with at least 10 counts in at least 5 samples
keep <- rowSums(assay(data_coad) >= 10) >= 5
data_filtered <- data_coad[keep, ]

# ---- Summary statistics ----
initial_genes  <- nrow(data_coad)
retained_genes <- sum(keep)
excluded_genes <- initial_genes - retained_genes

# ---- Percentages ----
percent_initial  <- round((initial_genes / initial_genes) * 100, 1)
percent_excluded <- round((excluded_genes / initial_genes) * 100, 1)
percent_retained <- round((retained_genes / initial_genes) * 100, 1)

# ---- Construct Table 2 ----
table2 <- data.frame(
  Filtering_Step      = c("Initial dataset", "Excluded genes", "Retained genes"),
  Description         = c(
    "All annotated genes obtained from HTSeq-Counts",
    "Genes with fewer than 10 counts in fewer than 5 samples",
    "Genes passing expression threshold for analysis"
  ),
  Number_of_Genes     = c(initial_genes, excluded_genes, retained_genes),
  Percentage_of_Total = c(percent_initial, percent_excluded, percent_retained),
  stringsAsFactors = FALSE,
  check.names = TRUE
)

# ---- Save as CSV ----
write.csv(table2, "Table2_FilteringSummary.csv", row.names = FALSE)

# ---- Print short summary ----
cat("Genes before filtering:", initial_genes, "\n")
cat("Genes excluded:", excluded_genes, "\n")
cat("Genes retained:", retained_genes, "\n")

# ===================== 5. Normalization =====================
dds <- DESeqDataSet(data_filtered, design = ~shortLetterCode)
dds <- DESeq(dds)

# Variance-stabilizing transformation
vsd <- vst(dds, blind = FALSE)

library(reshape2)

# Extract raw and VST counts
raw_counts <- assay(data_filtered)
vst_counts <- assay(vsd)

# Subsample to reduce plotting size
set.seed(42)
sample_genes <- sample(rownames(raw_counts), 1000)

# Convert to long format for ggplot
raw_melt <- melt(log10(raw_counts[sample_genes, ] + 1))
vst_melt <- melt(vst_counts[sample_genes, ])

# Add a column indicating data type
raw_melt$Type <- "Raw"
vst_melt$Type <- "Normalized"

# Combine both data frames
combined_df <- rbind(raw_melt, vst_melt)

# Plot density
library(ggplot2)
ggplot(combined_df, aes(x = value, fill = Type, color = Type)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("Raw" = "#1F78B4", "Normalized" = "#E69F00")) +
  scale_color_manual(values = c("Raw" = "#1F78B4", "Normalized" = "#E69F00")) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Distribution of Raw vs Normalized Counts (DESeq2)",
    x = "Log10(Expression + 1)",
    y = "Density"
  )

# ===================== 6. Quality Control =====================

## PCA Plot
plotPCA(vsd, intgroup = "shortLetterCode") +
  ggtitle("PCA of TCGA-COAD Samples (VST normalized)")

## Sample Distance Heatmap

ann_colors <- list(SampleType = c(
  "NT" = "#33A02C",  # Normal Tissue
  "TP" = "#1F78B4",  # Primary Tumor
  "TR" = "#E31A1C",  # Recurrent Tumor
  "TM" = "#FF7F00"   # Metastatic Tumor
))

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
annotation_col <- data.frame(SampleType = vsd$shortLetterCode)
rownames(annotation_col) <- colnames(vsd)

# pheatmap(
#   sample_dist_matrix,
#   annotation_col = annotation_col,
#   main = "Sample-to-Sample Distance Heatmap"
# )
pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)
)

# ===================== 7. Gene Expression Patterns =====================
# Top 50 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat <- assay(vsd)[topVarGenes, ] 
mat <- mat - rowMeans(mat) 
# pheatmap( mat, 
#           annotation_col = annotation_col, 
#           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255), 
#           show_rownames = TRUE, 
#           main = "Heatmap of Top 50 Variable Genes in TCGA-COAD" 
# )
pheatmap(
  mat,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize = 10
)

# ===================== 8. Save Workspace =====================
save(dds, vsd, data_filtered, file = "TCGA_COAD_data_acquisition_Processed.RData")

