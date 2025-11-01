# Load required libraries
library(DESeq2)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(AnnotationDbi)
library(org.Hs.eg.db)

load("TCGA_COAD_data_acquisition_Processed.RData") 

#--- 1. VST transformation ---
if (!exists("vst_mat")) {
  vst_obj <- vst(dds, blind = FALSE)
  vst_mat <- assay(vst_obj)
}

#--- 2. Create phenotype vector (Normal = 0, Tumor = 1) ---
pheno_vec <- as.data.frame(colData(dds))$sample_type
pheno_num <- ifelse(pheno_vec == "Primary Tumor", 1, 0)
names(pheno_num) <- colnames(vst_mat)

#--- 3. Ensure dimensions and order match ---
pheno_num <- pheno_num[colnames(vst_mat)]

#--- 4. Gene annotation: map ENSEMBL IDs to gene symbols ---
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gsub("\\..*", "", rownames(vst_mat)),  # remove version numbers
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Add gene symbols as rownames
rownames(vst_mat) <- gene_symbols

#--- 5. Compute Spearman correlation for each gene vs phenotype ---
gene_corr <- apply(vst_mat, 1, function(x) cor(x, pheno_num, method = "spearman", use = "complete.obs"))

# Create dataframe of gene correlations
gene_corr_df <- data.frame(
  SYMBOL = names(gene_corr),
  corr = gene_corr,
  stringsAsFactors = FALSE
)

# Remove rows with NA SYMBOLs
gene_corr_df <- gene_corr_df[!is.na(gene_corr_df$SYMBOL), ]

#--- 6. Prepare ranked list for GSEA ---
gscore_rank <- gene_corr_df$corr
names(gscore_rank) <- gene_corr_df$SYMBOL

# Remove duplicates (keep most correlated one)
gscore_rank <- tapply(gscore_rank, names(gscore_rank), function(x) x[which.max(abs(x))])
gscore_rank <- sort(unlist(gscore_rank), decreasing = TRUE)

#--- 7. Load Hallmark gene sets ---
m_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(m_df$gene_symbol, m_df$gs_name)

set.seed(42)

fgsea_res_gscore <- fgsea(
  pathways = hallmark_list,
  stats = gscore_rank,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

#--- 9. Sort and keep all results ---
fgsea_res_full <- fgsea_res_gscore %>%
  arrange(padj) %>%
  select(pathway, NES, pval, padj)



# Save top 50 for report
fgsea_table50 <- fgsea_res_full %>% head(50)
write.csv(fgsea_table50, "Gscore_FGSEA_Top50.csv", row.names = FALSE)

# Optional: short top 10 just for quick preview
fgsea_res_top10 <- fgsea_res_full %>% head(10)
# Order by adjusted p-value (FDR)

write.csv(fgsea_res_top10, "Table 4. Gscore_FGSEA_Top10.csv", row.names = FALSE)
knitr::kable(fgsea_res_top10, caption = "Table 4. Gscore_FGSEA_Top10.csv")
view(fgsea_res_top10)

# 11. Visualization
library(ggplot2)

sig_fgsea <- fgsea_res_gscore %>%
  filter(padj < 0.05) %>%
  arrange(NES)

ggplot(sig_fgsea, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Significantly Enriched Hallmark Pathways (FGSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  )

# Example: Plot a single pathway
plotEnrichment(hallmark_list[["HALLMARK_KRAS_SIGNALING_UP"]], gscore_rank)

# pos_top <- fgsea_res_full %>% filter(NES > 0) %>% arrange(padj) %>% slice_head(n = 1) %>% pull(pathway)
p2 <- plotEnrichment(hallmark_list[["HALLMARK_KRAS_SIGNALING_UP"]], gscore_rank) +
  # ggtitle(paste0("GSCORE (fgsea) enrichment: ", pos_top)) +
  theme_minimal(base_size = 13)
print(p2)
ggsave(filename = "GSCORE_enrichment.png", plot = p2, width = 7, height = 5, dpi = 300)

#Save results
saveRDS(fgsea_res, file = "GSEA_fgsea_results.rds")
saveRDS(gene_corr_df, file = "GSEA_gene_correlation_results.rds")
