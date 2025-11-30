
# Install / load packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_cran <- c("data.table", "dplyr", "ggplot2", "DT")
for (p in required_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
bioc_pkgs <- c("fgsea", "msigdbr")
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p)

library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)
library(ggplot2)
library(DT)
load("TCGA_COAD_data_acquisition_Processed.RData") 

# --- Prepare pathways: Hallmark (H) collection from MSigDB
halls <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# --- Prepare ranked vector (names = SYMBOL)
# use log2FoldChange as ranking metric (descending -> upregulated at top)
ranks_df <- geneTable %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  distinct(SYMBOL, .keep_all = TRUE)  

ranks <- ranks_df$log2FoldChange
names(ranks) <- ranks_df$SYMBOL


# --- Run fgsea
set.seed(42)
fgsea_res <- fgsea::fgsea(pathways = halls,
                          stats = ranks,
                          minSize = 15,
                          maxSize = 500,
                          nperm = 10000)

# Tidy results
fgsea_df <- as.data.frame(fgsea_res) %>%
  arrange(padj)

# Convert fgsea results to a data frame and simplify list columns
fgsea_df <- as.data.frame(fgsea_res)

# Collapse list columns into comma-separated strings
fgsea_df$leadingEdge <- sapply(fgsea_df$leadingEdge, paste, collapse = ",")

write.csv(fgsea_df, "GSEA_fgsea_hallmark_results.csv", row.names = FALSE)

#Truncated table
fgsea_table_short <- fgsea_df %>%
  arrange(padj) %>%
  select(pathway, NES, pval, padj) %>%
  head(10)

write.csv(fgsea_table_short, "Table 3. GSEA_fgsea_table_short.csv", row.names = FALSE)

knitr::kable(fgsea_table_short, caption = "Table 3. Top 10 enriched Hallmark pathways (GSEA).")


# --- Plot: enrichment plot for the top positively enriched pathway
if (nrow(fgsea_df) > 0) {
  pos_top <- fgsea_df %>% filter(NES > 0) %>% arrange(padj) %>% slice_head(n = 1) %>% pull(pathway)
  if (!is.na(pos_top)) {
    p1 <- fgsea::plotEnrichment(halls[[pos_top]], ranks) +
      ggtitle(paste0("GSEA (fgsea) enrichment: ", pos_top)) +
      theme_minimal(base_size = 13)
    print(p1)
    ggsave(filename = "Figure11_GSEA_enrichment_top_pos.png", plot = p1, width = 7, height = 5, dpi = 300)
  }
  
  # top negatively enriched pathway
  neg_top <- fgsea_df %>% filter(NES < 0) %>% arrange(padj) %>% slice_head(n = 1) %>% pull(pathway)
  if (!is.na(neg_top)) {
    p2 <- fgsea::plotEnrichment(halls[[neg_top]], ranks) +
      ggtitle(paste0("GSEA (fgsea) enrichment: ", neg_top)) +
      theme_minimal(base_size = 13)
    print(p2)
    ggsave(filename = "GSEA_enrichment_top_neg.png", plot = p2, width = 7, height = 5, dpi = 300)
  }
}

# --- Barplot of top pathways by NES (top 15 by padj)
topN <- 15
topbar <- fgsea_df %>% dplyr::slice_head(n = topN)
pbar <- ggplot(topbar, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
  geom_col() + coord_flip() +
  theme_minimal(base_size = 12) +
  labs(x = "", y = "Normalized Enrichment Score (NES)",
       title = "Top Hallmark pathways (GSEA - fgsea)") +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2"), guide = FALSE)
print(pbar)
ggsave(filename = "Figure12_GSEA_topN_NES_barplot.png", plot = pbar, width = 8, height = 6, dpi = 300)


# saveRDS(fgsea_df, file = "GSEA_fgsea_results.rds")
