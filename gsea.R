
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(enrichplot)

# Load your processed data
load("TCGA_COAD_data_acquisition_Processed.RData")

# Prepare the Ranked Gene List
gene_list <- geneTable$log2FoldChange
names(gene_list) <- geneTable$SYMBOL
gene_list <- sort(gene_list, decreasing = TRUE)

# Download the entire C2 (Curated Pathways) collection
all_c2 <- msigdbr(species = "Homo sapiens", category = "C2")

kegg_msig <- all_c2 %>% 
  filter(gs_subcat == "CP:KEGG" | grepl("KEGG", gs_name))

# Check if it worked
print(paste("Number of KEGG pathways found:", length(unique(kegg_msig$gs_name))))

#  Convert to list for GSEA
kegg_pathways <- split(x = kegg_msig$gene_symbol, f = kegg_msig$gs_name)


gene_df <- data.frame(Symbol = names(gene_list), Log2FC = gene_list)


gene_df_unique <- gene_df %>%
  group_by(Symbol) %>%
  summarize(Log2FC = mean(Log2FC)) %>%
  filter(Symbol != "" & !is.na(Symbol)) # Remove any empty symbols


gene_list_final <- gene_df_unique$Log2FC
names(gene_list_final) <- gene_df_unique$Symbol


gene_list_final <- sort(gene_list_final, decreasing = TRUE)


gsea_results <- fgsea(
  pathways = kegg_pathways, 
  stats    = gene_list_final, 
  minSize  = 15,
  maxSize  = 500
)

gsea_results_for_csv <- gsea_results %>%
  as.data.frame() %>%
  dplyr::select(-leadingEdge) 

# 2. Save to CSV
write.csv(gsea_results_for_csv, "GSEA_KEGG_Results_Final.csv", row.names = FALSE)

gsea_table_short <- gsea_results_for_csv %>%
  arrange(padj) %>%
  select(pathway, NES, pval, padj) %>%
  head(10)

knitr::kable(gsea_table_short, caption = "Table 3. Top 10 enriched KEGG pathways.")

# ==============================================================================
# VISUALIZATION FOR REPORT
# ==============================================================================

#Top 15 KEGG Pathways Bar Plot
top_15_gsea <- head(gsea_results, 15)

ggplot(top_15_gsea, aes(x = reorder(pathway, NES), y = NES)) +
  geom_bar(stat = "identity", aes(fill = padj < 0.05)) +
  coord_flip() +
  scale_fill_manual(values = c("grey", "#E69F00")) +
  theme_minimal() +
  labs(
    title = "Figure 11. Top Enriched KEGG Pathways (GSEA)",
    x = "Pathway",
    y = "Normalized Enrichment Score (NES)"
  )

# Specific Enrichment Plot (e.g., for the top pathway)
top_pathway_name <- gsea_results$pathway[1]



plotEnrichment(kegg_pathways[[top_pathway_name]], gene_list) +
  labs(title = paste("Figure 10. GSEA Enrichment Plot:", top_pathway_name))

# Save workspace for the final comparison step
save(gsea_results, kegg_pathways, file = "GSEA_Final_Results.RData")