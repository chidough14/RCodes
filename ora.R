
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(enrichplot)

# ===== 1. Define gene sets =====
load("TCGA_COAD_data_acquisition_Processed.RData") 
sig_genes <- geneTable %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>%
  pull(SYMBOL) %>% 
  unique()

gene_universe <- unique(geneTable$SYMBOL)

# ===== 2. Perform ORA (GO Biological Process) =====
ego <- enrichGO(
  gene          = sig_genes,
  universe      = gene_universe,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# ===== 3. Visualize results =====

dotplot(ego, showCategory = 15, title = "GO Biological Processes (ORA)") +
  theme_minimal(base_size = 13)


barplot(ego, showCategory = 15, title = "Top GO Biological Processes (ORA)") +
  theme_minimal(base_size = 13)

library(kableExtra)
go_table <- as.data.frame(ego) %>%
  arrange(p.adjust) %>%
  head(10) %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)

go_table %>%
  kbl(caption = "Table 4. Top 10 enriched GO Biological Processes (ORA)") %>%
  kable_styling(full_width = FALSE, position = "center", bootstrap_options = "striped")


# Step 4: KEGG pathway enrichment
ekegg <- enrichKEGG(
  gene          = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
  universe      = bitr(gene_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
  organism      = 'hsa',
  pvalueCutoff  = 0.05
)

# Step 5: Visualize KEGG pathways
barplot(ekegg, showCategory = 10, title = "Top KEGG Pathways (ORA)")
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment (ORA)")



kegg_table <- as.data.frame(ekegg) %>%
  arrange(p.adjust) %>%
  head(10) %>%
  dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count)

kegg_table %>%
  kbl(caption = "Table 5. Top 10 enriched KEGG Pathways (ORA)") %>%
  kable_styling(full_width = FALSE, position = "center", bootstrap_options = "striped")

# Step 6: Save results
write.csv(as.data.frame(ego), "ORA_GO_results.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "ORA_KEGG_results.csv", row.names = FALSE)

# Appendix Figure A6 — GO Enrichment Map
ego_sim <- pairwise_termsim(ego)
png("Appendix_Figure_A6—GO_Enrichment_Map.png", 
    width = 3000,    # width in pixels
    height = 2000,   # height in pixels
    res = 300)       # resolution in dpi for clarity
# Plot the emap
emapplot(ego_sim, showCategory = min(30, length(ego$Description)))  # ensure showCategory <= #terms
dev.off()  # close the PNG device

#Appendix Table B5. Complete GO Biological Process enrichment results from ORA.
write.csv(as.data.frame(ego),
          "Appendix_Table_B5_ORA_GO_All.csv",
          row.names = FALSE)
#TableB5 preview
selected_metadata <- ego@result %>%
  dplyr::select(ID, Description, pvalue, p.adjust, Count)

preview_B5 <- head(selected_metadata, 10)
write.csv(as.data.frame(preview_B5),
          "Appendix_Table_B5_ORA_GO_Preview.csv",
          row.names = FALSE)

#Appendix Table B6. Complete KEGG pathway enrichment results from ORA.
write.csv(as.data.frame(ekegg),
          "Appendix_Table_B6_ORA_KEGG_All.csv",
          row.names = FALSE)
#TableB6 preview
selected_metadata_ekegg <- ekegg@result %>%
  dplyr::select(ID, Description, pvalue, p.adjust, Count)

preview_B6 <- head(selected_metadata_ekegg, 10)
write.csv(as.data.frame(preview_B6),
          "Appendix_Table_B6_ORA_KEGG_Preview.csv",
          row.names = FALSE)


