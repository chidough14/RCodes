###############################################################
# COMPARISON OF ENRICHMENT METHODS (ORA, GSEA, Gscore)
# Author: Chidozie Ohiaeri
# Project: Comparison of Enrichment Methods in CRC
# Date: October 2025
###############################################################

# ---- Load libraries ----
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")

library(tidyverse)
library(VennDiagram)
library(ggvenn)
library(pheatmap)

# ---- 1. Load data ----
ora_go <- read.csv("ORA_GO_results.csv")
ora_kegg <- read.csv("ORA_KEGG_results.csv")
gsea <- read.csv("GSEA_fgsea_hallmark_results.csv")
gscore <- read.csv("Gscore_FGSEA_Top50.csv")

# ---- 2. Combine ORA GO and KEGG ----
ora_combined <- bind_rows(
  ora_go %>% select(ID, Description, p.adjust) %>%
    rename(Pathway = Description, padj = p.adjust),
  ora_kegg %>% select(ID, Description, p.adjust) %>%
    rename(Pathway = Description, padj = p.adjust)
) %>%
  distinct(Pathway, .keep_all = TRUE)

# ---- 3. Clean pathway names ----
clean_names <- function(x) gsub("HALLMARK_|_", " ", x)
gsea$Pathway <- clean_names(gsea$pathway)
gscore$Pathway <- clean_names(gscore$pathway)
ora_combined$Pathway <- clean_names(ora_combined$Pathway)

# ---- 4. Identify significant pathways ----
sig_ora <- ora_combined %>% filter(padj < 0.05)
sig_gsea <- gsea %>% filter(padj < 0.05)
sig_gscore <- gscore %>% filter(padj < 0.05)

# ---- 5. Venn Diagram (shared pathways) ----
pathway_lists <- list(
  ORA = sig_ora$Pathway,
  GSEA = sig_gsea$Pathway,
  Gscore = sig_gscore$Pathway
)

# Using ggvenn (cleaner)
ggvenn(pathway_lists,
       fill_color = c("#FF9999", "#99CCFF", "#99FF99"),
       stroke_size = 0.5,
       set_name_size = 4,
       text_size = 4) +
  ggtitle("Venn Diagram of Shared Significant Pathways (ORA, GSEA, Gscore)")

ggsave("Figure13_Venn_PathwayOverlap.png", width = 6, height = 5, dpi = 300)

# ---- 6. Compute Jaccard Similarity ----
# Jaccard = intersection / union for each pair
jaccard <- function(a, b) {
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  return(inter / union)
}

jaccard_mat <- matrix(NA, nrow = 3, ncol = 3,
                      dimnames = list(names(pathway_lists), names(pathway_lists)))

for (i in 1:3) {
  for (j in 1:3) {
    jaccard_mat[i, j] <- jaccard(pathway_lists[[i]], pathway_lists[[j]])
  }
}

# Heatmap for Jaccard
png("Figure14_JaccardSimilarity.png", width = 1200, height = 1000, res = 150)
pheatmap(
  jaccard_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Jaccard Similarity of Enriched Pathways",
  display_numbers = TRUE
)
dev.off()

# 4. Prepare all methods as ranked lists ----
# Rank by adjusted p-value or NES (depending on availability)
get_ora_table <- function(df) {
  df <- df %>%
    rename_with(~ gsub("\\.", "", .x))  # remove dots in column names
  path_col <- grep("Description|Term|Pathway", names(df), value = TRUE)[1]
  pval_col <- grep("^pvalue$|pvalue", names(df), value = TRUE)[1]
  if (is.na(path_col) | is.na(pval_col)) return(NULL)
  df %>%
    select(Pathway = all_of(path_col), pval = all_of(pval_col)) %>%
    filter(!is.na(pval) & is.finite(pval))
}

ora_go_tidy <- get_ora_table(ora_go)
ora_kegg_tidy <- get_ora_table(ora_kegg)
ora_combined <- bind_rows(ora_go_tidy, ora_kegg_tidy)

rank_ora <- setNames(-log10(ora_combined$pval + 1e-10), ora_combined$Pathway)
rank_gsea <- setNames(gsea$NES, gsea$pathway)
rank_gscore <- setNames(gscore$NES, gscore$pathway)

# --- Step 2: Combine into one data frame by pathway name ---
combined_df <- full_join(
  data.frame(Pathway = names(rank_ora), ORA = rank_ora, stringsAsFactors = FALSE),
  data.frame(Pathway = names(rank_gsea), GSEA = rank_gsea, stringsAsFactors = FALSE),
  by = "Pathway"
) %>%
  full_join(
    data.frame(Pathway = names(rank_gscore), GSCORE = rank_gscore, stringsAsFactors = FALSE),
    by = "Pathway"
  )

# --- Step 3: Remove rows with all NAs ---
combined_df <- combined_df %>% filter(!(is.na(ORA) & is.na(GSEA) & is.na(GSCORE)))

# --- Step 4: Compute correlation matrix (pairwise complete cases) ---
spearman_mat <- cor(combined_df[, -1], method = "spearman", use = "pairwise.complete.obs")

# --- Step 5: Plot heatmap ---
col_fun <- colorRampPalette(c("blue", "white", "red"))(99)
breaksList <- seq(-1, 1, length.out = 100)

png("Figure15_SpearmanRankCorrelation.png", width = 1200, height = 1000, res = 150)
pheatmap(
  spearman_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Spearman Rank Correlation of Pathway Rankings",
  display_numbers = TRUE,
  color = col_fun,
  breaks = breaksList
)
dev.off()

# ---- 8. Print summary results ----
cat("Jaccard Similarity Matrix:\n")
print(round(jaccard_mat, 2))

cat("\nSpearman Correlation Matrix:\n")
print(round(spearman_mat, 2))

###############################################################
# End of comparison section
###############################################################
