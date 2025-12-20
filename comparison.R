# 1. LOAD AND CLEAN ORA DATA
ora_df <- read.csv("ORA_KEGG_results.csv", stringsAsFactors = FALSE)
# We use Description because it's human-readable (e.g., "ABC transporters")
ora_df$CleanName <- tolower(ora_df$Description)

# 2. LOAD AND CLEAN GSEA DATA
gsea_df <- as.data.frame(gsea_results)
# GSEA names often look like "KEGG_ABC_TRANSPORTERS"
# We remove "KEGG_", replace underscores with spaces, and make lowercase
gsea_df$CleanName <- gsub("KEGG_", "", gsea_df$pathway)
gsea_df$CleanName <- gsub("_", " ", gsea_df$CleanName)
gsea_df$CleanName <- tolower(gsea_df$CleanName)

# 3. IDENTIFY COMMON PATHWAYS BY NAME
common_names <- intersect(ora_df$CleanName, gsea_df$CleanName)

# 4. CALCULATE JACCARD SIMILARITY
union_names <- union(ora_df$CleanName, gsea_df$CleanName)
jaccard_val <- length(common_names) / length(union_names)

message("Fixed Jaccard Similarity Index: ", round(jaccard_val, 3))

# 5. MERGE FOR PLOTTING
# If common_names is still 0, we will print an error message
if(length(common_names) > 0) {
  ora_sub <- ora_df[ora_df$CleanName %in% common_names, c("CleanName", "p.adjust")]
  gsea_sub <- gsea_df[gsea_df$CleanName %in% common_names, c("CleanName", "NES")]
  
  comp_df <- merge(ora_sub, gsea_sub, by = "CleanName")
  
  # 6. FIGURE 14: RANK COMPARISON PLOT
  library(ggplot2)
  ggplot(comp_df, aes(x = rank(p.adjust), y = rank(-abs(NES)))) +
    geom_point(color = "darkorange", size = 5, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
    labs(title = "Rank Consistency: ORA vs. GSEA",
         subtitle = paste("Jaccard Index =", round(jaccard_val, 3)),
         x = "ORA Significance Rank (p.adjust)",
         y = "GSEA Importance Rank (|NES|)") +
    theme_minimal(base_size = 16)
} else {
  message("STILL 0 OVERLAP. Trying fallback: checking first 5 names...")
  print(head(ora_df$Description))
  print(head(gsea_df$pathway))
}

# Calculate Spearman's Rank Correlation
# We use negative abs(NES) so that a higher magnitude is a "better" rank (Rank 1)
res <- cor.test(rank(comp_df$p.adjust), 
                rank(-abs(comp_df$NES)), 
                method = "spearman")

# Print the results
message("Spearman Correlation (rho): ", round(res$estimate, 3))
message("P-value: ", format.pval(res$p.value))


# Calculate rho for the plot subtitle
rho_val <- cor(rank(comp_df$p.adjust), rank(-abs(comp_df$NES)), method = "spearman")

# Generate Plot
ggplot(comp_df, aes(x = rank(p.adjust), y = rank(-abs(comp_df$NES)))) +
  geom_point(color = "darkorange", size = 5, alpha = 0.7) +
  geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  annotate("text", x = max(rank(comp_df$p.adjust))*0.2, 
           y = max(rank(-abs(comp_df$NES))), 
           label = paste("Spearman rho =", round(rho_val, 3)), 
           size = 3, fontface = "bold") +
  labs(title = "Rank Consistency: ORA vs. GSEA",
       x = "ORA Significance Rank (Lower p-value = Higher Rank)",
       y = "GSEA Importance Rank (Higher |NES| = Higher Rank)") +
  theme_bw(base_size = 12) 