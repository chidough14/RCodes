

# Required packages
needed <- c("msigdbr", "data.table", "dplyr", "stringr", "AnnotationDbi",
            "org.Hs.eg.db", "ggvenn", "pheatmap", "RColorBrewer", "ggplot2")
for (p in needed) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
# Bioconductor packages (if not installed)
if (!requireNamespace("msigdbr", quietly = TRUE)) BiocManager::install("msigdbr")

library(msigdbr)
library(data.table)
library(dplyr)
library(stringr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggvenn)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

#Load data
ora_go_file    <- "ORA_GO_results.csv"   
ora_kegg_file  <- "ORA_KEGG_results.csv"  
gsea_file      <- "GSEA_fgsea_hallmark_results.csv"  
gscore_file    <- "Gscore_FGSEA_Top50.csv" 

#Read ORA results and extract pathway -> gene list
read_ora_table <- function(path) {
  df <- fread(path, data.table = FALSE)
  genes_col <- NULL
  if ("geneID" %in% names(df)) genes_col <- "geneID"
  if (is.null(genes_col)) {
    possible <- c("geneID","gene_id","gene","genes")
    genes_col <- intersect(possible, names(df))[1]
  }
  if (is.na(genes_col) || is.null(genes_col)) stop("ORA file has no gene column named 'geneID' or similar.")
  list(df = df, gene_col = genes_col)
}

go_in  <- read_ora_table(ora_go_file)
kegg_in<- read_ora_table(ora_kegg_file)

# helper to convert geneID field to SYMBOL vector:
parse_gene_field <- function(gstring) {
  if (is.na(gstring) || gstring == "") return(character(0))
  pieces <- unlist(str_split(gstring, pattern = "[/;,]"))
  pieces <- str_trim(pieces)
  pieces <- pieces[pieces != ""]
  pieces
}

# Build named list: ORA_pathway_name -> character vector of gene symbols (or Entrez IDs)
build_ora_list <- function(ora_obj, source_label = "GO") {
  df <- ora_obj$df
  gene_col <- ora_obj$gene_col
  path_col <- ifelse("Description" %in% names(df), "Description",
                     ifelse("ID" %in% names(df), "ID", names(df)[1]))
  ora_list <- list()
  for (i in seq_len(nrow(df))) {
    pname <- paste0(source_label, ":", df[[path_col]][i])
    raw_genes <- parse_gene_field(df[[gene_col]][i])
    ora_list[[pname]] <- unique(raw_genes)
  }
  ora_list
}
ora_go_list   <- build_ora_list(go_in, "GO")
ora_kegg_list <- build_ora_list(kegg_in, "KEGG")
ora_all_list  <- c(ora_go_list, ora_kegg_list)


example_genes <- unlist(ora_all_list[1:min(5, length(ora_all_list))])

is_entrez_like <- all(grepl("^[0-9]+$", example_genes[ example_genes != "" ][1:min(length(example_genes),20) ]))

# If Entrez IDs, map to SYMBOL
if (is_entrez_like) {
  message("Detected numeric gene IDs (likely Entrez). Mapping to SYMBOL using org.Hs.eg.db ...")
  ora_all_list <- lapply(ora_all_list, function(vec) {
    mapped <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys = as.character(vec),
                                    column = "SYMBOL",
                                    keytype = "ENTREZID",
                                    multiVals = "first")
    unique(na.omit(as.character(mapped)))
  })
} else {
  ora_all_list <- lapply(ora_all_list, function(vec) unique(toupper(str_trim(vec[vec!=""]))))
}


#Load Hallmark gene sets (msigdbr)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(toupper(m_df$gene_symbol), m_df$gs_name)  # uppercase symbols for consistent matching

#Map each ORA pathway to Hallmark by gene-overlap (Jaccard > 0)
jaccard <- function(a,b){
  a <- unique(a); b <- unique(b)
  if (length(a)==0 || length(b)==0) return(0)
  inter <- length(intersect(a,b))
  union  <- length(union(a,b))
  inter/union
}

#For each ORA pathway, compute Jaccard to each hallmark and keep matches > 0
ora_to_hallmarks <- lapply(names(ora_all_list), function(pn){
  pgenes <- ora_all_list[[pn]]
  scores <- sapply(hallmark_list, function(hg) jaccard(pgenes, hg))
  names(scores)[which(scores > 0)]
})
names(ora_to_hallmarks) <- names(ora_all_list)


ora_hallmarks_mapped <- sort(unique(unlist(ora_to_hallmarks)))
length(ora_hallmarks_mapped)


#Read GSEA and GSCORE Hallmark results
fgsea_gsea <- fread(gsea_file, data.table = FALSE)
fgsea_gscore <- fread(gscore_file, data.table = FALSE)

# Keep significant Hallmarks (padj < 0.05)
gsea_hallmarks_sig <- fgsea_gsea$pathway[which(fgsea_gsea$padj < 0.05)]
gscore_hallmarks_sig <- fgsea_gscore$pathway[which(fgsea_gscore$padj < 0.05)]

# Make sure names align to hallmark_list keys
gsea_hallmarks_sig <- intersect(names(hallmark_list), gsea_hallmarks_sig)
gscore_hallmarks_sig <- intersect(names(hallmark_list), gscore_hallmarks_sig)

#Build three sets at Hallmark level for comparison
set_ORA   <- ora_hallmarks_mapped
set_GSEA  <- sort(unique(gsea_hallmarks_sig))
set_GSCORE<- sort(unique(gscore_hallmarks_sig))

# Save sets for later export
writeLines(set_ORA, "mapped_ORA_hallmarks.txt")
writeLines(set_GSEA, "GSEA_hallmarks_sig.txt")
writeLines(set_GSCORE, "GSCORE_hallmarks_sig.txt")

#Jaccard similarity at Hallmark-level (pairwise)
sets_list <- list(ORA = set_ORA, GSEA = set_GSEA, GSCORE = set_GSCORE)

pairwise_jaccard <- function(A,B){
  if (length(A)==0 && length(B)==0) return(1)
  if (length(union(A,B))==0) return(0)
  length(intersect(A,B)) / length(union(A,B))
}
methods <- names(sets_list)
jaccard_mat <- matrix(0, nrow=length(methods), ncol=length(methods), dimnames=list(methods,methods))
for (i in methods) for (j in methods) jaccard_mat[i,j] <- pairwise_jaccard(sets_list[[i]], sets_list[[j]])
round(jaccard_mat,2)
write.csv(jaccard_mat, "Jaccard_matrix_hallmark_level.csv", row.names = TRUE)

# Spearman rank correlation of pathway scores (Hallmark-level)
all_hallmarks_union <- sort(unique(c(names(hallmark_list))))  # use the full hallmark universe

# ORA: compute min pval among ORA pathways that mapped to each hallmark
ora_minp <- setNames(rep(NA_real_, length(all_hallmarks_union)), all_hallmarks_union)
for (h in all_hallmarks_union) {
  mapped_paths <- names(which(sapply(ora_to_hallmarks, function(x) h %in% x)))
  if (length(mapped_paths) == 0) {
    ora_minp[h] <- NA
  } else {
    combined_ora_df <- dplyr::bind_rows(go_in$df, kegg_in$df)
    matched_pvals <- numeric(0)
    for (pn in mapped_paths) {
      core_name <- sub("^(GO:|KEGG:)", "", pn)
      hit_rows <- which(combined_ora_df$Description == core_name | as.character(combined_ora_df$ID) == core_name)
      if (length(hit_rows)>0) matched_pvals <- c(matched_pvals, combined_ora_df$pvalue[hit_rows])
    }
    if (length(matched_pvals) > 0) ora_minp[h] <- min(matched_pvals, na.rm = TRUE) else ora_minp[h] <- NA
  }
}
score_ORA <- -log10(ora_minp + 1e-300)  # NA remain NA

# GSEA and GSCORE: build score by -log10(padj)
gsea_padj <- setNames(rep(NA_real_, length(all_hallmarks_union)), all_hallmarks_union)
if (nrow(fgsea_gsea) > 0) {
  idx <- match(fgsea_gsea$pathway, all_hallmarks_union)
  gsea_padj[fgsea_gsea$pathway] <- fgsea_gsea$padj
}
score_GSEA <- -log10(gsea_padj + 1e-300)

gscore_padj <- setNames(rep(NA_real_, length(all_hallmarks_union)), all_hallmarks_union)
if (nrow(fgsea_gscore) > 0) {
  gscore_padj[fgsea_gscore$pathway] <- fgsea_gscore$padj
}
score_GSCORE <- -log10(gscore_padj + 1e-300)

scores_df <- data.frame(ORA = score_ORA, GSEA = score_GSEA, GSCORE = score_GSCORE)

#compute pairwise Spearman correlation
spearman_mat <- cor(scores_df, method = "spearman", use = "pairwise.complete.obs")
round(spearman_mat, 2)
write.csv(spearman_mat, "Spearman_matrix_hallmark_scores.csv", row.names = TRUE)



# Jaccard heatmap
png("Figure15_Jaccard_heatmap.png", width = 600, height = 500, res = 150)
pheatmap::pheatmap(jaccard_mat,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("navy","white","firebrick3"))(50),
                   main = "Jaccard Similarity (Hallmark-level)",
                   display_numbers = TRUE)
dev.off()

# Spearman heatmap (ensure finite breaks)
png("Figure16_Spearman_heatmap.png", width = 600, height = 500, res = 150)
pheatmap::pheatmap(spearman_mat,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("blue","white","red"))(50),
                   main = "Spearman Correlation of Hallmark-level Scores",
                   display_numbers = TRUE)
dev.off()

# Print summary numbers
cat("Number of Hallmarks mapped from ORA:", length(set_ORA), "\n")
cat("Number of significant Hallmarks (GSEA):", length(set_GSEA), "\n")
cat("Number of significant Hallmarks (GSCORE):", length(set_GSCORE), "\n")
cat("Jaccard matrix (rounded):\n"); print(round(jaccard_mat,2))
cat("Spearman matrix (rounded):\n"); print(round(spearman_mat,2))

# Save workspace variables for reproducibility
# save(sets_list, jaccard_mat, spearman_mat, scores_df, file = "mapping_Jaccard_Spearman_results.RData")
