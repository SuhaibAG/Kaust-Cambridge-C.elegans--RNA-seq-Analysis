library(clusterProfiler)
library(org.Ce.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(dplyr)

# Get unique developmental stages from your data
# Assuming you have a dds object with colData containing developmental_stage
developmental_stages <- unique(colData(dds)$developmental_stage)
developmental_stages <- sort(developmental_stages)  # Sort chronologically if needed

# 3. For each developmental stage, identify differentially expressed genes
go_results <- list()

for (stage in developmental_stages) {
  # Get samples for this developmental stage
  stage_samples <- rownames(colData(dds))[colData(dds)$developmental_stage == stage]
  other_samples <- setdiff(colnames(plot_data), stage_samples)
  
  cat(paste0("\nProcessing ", stage, " with ", length(stage_samples), " samples...\n"))
  
  # Calculate mean expression and fold change
  stage_mean <- rowMeans(plot_data[, stage_samples, drop = FALSE])
  other_mean <- rowMeans(plot_data[, other_samples, drop = FALSE])
  fold_change <- stage_mean - other_mean
  
  # Select top upregulated genes in this developmental stage
  top_genes <- names(sort(fold_change, decreasing = TRUE))[1:51]
  
  # Convert gene names to Entrez IDs
  gene_to_entrez <- res_shrink_annot %>%
    dplyr::filter(gene_name %in% top_genes) %>%
    dplyr::select(gene_name, entrezid) %>%
    distinct()
  
  entrez_ids <- gene_to_entrez$entrezid[!is.na(gene_to_entrez$entrezid)]
  
  if (length(entrez_ids) > 5) {
    cat(paste0("Analyzing ", length(entrez_ids), " genes for ", stage, "\n"))
    
    # Biological Process
    go_bp <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Ce.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    # Cellular Component
    go_cc <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Ce.eg.db,
      keyType = "ENTREZID",
      ont = "CC",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    # Molecular Function
    go_mf <- enrichGO(
      gene = entrez_ids,
      OrgDb = org.Ce.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    
    go_results[[stage]] <- list(
      BP = go_bp,
      CC = go_cc,
      MF = go_mf,
      gene_count = length(entrez_ids),
      samples = stage_samples
    )
  } else {
    cat(paste0("Not enough genes for enrichment in ", stage, "\n"))
  }
}

# 4. Visualize results for each developmental stage
for (stage in developmental_stages) {
  results <- go_results[[stage]]
  
  if (!is.null(results)) {
    cat(paste0("\n=== ", stage, " (", results$gene_count, " genes, ", 
               length(results$samples), " samples) ===\n"))
    cat("Samples:", paste(results$samples, collapse = ", "), "\n")
    
    # Plot top BP terms
    if (nrow(results$BP) > 0) {
      print(dotplot(results$BP, showCategory=15, title=paste(stage, "- Biological Process")))
    } else {
      cat("No significant BP terms found\n")
    }
    
    # Plot top CC terms
    if (nrow(results$CC) > 0) {
      print(dotplot(results$CC, showCategory=15, title=paste(stage, "- Cellular Component")))
    } else {
      cat("No significant CC terms found\n")
    }
    
    # Plot top MF terms
    if (nrow(results$MF) > 0) {
      print(dotplot(results$MF, showCategory=15, title=paste(stage, "- Molecular Function")))
    } else {
      cat("No significant MF terms found\n")
    }
  }
}

# 5. Print developmental stage information
cat("\n=== Developmental Stage Information ===\n")
for (stage in developmental_stages) {
  samples <- rownames(colData(dds))[colData(dds)$developmental_stage == stage]
  
  cat(paste0("\n", stage, ":\n"))
  cat("  Samples (", length(samples), "): ", paste(samples, collapse = ", "), "\n")
}

# 6. Compare developmental stages
all_entrez_ids <- lapply(go_results, function(x) {
  if (!is.null(x)) unique(unlist(lapply(c(x$BP, x$CC, x$MF), function(y) y@gene)))
})

# Compare developmental stages for BP
compare_bp <- compareCluster(
  all_entrez_ids,
  fun = "enrichGO",
  OrgDb = org.Ce.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05
)

if (nrow(compare_bp) > 0) {
  print(dotplot(compare_bp, showCategory=4, title="GO BP Comparison Between Developmental Stages", font.size = 10, group = F))
}

# Compare developmental stages for MF
compare_MF <- compareCluster(
  all_entrez_ids,
  fun = "enrichGO",
  OrgDb = org.Ce.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pvalueCutoff = 0.05
)

if (nrow(compare_MF) > 0) {
  print(dotplot(compare_MF, showCategory=5, title="GO MF Comparison Between Developmental Stages", font.size = 9, group = F))
}

# Compare developmental stages for CC
compare_CC <- compareCluster(
  all_entrez_ids,
  fun = "enrichGO",
  OrgDb = org.Ce.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 0.05
)

if (nrow(compare_CC) > 0) {
  print(dotplot(compare_CC, showCategory=5, title="GO CC Comparison Between Developmental Stages", font.size = 9, group = F))
}

# 7. Show top terms summary
cat("\n=== Top GO Terms Summary ===\n")
for (stage in developmental_stages) {
  results <- go_results[[stage]]
  
  if (!is.null(results) && nrow(results$BP) > 0) {
    cat(paste0("\n", stage, " - Top 5 BP terms:\n"))
    print(head(results$BP$Description, 5))
  }
  
  if (!is.null(results) && nrow(results$MF) > 0) {
    cat(paste0("\n", stage, " - Top 5 MF terms:\n"))
    print(head(results$MF$Description, 5))
  }
  
  if (!is.null(results) && nrow(results$CC) > 0) {
    cat(paste0("\n", stage, " - Top 5 CC terms:\n"))
    print(head(results$CC$Description, 5))
  }
}

