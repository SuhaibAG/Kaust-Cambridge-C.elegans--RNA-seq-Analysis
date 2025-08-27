# Extract significant genes (padj < 0.01)
sig_genes_df <- res_shrink_annot %>%
  drop_na(entrezid, padj) %>% 
  dplyr::filter(padj < 0.01) %>% 
  dplyr::select(gene_id, gene_name)

plot_data <- assay(vst(dds)[sig_genes_df$gene_id, ])
rownames(plot_data) <- sig_genes_df$gene_name

# Define your exact column order
custom_order <- c(
  "1_cell", "mpfc_010", "2_cell", "mpfc_050", "mpfc_100", "mpfc_150",
  "mpfc_200", "mpfc_220", "mpfc_260", "mpfc_300",
  "mpfc_340", "mpfc_400", "mpfc_410", "mpfc_450",
  "mpfc_490", "mpfc_510", "mpfc_620",
  "mpfc_750", "mpfc_800", "mpfc_810", "mpfc_830"
)

# Create and reorder the clustering
# First create hclust object for sample
sample_hclust <- plot_data %>% 
  t() %>% 
  scale() %>% 
  dist() %>% 
  hclust()

# Then convert to dendrogram and reorder
sample_dend <- as.dendrogram(sample_hclust)
reordered_dend_sample <- rotate(sample_dend, order = custom_order)

# Convert back to hclust for pheatmap
reordered_hclust_sample <- as.hclust(reordered_dend_sample)

plot(sample_hclust)
plot(reordered_hclust_sample)

# modified the hierarchical cluster based on the row (genes)
gene_hclust <- plot_data %>% 
  t() %>%
  scale() %>%
  t() %>%
  dist() %>% 
  hclust()

plot(gene_hclust)

# Prepare annotation data
existing_samples <- intersect(custom_order, colnames(plot_data))

annot_df <- data.frame(
  developmental_stage = colData(dds)[existing_samples, "developmental_stage"],
  row.names = existing_samples
)

# Plot final heatmap
pheatmap(
  plot_data[, existing_samples],
  cluster_cols = reordered_hclust_sample,  # Use the reordered hclust object
  cluster_rows = gene_hclust,
  scale = "row",
  fontsize_row = 15,
  fontsize_col = 15,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  cutree_rows = 7,
  cutree_cols = 5,
  treeheight_row = F,
  annotation_col = annot_df,
  main = "Significant genes (padj < 0.01)",
  angle_col = 45,
    )