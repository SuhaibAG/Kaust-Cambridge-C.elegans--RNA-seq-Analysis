## PCA
### After differential expression ###
# Get transformed counts (VST)
vsd <- vst(dds, blind = TRUE)
counts_transformed <- assay(vsd)


# Simple PCA plot using DESeq2's built-in function
plotPCA(vsd, intgroup = c("developmental_stage")) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Paired")

# First, get the PCA data
pca_data <- plotPCA(vsd, intgroup = c("developmental_stage"), returnData = TRUE)

# Then create the plot manually with the same aesthetics
ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 6) +
  geom_text_repel(aes(label = name), size = 4, box.padding = 0.5) +
  xlab(paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100), "% variance")) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "right") +
  scale_color_brewer(palette = "Paired") +
  labs(color = "developmental_stage")
