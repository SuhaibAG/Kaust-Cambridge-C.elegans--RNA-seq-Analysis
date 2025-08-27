library(pheatmap)
library(dendextend)
library(ggvenn)
library(RColorBrewer)
#### Annotation ####

# Open AnnotationHub
ah <- AnnotationHub()

# Query for C. elegans EnsDb (replace 'Caenorhabditis elegans' with the correct string if needed)
ce_ensdb <- query(ah, c("EnsDb", "Caenorhabditis elegans"))[[1]]

# Extract gene info as a data frame
annot <- genes(ce_ensdb, return.type = "data.frame")

# Keep only relevant columns
annot <- annot %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype)

# Filter to genes that are in DESeq2 results
annot <- annot %>%
  dplyr::filter(gene_id %in% rownames(res))

# (Optional) handle duplicate entrez IDs like in course
dup_entrez <- annot %>%
  dplyr::filter(!is.na(entrezid)) %>%
  add_count(entrezid) %>%
  arrange(entrezid) %>%
  dplyr::filter(n > 1)

# Annotate results
res_annotated <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  left_join(annot, by = "gene_id")

# Check
head(res_annotated)


#### visualization ####
# check by histogram
hist(res_annotated$pvalue)

# Shrink the result 
res_shrink <- lfcShrink(dds, 
                            res = res,
                            type = "ashr")

# Re annotate the shrinked result 
res_shrink_annot <- as.data.frame(res_shrink) %>%
  rownames_to_column("gene_id") %>% 
  left_join(annot, "gene_id")

# MA plot
par(mfrow = c(1,2))
plotMA(res, alpha = 0.05)
plotMA(res_shrink, alpha = 0.05)


# Venn diagrame 
get_genes <- function(res_shrink, direction) {
  sign <- ifelse(direction == "up", 1, -1)
  res_shrink %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::filter(sign * log2FoldChange > 0) %>%
    pull("gene_id")
}
vennList <- list(Upregulated_res = get_genes(res_shrink_annot, "up"),
                 Downregulated_res = get_genes(res_shrink_annot, "down"))
str(vennList)

ggvenn(vennList, set_name_size = 4)


