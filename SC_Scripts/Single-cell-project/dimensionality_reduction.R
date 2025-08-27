library(PCAtools)
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)


gene_var <- modelGeneVar(sce)

gene_var

gene_var %>% 
  as.data.frame() %>% 
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")


hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs)

plotExpression(sce, features = hvgs[1:20], point_alpha = 0.05)

sce <- runPCA(sce, subset_row = hvgs)
percent.var <- attr(reducedDim(sce), "percentVar")
chosen_elbow <- findElbowPoint(percent.var)
chosen_elbow
sce.denoised <- denoisePCA(sce, technical = gene_var, subset.row = hvgs)
ncol(reducedDim(sce.denoised, "PCA"))

set.seed(100)
sce <- runUMAP(sce)
plotUMAP(sce, colour_by = "Sample")


