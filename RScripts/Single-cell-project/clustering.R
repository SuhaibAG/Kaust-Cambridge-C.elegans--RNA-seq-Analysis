library(scater) 
library(scran)
library(PCAtools)
library(tidyverse)
library(patchwork)
library(igraph)
library(bluster)


out <- clusterSweep(reducedDim(sce, "PCA"),
                    BLUSPARAM = NNGraphParam(),
                    k = as.integer(c(5, 10, 15, 20, 25, 30)),
                    cluster.fun = "louvain",
                    BPPARAM=BiocParallel::MulticoreParam(7))


df <- as.data.frame(out$parameters)

# get the number of clusters
df$num.clusters <- apply(out$clusters, 2, max)

getMeanSil <- function(cluster) {
  sil <- approxSilhouette(reducedDim(sce, "UMAP"), cluster)
  mean(sil$width)
}


df$silhouette <- map_dbl(as.list(out$clusters), getMeanSil)

nclPlot <- ggplot(df, aes(x = k, y = num.clusters)) + 
  geom_line(lwd=2)
silPlot <- ggplot(df, aes(x = k, y = silhouette)) + 
  geom_line(lwd=2)
nclPlot + silPlot

colData(sce) <- cbind(colData(sce), DataFrame(out$clusters))
colLabels(sce) <- sce$k.25_cluster.fun.louvain

plotReducedDim(sce, 
               dimred = "UMAP",
               colour_by = "label", 
               text_by = "label") +
  ggtitle("Leiden k=25 clusters")


plotReducedDim(sce, 
               dimred = "UMAP",
               by_exprs_values = "logcounts",
               colour_by = "timepoint",
               text_by = "timepoint")


saveRDS(sce, "/home/participant/Course_Materials/week6/single_cell/sce-corrected-clustered")
