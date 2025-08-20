library(scater)
library(scran)
library(tidyverse)
library(BiocParallel)
library(DESeq2)

bpp <- MulticoreParam(7)
set.seed(100) 

##clustering
clust <- quickCluster(sce, BPPARAM = bpp)

rowData(sce)
##computing pooled factors
sce <- computePooledFactors(sce,
                            clusters = clust,
                            min.mean = 0.1,
                            BPPARAM = bpp)
assayNames(sce)
sce <- logNormCounts(sce)
counts <- counts(sce)

deconv.sf <- sizeFactors(sce)
summary(deconv.sf)
lib.sf <- librarySizeFactors(sce)

sce$sizeFactor
assay(sce, "logcounts")

