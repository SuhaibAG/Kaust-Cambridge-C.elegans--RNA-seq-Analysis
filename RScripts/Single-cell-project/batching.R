library(scater)
library(scran)
library(batchelor)
library(bluster)
library(pheatmap)
library(magrittr)
sce <- readRDS("/home/participant/Course_Materials/week6/single_cell/sce")
# obtain a batch-corrected SCE object
sce_all_corrected <- quickCorrect(sce, batch = sce$Sample)$corrected

# add the corrected matrix to the original object - to keep it all together
reducedDim(sce, "corrected") <- reducedDim(sce_all_corrected, "corrected")

# visualise both corrected and uncorrected
plotReducedDim(sce, dimred = "corrected", colour_by = "Sample")
