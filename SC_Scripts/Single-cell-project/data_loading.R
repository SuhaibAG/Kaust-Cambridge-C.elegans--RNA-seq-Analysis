library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)


sample_info <- read_csv("~/Desktop/Course_Materials/week6/single_cell/sample_info.csv")
bp.params <- MulticoreParam(workers = 7)

samples_list <- sample_info$sample

list_of_files <- str_c("~/Desktop/Course_Materials/week6/single_cell/preprocessed/scrnaseq/star/",
                         sample_info$sra_run,"_t",sample_info$timepoint,"/",
                         sample_info$sra_run,"_t",sample_info$timepoint,".Solo.out/Gene/filtered") 

names(list_of_files) <- samples_list
sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)

colData(sce) 
rowData(sce)
dim(sce) #46926 72591

table(rowData(sce)$Symbol)
 