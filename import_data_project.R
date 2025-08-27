# load library 
library(tximport)
library(DESeq2)
library(tidyverse)

# load the metadata 
sample_info <- read_csv("sample_info.csv", col_types = 'ccccc')

# Ensure developmental_stage follows CSV order
sample_info$developmental_stage <- factor(sample_info$developmental_stage, 
                                          levels = unique(sample_info$developmental_stage))

# Remove sample "mpfc_700" from sample_info
sample_info <- sample_info %>%
  dplyr::filter(sample != "mpfc_700")

# Read the count data 
files <- file.path("preprocessed/rnaseq/star_salmon/", sample_info$sample, "quant.sf")
names(files) <- sample_info$sample
tx2gene <- read_tsv("references/tx2gene.tsv")

# import Salmon data and summarize to gene-level
txi <- tximport(files = files, type = 'salmon', tx2gene = tx2gene)
str(txi) 
head(txi$counts)
head(txi$abundance)
head(txi$length)

dim(txi$counts)

## save txi for future sessions
### create directory
dir.create('results/r_objects/', recursive = T)
### Save the object
saveRDS(txi, file = 'results/r_objects/txi.rds')


## filter raw counts matrix ##
rawCounts <- round(txi$counts, 0)
rawCounts

dim(rawCounts)

keep <- rowSums(rawCounts) > 5
keep

table(keep, useNA = 'always')

fillCounts <- rawCounts[keep, ] # filtering lowly expressed genes
dim(fillCounts)
summary(fillCounts)

## visualize raw counts
boxplot(fillCounts, min= 'Raw counts', las=3)

## visualize mean vs sd
plot(rowMeans(fillCounts), rowSds(fillCounts),
     main= 'Raw counts: sd vs mean',
     xlin=c(0, 10000),
     ylin= c(0,5000))

# log2 transform data 
logcount <- log2(fillCounts +1)
plot(rowMeans(logcount), rowSds(logcount),
     main= 'Log2Counts: sd vs mean')

rlogcounts <- rlog(fillCounts)
boxplot(rlogcounts, main= 'Rlog counts', las=2)

plot(rowMeans(rlogcounts), rowSds(rlogcounts),
     main= 'Rlog counts: sd vs mean')
