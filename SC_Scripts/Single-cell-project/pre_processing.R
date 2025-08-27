colData(sce) %>%
  as.data.frame() %>% 
  select(Sample) %>% 
  distinct()

##adding row names to the coldata using the barcode 
sample_info$Sample <- sample_info$sample
sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), sample_info, by="Sample", sort=FALSE)
colData(sce) <- subset(colData(sce), select = c("Sample", "Barcode", "sra_run", "timepoint"))
rownames(colData(sce)) <- sce$Barcode

##removing undetected genes
detected_genes <- rowSums(counts(sce)) > 0
sce <- sce[detected_genes,]
dim(sce) #21547 72591


 ##Annotating
ah <- AnnotationHub()
ens.hs.113<- query(ah, c("Caenorhabditis elegans", "EnsDb"))[[1]] 

genes <- rowData(sce)$ID
gene_annot <- AnnotationDbi::select(ens.hs.113, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME")) %>%
              set_names(c("ID", "Chromosome"))

rowData(sce) <- merge(rowData(sce), gene_annot, by = "ID", sort=FALSE)
rownames(rowData(sce)) <- rowData(sce)$ID


##removign mitochondria cells
table(rowData(sce)$Chromosome)
is.mito <- which(rowData(sce)$Chromosome=="MtDNA")
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)
dim(sce)
sce$sra_run
cell_qc_filters <- quickPerCellQC(colData(sce),
                                  sub.fields = TRUE,
                                  batch=sce$Sample)

as.data.frame(cell_qc_filters) %>% summarise(across(everything(), sum))

colData(sce) <- cbind(colData(sce), cell_qc_filters)
dim(sce) #21547 72591
sce <- sce[, !sce$discard]
dim(sce) #21547 70055

