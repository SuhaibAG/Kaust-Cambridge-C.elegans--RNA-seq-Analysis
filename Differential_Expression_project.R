# Check sample matching between counts and metadata
all(colnames(txi$counts) == sample_info$sample)

# Create DESeq2 object and pre-filter genes
dds_raw <- DESeqDataSetFromTximport(txi = txi, 
                                    colData = sample_info, 
                                    design = ~ developmental_stage)

keep <- rowSums(counts(dds_raw)) > 5
dds <- dds_raw[keep, ]

# Run differential expression with LRT test (any change across stages)
dds <- DESeq(dds, test="LRT", reduced=~1)

# Extract results and check dispersion estimates
res <- results(dds, alpha = 0.05)
plotDispEsts(dds) # Validate model fit