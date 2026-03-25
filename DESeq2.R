
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")

rownames(GSE310929_shock_metadata) <- GSE310929_shock_metadata$Sample.ID

GSE310929_shock_metadata_only_t1 <- GSE310929_shock_metadata_only_t1[-grep("_t3",GSE310929_shock_metadata_only_t1$Sample.ID),]

colnames(GSE310929_shock_metadata_only_t1)[colnames(GSE310929_shock_metadata_only_t1) == "Disease.Simplified"] <- "condition"

dds <- DESeqDataSetFromMatrix(countData = metanalysis_braga_healthy_only_t1_final, colData = GSE310929_shock_metadata_only_t1, design= ~condition)

str(dds)

# removal of low count reads (optional)

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

# setting reference for DEG analysis

dds$condition <- relevel(dds$condition, ref = "Healthy")

deg <- DESeq(dds)

res <- results(deg)
