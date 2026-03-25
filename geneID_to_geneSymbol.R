
# convert NCBI gene ID to gene Symbol with the function 
# getGenes from the package mygene

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
BiocManager::install("mygene")
library(mygene)
library(biomaRt)


df <- AD_TPMs
rownames(df) <- df$GeneID
df <- df[,-1]
df <- df[c("GSM7552864","GSM7552867","GSM7552869","GSM7552870","GSM7552874","GSM7552876","GSM7552878","GSM7552880")]
colnames(df) <- c("control01", "control02", "control03", "control04", "alzheimer01", "alzheimer02", "alzheimer03", "alzheimer04")

result <- getGenes(geneid = rownames(df), fields = c("symbol"))

symbols <- result$symbol

rownames(df) = make.names(symbols, unique=TRUE)

# the code below is for ensembl IDs

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_names_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                            filters = "ensembl_gene_id",
                            values = raw_autopsias$row.names,
                            mart = ensembl)





# bind column with the Symbol IDs
df <- cbind(symbols, df)
df <- df[,-2]

#order columns alphabetically
df <- df[,order(colnames(df))]



#save final count matrix
write.table(df, "~/Desktop/AD_TPMs_final.tsv", sep="\t")


