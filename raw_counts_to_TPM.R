getwd()
install.packages("tidyverse")
install.packages("tidyr")
library(tidyverse)
library(tidyr)
BiocManager::install(version = "3.22")
# Install and load the biomaRt package
if(!requireNamespace("biomaRt", quietly = TRUE)) install.packages("biomaRt")
library(biomaRt)

BiocManager::install("mygene")
library(mygene)

df <- metanalysis_braga_healthy_only_t1_final

result <- getGenes(geneid = df$GeneID, fields = c("symbol"))

symbols <- result$symbol

# bind column with the Symbol IDs
df <- cbind(symbols, df)

# remove rows with NAs in the column "symbols"

df <- df |>
  drop_na(symbols)

# convert gene symbol to ensembl ID

# Example for biomaRt package


mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

annotLookup <- getBM(mart = mart, attributes = c('hgnc_symbol', 'ensembl_gene_id'),
  values = df$symbols, uniqueRows = TRUE)

colnames(annotLookup) <- c("symbols", "ensembl_ID")

colnames(df)[colnames(df) == 'row.names'] <- 'symbols'
df_final <- inner_join(annotLookup, df)

# search gene transcripts length

ensembl_list <- df_final$ensembl_ID

gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id", "transcript_length"), filters="ensembl_gene_id", 
                  values=ensembl_list, mart=mart)

# keep only maximum transcript lengths

colnames(gene_coords) <- c("symbols", "ensembl_ID", "transcript_length")

gene_coords <- gene_coords %>% 
  group_by(symbols) %>% 
  slice_max(transcript_length)

df_final <- inner_join(gene_coords, df_final)


countsMatrix <- df_final[,5:81]

class(countsMatrix)
countsMatrix <- as.matrix(countsMatrix)

install.packages("DGEobj.utils")
library(DGEobj.utils)

pneumonia_TPMs <- convertCounts(
  countsMatrix = countsMatrix,
  unit = "tpm",
  geneLength = df_final$transcript_length,
  log = FALSE,
  normalize = "none",
  prior.count = NULL
)

pneumonia_TPMs <- as.data.frame(pneumonia_TPMs)
rownames(pneumonia_TPMs) = make.names(df_final$symbols, unique=TRUE)

#save normalized table (TPMs)
write.table(pneumonia_TPMs, "~/Desktop/zhang_pneumonia_TPMs.tsv", sep="\t")


