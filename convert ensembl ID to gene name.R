

rm(list = ls())

# https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html

getwd()
setwd("/Users/fabiano/Desktop/projetos_em_andamento/LL37_pancreas/FPKMs")
library(tidyverse)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library("biomaRt")

# to have more information: browseVignettes("biomaRt")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

samples <- c("control1", "control2", "control3", "control4", "cancer1", "cancer2", "cancer3", "cancer4", "cancer5", "cancer6", "cancer7", "cancer8", "cancer9", "cancer10", "cancer11", "cancer12", "cancer13", "cancer14")

for (sample in samples) {
  
pancreas <- read_delim(paste0("pancreas.", sample, ".FPKM.txt"), col_names = F, delim = "\t")

names(pancreas) <- c("ensembl_ID", "FPKM")

# in case it's gencode, this mostly works
#if ensembl, will leave it alone

pancreas$ensembl_ID <- sub("[.][0-9]*","", pancreas$ensembl_ID)


gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = pancreas$ensembl_ID, mart= mart)




names(gene_IDs) <- c("ensembl_ID", "gene_ID")


pancreas <- left_join(pancreas, gene_IDs, by = "ensembl_ID")

write_delim(pancreas, paste0("pancreas.", sample, ".FPKM_with_geneID.tsv"), delim = "\t")

}
