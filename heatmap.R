

install.packages("pheatmap")

library(pheatmap)

setwd("/Users/fabianopinheiro/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results")

data <- read.table(file = "log2cpm_with_variance_per_tissue.tsv", header = T, sep = "\t")

pheatmap(mat = data)

