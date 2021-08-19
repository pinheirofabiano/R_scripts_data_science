
rm(list = ls())

# browseVignettes("fgsea")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")

install.packages("GSA")
install.packages("data.table")
install.packages("gtools")

library(fgsea)
library(data.table)
library(GSA)
library(gtools)

getwd()
setwd("/Users/fabiano/Desktop/projetos_em_andamento/LL37_cancer/data_analysis/fgsea")

#read files

breast_module <- fread("breast_M7.tsv", header = TRUE)
lung_module <- fread("lung_M1.tsv", header = TRUE)



# rank <- data.frame(Symbol = DEG$Symbol, Log2FC = DEG$metalog2fc)
# rank <- DEG_rank[order(DEG_rank$Log2FC, decreasing = T),]

# write.table(x = rank, file = "rank.rnk", sep = "\t", row.names = F)

# fGSEA_result_Reactome <- run_fgsea(pathwaysDatabase = "Reactome.gmt", ranksOfGenes = "GSEA_simple/colon/DEG_rank.rnk")

# fGSEA_result_Reactome <- fGSEA_result_Reactome[order(fGSEA_result_Reactome$pval, decreasing = F, na.last = NA),]
# fGSEA_result_Reactome <- fGSEA_result_Reactome[fGSEA_result_Reactome$pval <= 0.05,]



# fwrite(x = fGSEA_result_Reactome, file = ("GSEA_simple/colon/Reactome_colon.tsv"), sep = "\t", sep2 = c("", " ", ""))
# gmt.file <- system.file("extdata", "human.reactome.gmt", package="fgsea")



# and running fgsea:

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)

head(fgseaRes)

#-------------------








