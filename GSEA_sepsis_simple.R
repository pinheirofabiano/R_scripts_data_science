rm(list = ls())

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

setwd("/Users/fabianopinheiro/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results")

source("https://raw.githubusercontent.com/nicolau/code-R/master/run_fgsea_function.R")

if(!dir.exists("GSEA_simple/colon")) {
   dir.create("GSEA_simple/colon", recursive = T)
}

DEG <- fread("Meta_Volcano/Metacombs/DEGs_colon.tsv", header = TRUE)
DEG$Symbol <- gsub("\\|.*", "", DEG$gene_ID)
DEG$metalog2fc <- foldchange2logratio(DEG$metafc, base = 2)

DEG_rank <- data.frame(Symbol = DEG$Symbol, Log2FC = DEG$metalog2fc)
DEG_rank <- DEG_rank[order(DEG_rank$Log2FC, decreasing = T),]

write.table(x = DEG_rank, file = "GSEA_simple/colon/DEG_rank.rnk", sep = "\t", row.names = F)

fGSEA_result_Reactome <- run_fgsea(pathwaysDatabase = "Reactome.gmt", ranksOfGenes = "GSEA_simple/colon/DEG_rank.rnk")
fGSEA_result_Reactome <- fGSEA_result_Reactome[order(fGSEA_result_Reactome$pval, decreasing = F, na.last = NA),]
fGSEA_result_Reactome <- fGSEA_result_Reactome[fGSEA_result_Reactome$pval <= 0.05,]



fwrite(x = fGSEA_result_Reactome, file = ("GSEA_simple/colon/Reactome_colon.tsv"), sep = "\t", sep2 = c("", " ", ""))

fGSEA_result_KEGG <- run_fgsea(pathwaysDatabase = "KEGG.gmt", ranksOfGenes = "GSEA_simple/colon/DEG_rank.rnk")
fGSEA_result_KEGG <- fGSEA_result_KEGG[order(fGSEA_result_KEGG$pval, decreasing = F, na.last = NA),]
fGSEA_result_KEGG <- fGSEA_result_KEGG[fGSEA_result_KEGG$pval <= 0.05,]


fwrite(x = fGSEA_result_KEGG, file = ("GSEA_simple/colon/KEGG_colon.tsv"), sep = "\t", sep2 = c("", " ", ""))
