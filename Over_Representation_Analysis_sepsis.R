rm(list = ls())


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")

library(clusterProfiler)

setwd("/Users/fabianopinheiro/Desktop/projetos_em_andamento/INOVA_USP/Projeto_AUTOPSIAS/data_analysis/results")

phenodata <- read.table(file = "phenodata.tsv", header = T, sep = "\t")

install.packages("GSA")
install.packages("reshape2")
install.packages("tidyverse")

library(GSA)
library(reshape2)
library(tidyverse)

geneSet1 <-  GSA.read.gmt("Reactome.gmt")
geneSet2 <-  GSA.read.gmt("KEGG.gmt")

pathways_Reactome <- data.frame(matrix(unlist(geneSet1$geneset.names)),stringsAsFactors=FALSE)
pathways_KEGG <- data.frame(matrix(unlist(geneSet2$geneset.names)),stringsAsFactors=FALSE)

n.obs <- sapply(geneSet1[1], length)
seq.max <- seq_len(max(n.obs))
genes_Reactome <- as.data.frame(sapply(geneSet1[1], "[", i = seq.max))

n.obs <- sapply(geneSet2[1], length)
seq.max <- seq_len(max(n.obs))
genes_KEGG <- as.data.frame(sapply(geneSet2[1], "[", i = seq.max))

term_to_gene_Reactome <- cbind(pathways_Reactome, genes_Reactome)
term_to_gene_KEGG <- cbind(pathways_KEGG, genes_KEGG)


for(direction in c("up", "down")) {
  
    message(direction)
    geneData <- read.table(file = paste0("DEGs_pvalue/lung.sepsis_vs_lung.control_pvalue_0.05_", direction, "-regulated.tsv"), header = T, sep = "\t")
    print(geneData)
    class(geneData)
    geneList <- as.vector(geneData[, "Symbol"])
    class(geneList)
    geneList <- gsub(pattern = "\\|.*", replacement = "", x = geneList)
    print(geneList)

    
   
    df_Reactome      <- enricher (gene = geneList,
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "none",
                       universe      = unique(unlist(geneSet1$genesets)),
                       minGSSize     = 15,
                       maxGSSize     = 500,
                       qvalueCutoff  = 0.02,
                       TERM2GENE     = term_to_gene_Reactome,
                       TERM2NAME = NA
    )

    print(df_Reactome)
    
    df_Reactome <- as.data.frame(df_Reactome@result)
    
    # calculate overlap
    x <- colsplit(df_Reactome[, 3], "\\/", c("a", "b"))
    y <- colsplit(df_Reactome[, 4], "\\/", c("a", "b"))
    
    df_Reactome[, "overlap"] <- round(100 * x[,  1]/y[, 1], 1)
    
    # reorder columns
    colsOrder <- c("ID", "Description", "GeneRatio", "BgRatio",
                   "overlap", "pvalue", "p.adjust", "qvalue",
                   "geneID", "Count")
    
    
    df_Reactome <- df_Reactome[, colsOrder]
    print(df_Reactome)
    df_Reactome_sub <- df_Reactome[df_Reactome$pvalue < 0.05,]
    df_Reactome_sub <- df_Reactome_sub[df_Reactome_sub$Count >= 5,]
    print(df_Reactome_sub)
    
    if(dim(df_Reactome_sub)[1] > 0){
      
      dirOut <- "ORA"
      if(!dir.exists(dirOut)) {
        dir.create(dirOut, recursive = T)
      }
    }
      library(ggplot2)
      ggplot(data = df_Reactome_sub, aes(x = reorder(ID, -log10(pvalue)), y = -log10(pvalue), fill = -log10(pvalue))) +
        geom_bar(stat = "identity") +
        geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") +
        xlab("Pathways") +
        ylab("Pvalue") +
        coord_flip() +
        theme_bw() +
        theme(text = element_text(size = 8))
      # dev.off()
      ggsave(filename = paste0("barplot_Reactome_pvalue_0.05_", direction, ".pdf"), path = "ORA", width = 8, height = 1.7)
      # ggsave(filename = "ORA/barplot_Reactome_pvalue_0.05_.png", width = 6, height = 3)
      

      
        # write output
        outname <- paste0("ORA/clusterProfiler_Reactome_testGMP_pvalue_0.05_", direction, ".tsv")
        # outname <- "ORA/clusterProfiler_Reactome_testGMP_pvalue_0.05.tsv"
        
      write.table(x = df_Reactome, file = outname, row.names = F, col.names = T, quote = F, sep = "\t") 
    

}
      

  for(direction in c("up", "down")) {
  
  message(direction)
  geneData <- read.table(file = paste0("DEGs_pvalue/lung.sepsis_vs_lung.control_pvalue_0.05_", direction, "-regulated.tsv"), header = T, sep = "\t")
  print(geneData)
  class(geneData)
  geneList <- as.vector(geneData[, "Symbol"])
  class(geneList)
  geneList <- gsub(pattern = "\\|.*", replacement = "", x = geneList)
  print(geneList)
  
      df_KEGG      <- enricher(gene = geneList,
                                   pvalueCutoff  = 0.05,
                                   pAdjustMethod = "none",
                                   universe      = unique(unlist(geneSet2$genesets)),
                                   minGSSize     = 15,
                                   maxGSSize     = 500,
                                   qvalueCutoff  = 0.02,
                                   TERM2GENE     = term_to_gene_KEGG,
                                   TERM2NAME = NA)
      
      df_KEGG <- as.data.frame(df_KEGG@result)
# calculate overlap
x <- colsplit(df_KEGG[, 3], "\\/", c("a", "b"))
y <- colsplit(df_KEGG[, 4], "\\/", c("a", "b"))

df_KEGG[, "overlap"] <- round(100 * x[,  1]/y[, 1], 1)

# reorder columns
colsOrder <- c("ID", "Description", "GeneRatio", "BgRatio",
               "overlap", "pvalue", "p.adjust", "qvalue",
               "geneID", "Count")


df_KEGG <- df_KEGG[, colsOrder]
df_KEGG_sub <- df_KEGG[df_KEGG$pvalue < 0.05,]
df_KEGG_sub <- df_KEGG_sub[df_KEGG_sub$Count >= 5,]

if(dim(df_KEGG_sub)[1] > 0){
  
  dirOut <- "ORA"
  if(!dir.exists(dirOut)) {
    dir.create(dirOut, recursive = T)
  }
}
library(ggplot2)
ggplot(data = df_KEGG_sub, aes(x = reorder(ID, -log10(pvalue)), y = -log10(pvalue), fill = -log10(pvalue))) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") +
  xlab("Pathways") +
  ylab("Pvalue") +
  coord_flip() +
  theme_bw() +
  theme(text = element_text(size = 8))
# dev.off()
ggsave(filename = paste0("barplot_KEGG_pvalue_0.05_", direction, ".pdf"), path = "ORA", width = 8, height = 1.7)
# ggsave(filename = "ORA/barplot_KEGG_pvalue_0.05_.png", width = 6, height = 3)



# write output
outname <- paste0("ORA/clusterProfiler_KEGG_testGMP_pvalue_0.05_", direction, ".tsv")
# outname <- "ORA/clusterProfiler_KEGG_testGMP_pvalue_0.05.tsv"


write.table(x = df_KEGG, file = outname, row.names = F, col.names = T, quote = F, sep = "\t") 

  }

