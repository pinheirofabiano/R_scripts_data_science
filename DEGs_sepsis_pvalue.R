rm(list = ls())


source("https://raw.githubusercontent.com/nicolau/code-R/master/DEG_analysis.R")

setwd("/Users/fabiano/Desktop/projetos_em_andamento/INOVA_USP/projeto_AUTOPSIAS/data_analysis/results")


treatedGroupArray <- c("colon.sepsis", "cortex.sepsis", "heart.sepsis", "hippocampus.sepsis", "kidney.sepsis", "lung.sepsis")
controlGroupArray <- c("colon.control", "cortex.control", "heart.control", "hippocampus.control", "kidney.control", "lung.control")

pvalue <- 0.05

dirOut <- "DEGs_pvalue"
if(!dir.exists(dirOut)) {
  dir.create(dirOut, recursive = T)
}


for(i in 1:length(treatedGroupArray)) { ################
  # i <- 1
  
  counts <- read.table(file = "counts.tsv", header = T, sep = "\t")
  phenodata <- read.table(file = "phenodata.tsv", header = T, sep = "\t")
  
  
  treatedGroup <- treatedGroupArray[i] ################
  controlGroup <- controlGroupArray[i] ################
  message(paste0(treatedGroup, " vs ", controlGroup))

  phenodata$Class <- paste0(phenodata$Tissue, ".", phenodata$Class) ################
  
  phenodata <- phenodata[ c( which( phenodata$Class == controlGroup ),
                                 which( phenodata$Class == treatedGroup ) ), ]
  
  phenodata <- phenodata[ phenodata$Sample %in% colnames( counts ), ]
  #head(samplesinfo)
  counts <- counts[ ,c( "Symbol", as.character( phenodata$Sample ) ) ]
  phenodata <- phenodata[ phenodata$Sample %in% colnames( counts ), ]
  
  rownames( counts ) <- as.character( counts$Symbol )
  counts$Symbol <- NULL
  
  message(paste0("Control: ", controlGroup))
  head(phenodata)
  
  results_DEG <- DEG_analysis(data = counts, samplesinfo = phenodata, nontreated = controlGroup,
                              treated = treatedGroup, method = "edgeR", disp = 0, verbose = T)
  
  write.table(x = data.frame(Symbol = rownames(results_DEG), results_DEG), file = paste0(dirOut, "/", treatedGroup, "_vs_", controlGroup, ".tsv"), quote = F, sep = "\t", row.names = F)
  
  DEGs <- results_DEG[results_DEG$PValue < pvalue,]
  write.table(x = data.frame(Symbol = rownames(DEGs), DEGs),
              file = paste0(dirOut, "/", treatedGroup, "_vs_", controlGroup, "_pvalue_", pvalue, "_all_degs.tsv"),
              quote = F, sep = "\t", row.names = F)
  
  DEGs <- results_DEG[results_DEG$PValue < pvalue & results_DEG$logFC > 0,]
  write.table(x = data.frame(Symbol = rownames(DEGs), DEGs),
              file = paste0(dirOut, "/", treatedGroup, "_vs_", controlGroup, "_pvalue_", pvalue, "_up-regulated.tsv"),
              quote = F, sep = "\t", row.names = F)
  
  DEGs <- results_DEG[results_DEG$PValue < pvalue & results_DEG$logFC < 0,]
  write.table(x = data.frame(Symbol = rownames(DEGs), DEGs),
              file = paste0(dirOut, "/", treatedGroup, "_vs_", controlGroup, "_pvalue_", pvalue, "_down-regulated.tsv"),
              quote = F, sep = "\t", row.names = F)
  
  message("")
}
