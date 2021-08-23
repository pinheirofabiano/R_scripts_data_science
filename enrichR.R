rm(list = ls())

getwd()
setwd("/Users/fabiano/Desktop/projetos_em_andamento/LL37_cancer/data_analysis/CEMiToll_pvalue005/breast/Tables_breast")

breast_modules <- read.delim("module.tsv", sep = "\t")
breast_M7 <- breast_modules[which(breast_modules$modules == "M7"), 1]

install.packages("enrichR")
library(enrichR)

setEnrichrSite("Enrichr")
websiteLive <- TRUE

dbs <- listEnrichrDbs()

dbs <- c("Reactome_2016", "KEGG_2021_Human")

if (websiteLive) {
  enriched <- enrichr(breast_M7, dbs)
}

#Reactome database
if (websiteLive) enriched[["Reactome_2016"]]

if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

#KEGG database

if (websiteLive) enriched[["KEGG_2021_Human"]]

if (websiteLive) plotEnrich(enriched[[2]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

