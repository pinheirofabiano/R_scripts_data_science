

rm(list = ls())

library(grid)
library(gridExtra)

getwd()

fgsea <- read.delim("Desktop/projetos_em_andamento/Projeto_Autopsias/article/final_format/fgsea_results_all_pathways_pval0.01.tsv", header = T, sep = "\t")

maxrow <-  35
npages <-  ceiling(nrow(fgsea)/maxrow)
pdf("Desktop/supplementary_table3.pdf", height=11, width=10)

idx <- seq(1, maxrow)

grid.table(fgsea[idx,], rows = NULL)

for (i in 2:npages){
  grid.newpage()
  if(i*maxrow <= nrow(fgsea)) {
    idx <-  seq(1 + ((i - 1)*maxrow), i*maxrow)
    }
  else {
    idx <-  seq(1 + ((i - 1)*maxrow), nrow(fgsea))
    }
grid.table(fgsea[idx,], rows = NULL)
}

dev.off()
