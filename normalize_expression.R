rm(list = ls())
library("jsonlite")

params <- fromJSON("~/sandbox/working/2020_Sandra/v4/analysis/scripts/config.json")

setwd(params$workdir)
raw.exp <- read.table(file = params$rawcounts, header = T, sep = "\t")
phenodata <- read.table(file = params$phenodata, header = T, sep = "\t")

library(edgeR)
raw.exp$meanG <- apply( raw.exp[ , 2:ncol( raw.exp ) ], 1, mean )
raw.exp <- raw.exp[ order( raw.exp[ , 'Symbol' ], raw.exp[ , 'meanG' ] ), ]
raw.exp <- raw.exp[ !duplicated( raw.exp$Symbol ), ]
rownames(raw.exp) <- raw.exp$Symbol
raw.exp <- raw.exp[ , 2:( ncol( raw.exp ) - 1 ) ]

write.table(x = data.frame(Symbol = rownames(raw.exp), raw.exp), file = "analysis/geneCounts/raw_expression.tsv", quote = F, sep = "\t", row.names = F)

if(params$normalizeMethod == "CPM") {
  cpm.exp <- normalizeQuantiles(cpm(raw.exp))
  write.table(x = data.frame(Symbol = rownames(cpm.exp), cpm.exp), file = "analysis/geneCounts/cpm_expression.tsv", quote = F, sep = "\t", row.names = F)
} else if(params$normalizeMethod == "Log2CPM") {
  cpm.exp <- normalizeQuantiles(log2(cpm(raw.exp)+1))
  write.table(x = data.frame(Symbol = rownames(cpm.exp), cpm.exp), file = "analysis/geneCounts/log2cpm_expression.tsv", quote = F, sep = "\t", row.names = F)
} else if(params$normalizeMethod == "TMM") {
  frame <- DGEList(raw.exp, genes = rownames(raw.exp), group = phenodata$Class)
  norm <- calcNormFactors(frame, method = "TMM")
  tmm.exp <- normalizeQuantiles(cpm(norm, log = F))
  write.table(x = data.frame(Symbol = rownames(tmm.exp), tmm.exp), file = "analysis/geneCounts/tmm_expression.tsv", quote = F, sep = "\t", row.names = F)
} else if(params$normalizeMethod == "Log2TMM") {
  frame <- DGEList(raw.exp, genes = rownames(raw.exp), group = phenodata$Class)
  norm <- calcNormFactors(frame, method = "TMM")
  tmm.exp <- normalizeQuantiles(log2(cpm(norm, log = F)+1))
  write.table(x = data.frame(Symbol = rownames(tmm.exp), tmm.exp), file = "analysis/geneCounts/log2tmm_expression.tsv", quote = F, sep = "\t", row.names = F)
}
