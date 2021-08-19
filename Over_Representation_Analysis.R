rm(list = ls())
library("jsonlite")

params <- fromJSON("~/sandbox/working/2020_Sandra/v4/analysis/scripts/config.json")

phenodata <- read.table(file = "analysis/geneCounts/phenodata.tsv", header = T, sep = "\t")
unique_classes <- unique(sort(phenodata$Class))
unique_classes <- unique_classes[unique_classes != params$control_label]

setwd(params$workdir)

# libraries
suppressMessages(suppressWarnings({
  library("clusterProfiler")
  library("GSA")
  library("reshape2")
}))

for(pathway_file in params$ORA$databases){
  geneSets   <- read.gmt(pathway_file)
  
  
  for(comp in unique_classes) {
    comp <- "La"
    message(comp)
    for(direction in c("up", "down")) {
      direction <- "up"
      message(direction)
      geneData <- read.table(file = paste0("analysis/DEG/", params$DEG$method, "/", comp, "_vs_Uninfected_FDR_", params$DEG$FDRcut, "_", direction, "-regulated.tsv"), header = T, sep = "\t")
      geneList <- as.vector(geneData[, "Symbol"])
      geneList <- gsub(pattern = "\\|.*", replacement = "", x = geneList)
      
      df      <- enricher(geneList,
                          minGSSize     = 15,
                          maxGSSize     = 500,
                          universe      = unique(unlist(geneSets$gene)),
                          pvalueCutoff  = params$ORA$FDRcut,
                          qvalueCutoff  = 0.02,
                          pAdjustMethod = "fdr",
                          TERM2GENE     = geneSets
      )
      
      #names(df@geneSets)
      df <- as.data.frame(df@result)
      
      # calculate overlap
      x <- colsplit(df[, 3], "\\/", c("a", "b"))
      y <- colsplit(df[, 4], "\\/", c("a", "b"))
      df <- as.data.frame(df)
      df[, "overlap"] <- round(100 * x[,  1]/y[, 1], 1)
      
      # reorder columns
      colsOrder <- c("ID", "Description", "GeneRatio", "BgRatio",
                     "overlap", "pvalue", "p.adjust", "qvalue",
                     "geneID", "Count")
      
      
      df <- df[, colsOrder]
      df_sub <- df[df$pvalue < params$ORA$FDRcut,]
      df_sub <- df_sub[df_sub$Count >= params$ORA$minOverlap,]
      
      if(dim(df_sub)[1] > 0){
        
        dirOut <- paste0(params$ORA$output_dir, "/", params$ORA$method)
        if(!dir.exists(dirOut)) {
          dir.create(dirOut, recursive = T)
        }
        
        library(ggplot2)
        ggplot(data = df_sub, aes(x = reorder(ID, -log10(pvalue)), y = -log10(pvalue), fill = -log10(pvalue))) +
          geom_bar(stat = "identity") +
          geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "dashed") +
          xlab("Pathways") +
          ylab("Pvalue") +
          coord_flip() +
          theme_bw() +
          theme(text = element_text(size = 8))
        # dev.off()
        # ggsave(filename = paste0(dirResult, "/barplot_pvalue_0.05_", comp, "_", direction, ".pdf"), width = 8, height = 1.7)
        ggsave(filename = paste0(dirOut, "/barplot_pvalue_0.05_.png"), width = 6, height = 3)
        
        
        
        # write output
        # outname <- paste0(dirResult, "/clusterProfiler_testGMP_pvalue_0.05_", comp, "_", direction, ".tsv")
        outname <- paste0(dirOut, "/clusterProfiler_testGMP_pvalue_0.05.tsv")
        write.table(x = df, file = outname, row.names = F, col.names = T, quote = F, sep = "\t")
        
        
        
        
        
        pathways_vs_genes_all <- NULL
        # colnames(pathways_vs_genes_all) <- c("Pathway", geneList)
        for(pathway in df_sub$ID) {
          # pathway <- df_sub$ID[1]
          message(pathway)
          genes <- strsplit(x = df_sub[which(df_sub$ID == pathway), "geneID"], split = "\\/", perl = TRUE)[[1]]
          DF_temp <- data.frame(Pathway = pathway, symbol = genes, value = 1, Direction = direction)
          
          if(is.null(pathways_vs_genes_all)) {
            pathways_vs_genes_all <- DF_temp
          } else {
            pathways_vs_genes_all <- rbind(pathways_vs_genes_all, DF_temp)
          }
          # genes <- geneSets[geneSets$ont == pathway, "gene"]
          
          norm.exp <- read.table(file = paste0("analysis/geneCounts/", stringi::stri_trans_tolower(params$normalizeMethod), "_expression.tsv"), header = T, sep = "\t")
          
          # all_genes_edgeR_sub <- all_genes_edgeR[all_genes_edgeR$Symbol %in% genes,]
          norm.exp_sub <- norm.exp[gsub("\\|.*", "", norm.exp$Symbol) %in% genes,]
          # exp.cpm.quantile_sub <- exp.cpm.quantile_sub[exp.cpm.quantile_sub$Symbol %in% geneList,]
          
          # rownames(all_genes_edgeR_sub) <- all_genes_edgeR_sub$Symbol
          # all_genes_edgeR_sub$Symbol <- NULL
          rownames(norm.exp_sub) <- norm.exp_sub$Symbol
          norm.exp_sub$Symbol <- NULL
          norm.exp_sub <- as.data.frame(t(scale(t(norm.exp_sub))))
          
          # samplesinfo_sub <- samplesinfo[match(as.vector(samplesinfo$Sample), colnames(exp.cpm.quantile_sub)),]
          phenodata_sub <- phenodata[order(phenodata$Class),]
          rownames(phenodata_sub) <- phenodata_sub$Sample
          phenodata_sub$Sample <- NULL
          
          # samplesinfo_sub <- samplesinfo[match(as.vector(samplesinfo$Sample), colnames(exp.cpm.quantile_sub)),]
          
          # dim(exp.cpm.quantile_sub)
          # dim(samplesinfo)
          
          pathway <- gsub(pattern = "\\/", replacement = " ", x = pathway)
          norm.exp_sub <- norm.exp_sub[, match(rownames(phenodata_sub), colnames(norm.exp_sub))]
          # pheatmap::pheatmap(mat = as.matrix(all_genes_edgeR_sub), cluster_cols = F, cluster_rows = F,
          #                    breaks = seq(-3, 3, length.out = 101)
          #                    , filename = paste0(dirResult, "/heatmap_of_", pathway, "_", comp, "_", direction, "group.png")
          # )
          pheatmap::pheatmap(mat = norm.exp_sub,
                             annotation_col = phenodata_sub, cluster_cols = F, cluster_rows = F,
                             breaks = seq(-3, 3, length.out = 101)
                             , filename = paste0(dirOut, "/heatmap_of_", pathway, "_", comp, "_", direction, ".png")
          )
        }
        
        library(reshape2)
        pathways_vs_genes_all_matrix <- dcast(pathways_vs_genes_all, Pathway ~ symbol, value.var = "value")
        pathways_vs_genes_all_matrix[is.na(pathways_vs_genes_all_matrix)] <- 0
        rownames(pathways_vs_genes_all_matrix) <- pathways_vs_genes_all_matrix$Pathway
        pathways_vs_genes_all_matrix$Pathway <- NULL
        
        overlap_pathways <- crossprod(as.matrix(t(pathways_vs_genes_all_matrix)))
        overlap_pathways_prop <- floor(t(overlap_pathways * 100 / diag(overlap_pathways)))
        # pheatmap::pheatmap(mat = overlap_pathways_prop, cluster_rows = F, cluster_cols = F)
        
        pheatmap::pheatmap(mat = pathways_vs_genes_all_matrix, color = c("gray90", "gray20"),
                           filename = paste0(dirOut, "/heatmap_of_overlap_genes_between_pathways", comp, "_", direction, ".png"),
                           width = 15, height = 5)
        # dev.off()
      }
    }
  }
}
