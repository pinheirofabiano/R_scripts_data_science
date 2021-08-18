

#BOXPLOTS

getwd()
setwd("/Users/fabiano")
install.packages("tidyverse")
library(tidyverse)
cytokines_kidneys <- read_delim("Desktop/projetos_em_andamento/lesao_renal/Frontiers_revised/NGAL.txt", delim = "\t")

str(ileum)
pairwise.wilcox.test(cytokines_kidneys$IL6, cytokines_kidneys$GROUP, paired = F, p.adj = "bonferroni", na.rm = T)
cytokines_kidneys$TNF_alpha <- as.numeric(cytokines_kidneys$TNF_alpha)
class(cytokines_kidneys)
str(cytokines_kidneys)
names(cytokines_kidneys) <- c("GROUP", "TNF_alpha", "IL1beta", "IL6")

bia_summary <- beatriz %>% 
  group_by(GROUP) %>% 
  summarise(median_potassium = median(potassium), median_calcium = median(calcium))

# put limits using ylim(0,1500); coord_cartesian()
# facet_wrap(~cor): subdivide o gráfico pela cor
# labs(x= "", title = ""): muda o nome das variáveis dos eixos x ou y e cria título
# coord_flip(): inverte eixo x pelo eixo y
# fct_reorder() para ordenar os nomes
# fill e color para colorir o gráfico

cytokines_kidneys$GROUP<-factor(cytokines_kidneys$GROUP, levels=c("KO sepsis", "WT sepsis", "KO rhabdo", "WT rhabdo", "KO control", "WT control"))

ggplot(cytokines_kidneys) + geom_boxplot(aes(x = GROUP, y = IL6), outlier.shape = NA) + 
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_text(color = "black", size = 16, angle = 0), 
        axis.title.x = element_text(color = "black", size =22, angle = 0),
        axis.text.x = element_text(color = "black", size = 22, angle = 0)) + ylab("IL-6") + ylim(0,30) + coord_flip() 










