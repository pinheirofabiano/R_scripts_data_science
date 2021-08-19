
which(pancreas_final$X1 == "DEFB1")

View(pancreas_final[14491, ])


# inicio
pancreas_final <- read_delim("Desktop/alunos/Lucilene/datasets/pancreas_final.txt", delim = "\t")

pancreas_final_t <- t(pancreas_final[,-1])
pancreas_final_t <- as.data.frame(pancreas_final_t)
names(pancreas_final_t) <- pancreas_final$X1
# pancreas_final_t
class(pancreas_final_t)
View(pancreas_final_t)


vec_cor <- cor(pancreas_final_t$DEFB1, pancreas_final_t)
df_cor <- data_frame(gene = attributes(vec_cor)$dimnames[[2]], cor = c(vec_cor))
str(df_cor)

library(tidyverse)
df_cor %>%
  arrange(cor)

df_cor %>%
  arrange(desc(cor)) %>% 
  head(n = 15)

lungs_final <- read_delim("Desktop/alunos/Lucilene/datasets/lungs_final.txt", delim = ",")

lungs_final_t <- t(lungs_final[,-1])
lungs_final_t <- as.data.frame(lungs_final_t)
names(lungs_final_t) <- lungs_final$gene_name
# lungs_final_t
class(lungs_final_t)
View(lungs_final_t)


vec_cor <- cor(lungs_final_t$CAMP, lungs_final_t)
df_cor <- data_frame(gene = attributes(vec_cor)$dimnames[[2]], cor = c(vec_cor))
str(df_cor)

library(tidyverse)
df_cor %>%
  arrange(cor)

df_cor %>%
  arrange(desc(cor)) %>% 
  head(n = 15)


