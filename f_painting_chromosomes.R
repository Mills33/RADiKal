setwd("/Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/matrix_files/First_tryAUG/50_obs_remov/")
library(mltools)
library(proxy)
library(textshape)
library(tidyr)
library(viridis)
df <- read.csv("All_matrix.csv", header = F)
df[, "binme"] <- bin_data(df$V7, bins=10, binType = "quantile")

bins <- c(1.72104808211956, 97.998528019456,104.317529345818,108.737395512572, 112.475550271692,
          115.962729433119, 119.512533005112,123.449863636102, 128.076652232464,134.929144600164, 
          197.308178866844)
tags <- c(0, 1, 2, 3, 4,5,6,7,8,9)

equal_bin <- cut(df$V7, 
                 breaks= bins,  
                 include.lowest=TRUE, 
                 right=FALSE,
                 labels = tags)
equal_bin <- as.data.frame(equal_bin)
df <- bind_cols(df, equal_bin)

chrom_1 <- subset(df, V9 == '1' )
chrom_1 <- chrom_1 %>% select(V1, V8, equal_bin)
tmp <- df %>% select(V1, V8, equal_bin)


chr1_mat <- pivot_wider(chrom_1, names_from = V8, values_from = equal_bin)
all_mat <- pivot_wider(tmp, names_from = V8, values_from = equal_bin)
all_mat$V1 <- NULL
chr1_mat$V1 <- NULL

re_ordered_matrix_all <- cluster_matrix(all_mat, dim = 'col', method ="complete")
col.order_all <- names(re_ordered_matrix_all)
re_ordered_matrix_1 <- cluster_matrix(chr1_mat, dim = 'col', method ="complete")
col.order_1 <- names(re_ordered_matrix_1)

df$ID <- factor(df$V8,levels = c("B2",  "C9",  "F1",  "A7" , "A8" , "B8" , "D9" , "H3" , "B4" , "C6" , "B11", "A11", "E5", 
                                 "G3" , "G6" , "A9" , "D8" , "A1" , "F9" , "F5" , "A3" , "F7" , "G9",  "A5" , "D6" , "A6", 
                                 "C10", "C4" , "F10", "C1", "D10", "F8",  "H8",  "D2",  "E11", "B3",  "E3",  "E10" ,"B10",
                                 "C2",  "D4",  "D1",  "F4",  "C5",  "H10","B6" , "E7",  "B1",  "E4",  "B5",  "E8",  "H7", 
                                 "D5",  "H6",  "B7",  "G11", "F11", "F6",  "D11" ,"D7" , "F2" , "H5" , "H9",  "G2" , "D3", 
                                 "E6" , "G4" , "F3",  "H11", "C11" ,"C7" , "C3",  "C8" , "A2" , "G7", "G8" , "H1",  "G5", 
                                 "A4",  "H4",  "E1",  "H2",  "E9" , "G1",  "A10", "B9" , "G10"))
df$V9 <- factor(df$V9, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                                  "15", "17", "18", "19","20", "21", "22", "23", "24", "25", "26",
                                  "27", "28","30", "32", "33", "Z" ))

df$equal_bin <- as.numeric(df$equal_bin)

ggplot(df, aes(x =ID, y =V5, fill= equal_bin)) + 
  geom_tile() + 
  scale_y_reverse() +
  facet_grid(. ~V9) +
  scale_fill_gradientn(colours = c('#801638','#027878','#FDB632','#F37338','#C22326'),trans = "reverse",
                       breaks=c(10,1),
                       labels=c("Low","High"), limits=c(10,1))+
  labs(x="", y ="Bin", fill = 'Score') +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())

###Draw density of rad sites###

#scale_fill_viridis() +
#scale_fill_brewer(palette = "RdYlBu") +  

ggplot(df, aes(x =ID, y =V5, fill= V6)) + 
  geom_tile() + 
  scale_y_reverse() +
  facet_grid(. ~V9) +
  scale_fill_gradientn(colours = c('#801638','#027878','#FDB632','#F37338','#C22326'))+
  labs(x="", y ="Bin", fill = '# Sites') +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank())


####Draw Chromosome 5 variation
chr5 <- subset(df, V9=="5")
ggplot(chr5, aes(x =ID, y =V5, fill= equal_bin)) + 
  geom_tile() + 
  scale_y_reverse() +
  scale_fill_gradientn(colours = c('#801638','#027878','#FDB632','#F37338','#C22326'),trans = "reverse",
                       breaks=c(10,1),
                       labels=c("Low","High"), limits=c(10,1))+
  labs(x="", y ="Bin", fill = 'Variation') +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(vjust =5),
        panel.background = element_rect(fill ="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)  


####Draw Chromosome 5 density
chr5 <- subset(df, df$V9 =="5")
ggplot(chr5, aes(x =ID, y =V5, fill= V6)) + 
  geom_tile() + 
  scale_y_reverse() +
  scale_fill_gradientn(colours = c('#801638','#027878','#FDB632','#F37338','#C22326'))+
  labs(x="", y ="Bin", fill = '# Sites') +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(vjust =5),
        panel.background = element_rect(fill ="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),)  








