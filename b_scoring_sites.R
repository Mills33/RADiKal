##Calculate the number of RAD sites used by RADiKal and how many in each scoring category/distribution

library(ggplot2)
library(tidyverse)
library(plotrix)

setwd("/Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_chr_results")

data <- read.csv("All_scores/All_RAD_scores.csv", header = F, sep = " ") # all sites 35332
falserad <- subset(data, V4 < 0) # 393
all_data<- data[!(data$V4 <= 0),] # only putative RAD (scored more than 0) 34940 

more175 <- subset(all_data, V4 > 175) # 11061
eighty8to175 <-  all_data %>% filter(V4 >88 &  V4 <= 175) #15464
less88 <- subset(all_data, V4 < 88) #8415

radornot <- c("n","r","r","r")
scores  <- c("<0","0-88", "88-175", ">175")
number_sites <- c(393,8415,15464,11061)
df <- as.data.frame(bind_cols(`# Sites`= number_sites, `Total score` = scores, col=radornot))

df$Scores <-  factor(df$Scores, levels = c("<0","0-88", "88-175", ">175"))

ggplot(df, aes(y =`# Sites` , x = `Total score`)) +
  geom_col(aes(fill = col),width =0.50) +
  scale_x_discrete(expand = c(0, 0)) +
  theme(legend.position ="none",
        panel.background = element_rect(fill ="white"),
        axis.line = element_line(colour = "black", size =0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("grey", "black")) + 
  scale_y_continuous(expand = c(0, 0)) 
