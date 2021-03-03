library(mltools)
library(proxy)
library(textshape)
library(tidyr)
setwd("/Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_scores_info/")
df <- read.csv("ALL_scores_info.csv", header = F, sep=" ")
df<- df[!(df$V4 <= 0),]
df[, "binme"] <- bin_data(df$V4, bins=10, binType = "quantile")

bins <- c(0.000214359356167648,35.4873191237245, 65.6608760600881, 87.4509361085066,
          102.441368790823, 115.551631313013, 137.093693977965, 158.511351118823,
  173.198888010529, 186.346588587725, 232.537542948543)
