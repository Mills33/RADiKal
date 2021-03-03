library(ggplot2)
setwd("/Users/ryanc/Documents/Pink_pigeon/Population_genetics/Heterozygosity/SeptFilteredDensity/bcftools_split_removehom/")

df <- read.table("All_snpden.csv", header = TRUE)

summary(df$VARIANTS.KB)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   2.831   3.700   3.625   4.548   8.398

df$BIN_START <- as.factor(df$BIN_START)
df <- na.omit(df)
df$CHROM <- factor(df$CHROM, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                                  "15", "17", "18", "19","20", "21", "22", "23", "24", "25", "26",
                                  "27", "28","30", "32", "33" ))

odd_list <- c("1","3","5","7","9","11","13","15", "18","20","22", "24","26","28","32")
size <- c("1","2","3","4","5","6","7", "8","9")

df$col <- 0
for (i in df$CHROM){
  if (i %in% odd_list){ df$col[df$CHROM == i] = 1 }}


samples <- c("S2B7805","S1272", "S3B7462", "S6W1688","S7W1687", "S4B7703","pp_novaseq")


####Population####
macro <- subset(df, CHROM %in% size)
micro <- subset(df, !CHROM %in% size)

chr5 <- subset(df, CHROM==5)

ggplot(macro, aes(x = as.numeric(as.character(BIN_START)), y = VARIANTS.KB)) +
  geom_col(aes(colour = factor(col)),position='identity', width=0.1) +
  geom_hline(yintercept=3.625,linetype="dashed") +
  theme(legend.position = "none",
        panel.background = element_rect(fill ="white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  labs(y="Variants/kb", x ="Chromosomes") +
  facet_grid(. ~as.factor(CHROM),scales = "free_x", space = "free") +
  scale_colour_manual(values = c("grey69", "grey69")) +
  scale_y_continuous(breaks=seq(0,20, 4), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(panel.spacing.x=unit(0, "lines")) 

ggplot(micro, aes(x = as.numeric(as.character(BIN_START)), y = VARIANTS.KB)) +
  geom_col(aes(colour = factor(col), fill=factor(col)),position='identity') +
  theme(legend.position = "none",
        panel.background = element_rect(fill ="white"),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  labs(y="Variants/kb", x ="Chromosomes") +
  facet_grid(. ~as.factor(CHROM),scales = "free_x", space = "free") +
  scale_colour_manual(values = c("steelblue", "lightblue")) +
  scale_fill_manual(values = c("steelblue", "lightblue")) +
  scale_y_continuous(breaks=seq(0,20, 4), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(panel.spacing.x=unit(0, "lines")) 

ggplot(chr5, aes(x = as.numeric(as.character(BIN_START)), y = VARIANTS.KB)) +
  geom_col(aes(colour = factor(col), fill=factor(col)),position='identity') +
  theme(legend.position = "none",
        panel.background = element_rect(fill ="white"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  labs(y="Variants/kb", x ="Position along chromosome") +
  facet_grid(. ~as.factor(CHROM),scales = "free_x", space = "free") +
  scale_colour_manual(values = c("steelblue", "lightblue")) +
  scale_fill_manual(values = c("steelblue", "lightblue")) +
  scale_y_continuous(breaks=seq(0,10, 2), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0), labels = comma) + 
  theme(panel.spacing.x=unit(0, "lines"),
        strip.text.x = element_blank()) 
