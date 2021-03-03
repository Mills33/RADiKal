library(ggplot2)
###global Spread scores ###
setwd("/Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_chr_results/try_RAD/")
rad <- read.table("All_RAD_scores.csv", header = F, sep = " ")
fake <- read.table("All_FAKE_scores.csv", header = F, sep = " ")
rad$type <- 'RAD'
rad$name <- 'rad'
fake$type <- 'Non-RAD'
fake$name <- 'Non-RAD'
both <- rbind(rad, fake)

both_scores <- ggplot(both, aes(x=V4, fill = type)) +
  geom_histogram(position="dodge", alpha=0.75, bins = 100) + 
  labs(x = "Scores", fill = 'Data type', y="Count") +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-250, 250, 25)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        legend.position = "none",
        panel.background = element_rect(fill ="white"),
        axis.line = element_line(colour = "black", size =0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_y_continuous(expand = c(0, 0)) 

fake_scores <- ggplot(fake, aes(x=V4)) +
  geom_histogram(position="dodge", alpha=0.75, bins = 100)+ 
  labs(x = "Scores", fill = 'Data', y = "Count") 
rad_scores <- ggplot(rad, aes(x=V4)) +
  geom_histogram(position="dodge", alpha=0.75, bins = 100) + 
  labs(x = "Scores", fill = 'Data')  

both_scores
fake_scores
rad_scores

ggsave('RAD_scores.png', rad_scores)
ggsave('FAKE_scores.png', fake_scores)
ggsave('radvsfake_scores.png', both_scores)
