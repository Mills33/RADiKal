##Before running this script need to get bins by using infor from all sites combined df[, "wt2"] <- bin_data(df$V3, bins=10, binType = "quantile")###
## For a single chromosome this script collates the rad scores from each individual and combines them into a
## To run on all chromosomes can be set off in a bash for loop.
## After matrix created for each chromosome they can be catted together to form the final matrix which is used to paint chromosome
### The 10 bins used in this script have to be manually entered having been calculated in a previouse script

##input##
#10 Bins - start values
#chromosome length file
#individual rad scores for a chromosome or contig

### Load libraries ###
library(tidyverse)
library(proxy)
library(textshape)

###Prepare to generate matrix###
chr = commandArgs(trailingOnly=TRUE) # This takes the value (x)after Rscript generate_paint.py x where x is the name/number of a chromosome. In this way this script can be efficiently looped
fileNames <- Sys.glob("*.csv") # produces a list of all the files ending *csv - this script is run in the scores directory for each chromosome so lists the scores for each individuals
len_file <- read.csv("/Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/pp_chr_len.csv", header = F) # file where V1 is chromosme number, V2 length
length_chr <- len_file$V2[len_file$V1 == chr] # assigns specific chromosomes length


###Function generate matrix###
## To amplify signal took 50 observation windows (50 rad sites)
## infile = ind.scores.csv where  V1 = Start RAD position, V2 = score fwd , V3= score reverse , V4 = total score
paint_chrom_mat <-function(infile){
  data <- read.csv(infile, sep = ' ', header = F) #ind_scores.csv where ind differs depending on individual
  name <- sub('\\_scores.csv$', '', infile) # extract individual id
  if (nrow(data) >= 50){
    Chr_bin <- cut(1:length(data$V1), # divide the scores into groups of 50 - excpet for the last group which will be the remainder
                   breaks = seq(from = 1,
                   to = length(data$V1), by = 50),
                   include.lowest = T,
                   right = F,
                   labels = F) # divides the chromosome into 50 bins (with remainder)
    Chr_bin <- as.data.frame(Chr_bin) # dataframe containing the bin for each site. The bins will not be evenly distrubted among sites as the bins are based on the chr_len and it may be that if you divide 200mb chromosome into 50 so 4mb each bin bin 1 may have 50 sites but bin 2 may have 67 etc.
    Chr_bin <- Chr_bin %>% mutate_if(is.integer, ~replace(., is.na(.), max(Chr_bin, na.rm=T) +1)) # unlikely that the chromosome length is divisble by 50 the remainder sites are assigned NA therefore this replaces NAs with the (maximum bin value +1) so creates a bin at the end
    df <- bind_cols(data, Chr_bin) # combine data and bin so now V5 is Chr_bin
    df$Chr_bin <- as.factor(df$Chr_bin)} else {data$Chr_bin <- as.factor(1) # makes Chr_bin factor and if there were not 50 sites in chromosomes assigns them all to Chr_bin 1 as a factor
    df <- data} # change name of data to df
  df<- df[!(data$V4 <= 0),] # Remove putative non RAD any site with a score less than 0 is deleted - this will inevitably remove some real rad sites however in the future RADiKal outputs those rads with scores less than 0 so you can manually check sites
  df <- df %>%
    group_by(Chr_bin) %>%
    mutate(count = n()) # count how many sites in each of the 50 chromosome bins (the bins are equal in bp length but not in number of sites contained)
  bins <- c(0.000214359356167648,35.4873191237245, 65.6608760600881, 87.4509361085066,
           102.441368790823, 115.551631313013, 137.093693977965, 158.511351118823,
           173.198888010529, 186.346588587725, 232.537542948543) # These bins are pre-calculated using df[, "wt2"] <- bin_data(df$V3, bins=10, binType = "quantile") and manually entered although this will be altered. These bins represent 10 levels of scoring bin1 = highly het low scoring, bin 10 = high hom. high scoring. Bins generated from x file and each level should contain equal number of data points
  tags <- c(0,1, 2, 3, 4,5,6,7,8,9) # tags needed for binning represent 10 levels in $bins
  equal_bin <- cut(df$V4,
                   breaks= bins,
                   include.lowest=TRUE,
                   right=FALSE,
                   labels = tags) # assigns the scores to one of 10 bins 
  equal_bin <- as.data.frame(equal_bin) # create data frame
  df <- bind_cols(df, equal_bin) # add equal_bin to dataframe df
  df <- df %>%
    group_by(Chr_bin) %>%
    mutate(mean = sum(as.numeric(as.character(equal_bin)))/count) # creates column called mean and for each chr_bin takes the mean of 'equal bin' which is which of the 10 levels the score landed in
  df <-   distinct(df, Chr_bin, .keep_all = TRUE) # Have redundeant information now so reduces data to those with distinct values CHR_bin so length depends on how many bins your chromosome was divided into
  df$ID <- name # add id column
  df$chromosome <- chr # add column for chromosome
  df<- as.data.frame(df) # make sure whole thing is a dataframe
}

### Dataframe produced
#V1 = start, V2 = fwd score, V3 = reverse score, V4 = Total score, V5 = Chr_bin, V6 = number sites in each Chr_bin, V7 = equal bin (which of 10 levels), V8 = mean equal bin per chr_bin, V9 = ID, V10 = chromosome
#NB because only kept distinct values 1 per Chr_bin V1 = the position of the first start in that bin, V2,V3,V4 again not very useufl as only represent the first site in that bin, instead will be using the mean of 10 levels ot proceed
###
combined_df <- do.call(rbind, lapply(fileNames, function(fn) paint_chrom_mat(fn))) # takes each file in filename list and runs function and combines output matrix to create final data which is all individuals for a singel chromosome
write.table(combined_df,file=paste(chr,"_50_obvs_minus_10level_matrix.csv", sep=""), sep =",", col.names =F, row.names=FALSE) # output to csv
