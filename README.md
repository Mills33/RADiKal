# RADiKal 
Software package that processes and visualises RADSeq data to show the genetic diversity of a population (currently not a finished package)

[Draft]

Currently the scripts that make up RADiKal take raw unprocessed RAD reads and create painted chromosomes (see below) which gives and indication about the genetic variation (labelled score in image where a higher score indicates higher variation compared to a reference genome) present in the population. RADiKal produces several plots which all contribute to understanding the genetic diversity in a population.

![Painted_chromosome_variation](https://user-images.githubusercontent.com/23661194/109863479-9776a200-7c59-11eb-88aa-a022b106a5c1.png)


RADiKal uses kmers.


## Run Radikal (run_radikal.py, classifier.py)

```bash
source lmod-6.1
ml restore sdg_swig
mapfile -t INPUTS  < all_chromosome_list.txt
cd /ei/workarea/users/ryanc/2020RADiKal/output/PP_chr_${INPUTS[${SLURM_ARRAY_TASK_ID}]}
mkdir classifier_images scores
chrm=PP_chr_${INPUTS[${SLURM_ARRAY_TASK_ID}]}.fasta

python ../../scripts/run_radikal.py $chrm

mv *png classifier_images
mv *csv scores
```



## Calculate chromosome length and relative size

```bash
pwd : /hpc-home/ryanc/scratch/2020_Chapter3_popgen/data

for x in $(cat all_chromosome_list.txt ):
do
grep ">Chr_${x}_RagTag" -A 1  pp_ragtag_chr.fasta | tail | wc -c
done

grep '>ChrU' -A 1  pp_ragtag_chr.fasta | tail | wc -c
```



## Create All_RAD_scores.csv 

```bash
for x in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,30,32,33,Z}
do
awk -v chr=$x '{print $0 " " chr}' /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/scores/RAD_scores.csv > ./All_chr_results/${x}_RAD_scores.csv
cp /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/scores/FAKE_scores.csv ./All_chr_results/${x}_FAKE_scores.csv
mv /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/scores/RAD_scores.csv /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/
mv /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/scores/FAKE_scores.csv /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/
done

cd /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_chr_results
cat *RAD_scores.csv > All_RAD_scores.csv
cat *FAKE_scores.csv > All_FAKE_scores.csv
```



## Look at spread of scores (Spread_scores.r)



## Number of sites in each of the three peaks (scoring_sites.r)





## Calculate number of sites and the mean number sites compared to the relative to length of each chromosome (sites_density.r)



#############################################################################################

**The above looked at the scores produced from the scoring matrix when the classifier was being created and all individuals were catted together to represent a population file. The classifier is now used to score each individual in the population and those scores combined to assess the levels of variation within the population. The level of variation can be seen at-a-glance by plotting painted chromosomes.**

#############################################################################################

## Create all_scores

```
wd:/Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/
mkdir All_scores_info

for x in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,30,32,33,Z}
do
cat /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/scores/*csv > /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_scores_info/${x}_all_ind_scores.csv
done


cat /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_scores_info/*csv > /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/All_scores_info/ALL_scores_info.csv
```



## Create bins (Create_bins.r)



## Generate matrix which is used to paint chromosomes (Generate_matrix_for_painted_chromosome.r)

```
for x in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,25,26,27,28,30,32,33,Z}
do
cd /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/PP_chr_${x}/scores/
Rscript /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/scripts/Generate_matrix_for_painted_chromosome.r $x
mv *matrix*.csv /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/matrix_files
done

cd  /Users/ryanc/Documents/Pink_pigeon/RADIKal_v2/2020Aug/FINAL_RADIKAL/matrix_files
cat *matrix* > All_matrix.csv
```



## Paint chromosomes - Variation at-a-glance (painting_chromosomes.r)

- This also contains code for plotting the density of sites across the genome and Chromosome 5 individually



## Compare with whole genome reseqeuncing variation (Whole_genome_comparison.r)
