import sys
import classifier
cl = classifier.Classifier() # initiates classifier creates empty SM
cl.create_classifier('GGACGTCC',sys.argv[1],'/ei/workarea/users/ryanc/2020RADiKal/data/PP_all_R1.fastq') # R1 from all individuals concatenated,  creates the classifier
file_ids = open("/ei/workarea/users/ryanc/2020RADiKal/scripts/file_ids.txt") # list of ids for every read file
for ind in file_ids:
    cl.score_individuals("/hpc-home/ryanc/camilla/PinkPigeon_data/RAD/140715_SN790_0359_BH9WVBADXX/Demultiplexed_lane1/RAW_data/radplex_out_"+str(ind).strip()+"_R1.fastq", ind)
    print("individual", ind, "complete. There are", len(cl.putative_sites), "remaining")
file_ids.close()
