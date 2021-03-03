import re
import weakref
import random
import sys
sys.path.append("/hpc-home/ryanc/sdg/build/pysdg/pysdg")
import pysdg
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
class Classifier:
    """
    The Classifier class is responsible for generating a forward strand and reverse strand relative frequency scoring matrix. Run using the function create_classifier().
    1. Formats cut site provided by user to ensure both possible motifs included in regex
    2. Loads reference genome and creates reverse complement using SDG
    3. Loads read ones from all individuals (prior to this they should have been catted into single file)
    4. Finds putative windows - this includes putative sites, ie sites with the cut motif, and the surrounding sequence. Each putative window becomes an instance of PutativeSite class which gets the site sequence string. projects kmer coverage and then creates a nomalised scoring scoring matrix
        each instance contains this information for both the forward strand and the reverse complement. For the reverse complement the coordinates of the sites have been altered such that the sites occur at the same coordinates.
    5. Two checks are performed on putative windows to ensure that they are RAD - these are if window all 0's discard and if window only has single value discard if one of these is true for either the forward or the reverse then the putative window is discarded.
    6. Add all the forward matrices an reverse matrices and create two relative site frequency matrices.
    """
    def __init__(self): # need cutter, bp_before and after depend on read length and buffer, K How to have defaults?
        # Settings / parameters
        self.cutter = None # Default.
        self.bp_before = 10 # buffer
        self.bp_after = 160 # read length + buffer
        self.coverage_k = 27 # K paramter, length of Kmers.
        self.cut_site_length = 8
        # Data / state
        self.fw_ws = pysdg.WorkSpace() # The SDG workspace containing the forward genome.
        self.rv_ws = pysdg.WorkSpace() # The SDG workspace containing the reverse of the genome.
        self.num_bins = 21
        self.putative_sites = [] # Will store locations of each cut-site +/- 250bp in the reference.
        self.ksl = ((self.bp_before + self.bp_after + self.cut_site_length) - self.coverage_k) + 1 #can i pass parameters from other class
        self.delete_site = []
        self.bins = []
        self.rel_freq_mat = []
        self.fake_rel_freq_mat = []
        self.fake_sites  = []
        self.fw_total_sm = []
        self.rv_total_sm = []
        self.fw_rel_freq_mat = []
        self.rv_rel_freq_mat = []
        self.global_scores = []
        self.fake_scores = []
    def set_cutter_motif(self, motif):
        """
        Users prove RADiKal with one version of the palindromic motif and this creates the complement and puts them in a regex fw_string
        required for the rest of the program.
        """
        cut_site_motif = motif +'|'+ motif[::-1]
        self.cut_site_length = len(motif)
        self.cutter = cut_site_motif
    # Load the reference genome in a forward and reverse orientation, in a
    # fw_ SDG workspace, and a rv_ SDG workspace.
    def load_reference_genome(self, file):
        """
        Uses SDG to load a genome into two workspaces, which were initialised as empty workspaces, and then use SDG function to get the reverse complement of the sequence in each node.
        """
        self.fw_ws.sdg.load_from_fasta(file)
        self.rv_ws.sdg.load_from_fasta(file)
        for n in range(self.rv_ws.sdg.nodes.size()):
            self.rv_ws.sdg.nodes[n].make_rc()
    def prep_rad_read_ones(self, fastq_file, main='main'):
        """
        Loads read ones into workspace and creates kmer counter and adds kmers counts for both the forward and reverse complemented workspaces
        """
        fw_rad_counts = self.fw_ws.add_kmer_counter(str(main)+"fw_pooled_rads_noncanon", self.coverage_k, pysdg.NonCanonical)
        fw_rad_counts.add_count(str(main)+"fw_coverage", [fastq_file])
        rv_rad_counts = self.rv_ws.add_kmer_counter(str(main)+"rv_pooled_rads_noncanon", self.coverage_k, pysdg.NonCanonical)
        rv_rad_counts.add_count(str(main)+"rv_coverage", [fastq_file])
    def save_workspace(self):
        self.fw_ws.dump_to_disk("RADiKal_fwd")
        self.rv_ws.dump_to_disk("RADiKal_rv")
    def create_empty_scoring_matrix(self):
        """
        Creates an empty scoring matrix used in several functions. The dimensions of the array will depend on the read length, buffer size and number of bins chosen
        """
        self.empty_sm = np.zeros([self.num_bins, self.ksl])
    def find_putative_windows(self):
        """
        Finds every instance in the forward genome that containes the cutter motif and add bases before (buffer) and after ( read length and buffer) to give  whole window size that allows the detection of the signal. Then for each putative site create an instance of
        PutatveSite which will get the string, project kmers, run checks and create normalised scoring matrix for each putative site.
        """
        del self.putative_sites[:]
        for nodeid in range(1, len(self.fw_ws.sdg.nodes)):
            for match in re.finditer(self.cutter, self.fw_ws.sdg.nodes[nodeid].sequence):
                coord = match.span()
                first = coord[0] #- self.bp_before
                last = coord[1] #+ self.bp_after
                kmer_end = last - 28 #adjust to kmer_counter
                self.putative_sites.append(PutativeSite(self, nodeid, first, last, kmer_end))
    def find_fake_windows(self, number):
        """
        Finds every instance in the forward genome that containes the cutter motif and add bases before (buffer) and after (read length and buffer) to give  whole window size that allows the detection of the signal. Then for each putative site create an instance of
        PutatveSite which will get the string, project kmers, run checks and create normalised scoring matrix for each putative site.
        """
        del self.fake_sites[:]
        num_per_node = number//(len(self.fw_ws.sdg.nodes) -1)
        for nodeid in range(1, len(self.fw_ws.sdg.nodes)):
            max = len(self.fw_ws.sdg.nodes[nodeid].sequence) - self.bp_after
            coord = random.sample(range(self.bp_before, max), num_per_node) # list of sites
            for match in coord:
                first = match - self.bp_before
                last = match + 8 + self.bp_after
                kmer_end = last - 29 #adjust to kmer_counter
                self.fake_sites.append(FakeSite(self, nodeid, first, last, kmer_end))
    def check_sites_1(self):
        """
        If true will mean the site is definitely not a RAD site. Condition 1 states that a site cannot be made of only 0s in either the forward or the reverse direction. This does not mean that all the remaining sites are RAD but this is a quick and obvious checks to run.
        """
        for site in reversed(range(0,len(self.putative_sites))): # change this
            condition_1 = sum(self.putative_sites[site].fw_counts) == 0 or sum(self.putative_sites[site].rv_counts) == 0
            if condition_1 == True:
                self.putative_sites.remove(self.putative_sites[site])
    def check_sites_2(self):
        """
        If true will mean the site is definitely not a RAD site. Condition 2 states that there cannot be a 1 in either the forward of reverse normalised kmer counts this would mean there would only be a single kmer in the window. This does not mean that all the remaining sites are RAD but this is a quick and obvious checks to run.
        """
        for site in reversed(range(0,len(self.putative_sites))):
            condition_2 = 1 in self.putative_sites[site].fw_normalised_count or 1 in self.putative_sites[site].rv_normalised_count
            if condition_2 == True:
                self.putative_sites.remove(self.putative_sites[site])
    def create_bins_for_scoring_matrix(self):
        """
        Creates bin values based on number of bins decided to have (default 20).
        """
        for x in range(0, 20):
            bin_value = (0.01/20)*x
            self.bins.append(bin_value)
        self.bins.append(0.01)
        self.bins.append(1) # can you initialise this whole function?
    def fw_add_scoring_matrix(self):
        """
        Adds every normalised forward scoring matrix from each instance of PutativeSite that passed the checks.
        """
        self.fw_total_sm = np.zeros([self.num_bins, self.ksl])
        for site in range(0, len(self.putative_sites)-1):
            self.fw_total_sm = self.putative_sites[site].fw_scoring_matrix + self.fw_total_sm
        return  self.fw_total_sm
    def rv_add_scoring_matrix(self):
        """
        Adds every normalised reverse scoring matrix from each instance of PutativeSite that passed the checks.
        """
        self.rv_total_sm = np.zeros([self.num_bins, self.ksl])
        for site in range(0, len(self.putative_sites)-1):
            self.rv_total_sm = self.putative_sites[site].rv_scoring_matrix + self.rv_total_sm
        return  self.rv_total_sm
    def fw_fake_add_scoring_matrix(self):
        """
        Adds every normalised forward scoring matrix from each instance of FakeSite that passed the checks.
        """
        self.fw_fake_total_sm = np.zeros([self.num_bins, self.ksl])
        for site in range(0, len(self.fake_sites)-1):
            self.fw_fake_total_sm = self.fake_sites[site].fw_scoring_matrix + self.fw_fake_total_sm
        return  self.fw_fake_total_sm
    def rv_fake_add_scoring_matrix(self):
        """
        Adds every normalised forward scoring matrix from each instance of FakeSite that passed the checks.
        """
        self.rv_fake_total_sm = np.zeros([self.num_bins, self.ksl])
        for site in range(0, len(self.fake_sites)-1):
            self.rv_fake_total_sm = self.fake_sites[site].rv_scoring_matrix + self.rv_fake_total_sm
        return  self.rv_fake_total_sm
    def fw_relative_freq(self):
        """
        Takes added scoring matrix from forward strand and creates a relative frequency matrix which is used as the forward classifier scoring matrix
        """
        fw_total_sm = self.fw_add_scoring_matrix()
        self.fw_rel_freq_mat = []
        for x in range(0, self.ksl):
            prop = []
            max = np.amax(fw_total_sm[:,x])
            for value in np.nditer(fw_total_sm[:,x]):
                if value < 0 :
                    prop.append(-1)
                else:
                    rel_freq = float(value)/max
                    prop.append(rel_freq) # convert the frequencies in each bin into raltive frequencies
            self.fw_rel_freq_mat.append(prop)
        self.fw_rel_freq_mat = np.transpose(np.array(self.fw_rel_freq_mat))
        #self.rel_freq_mat = np.column_stack((self.bins, transposed_rel_matrix))
    def rv_relative_freq(self):
        """
        Takes added scoring matrix from reverse complemented strand and creates a relative frequency matrix which is used as the reverse classifier scoring matrix
        """
        rv_total_sm = self.rv_add_scoring_matrix()
        self.rv_rel_freq_mat = []
        for x in range(0,self.ksl):
            prop = []
            for value in np.nditer(rv_total_sm[:,x]):
                if value < 0 :
                    prop.append(-1)
                else:
                    rel_freq = float(value)/np.amax(rv_total_sm[:,x])
                    prop.append(rel_freq) # convert the frequencies in each bin into raltive frequencies
            self.rv_rel_freq_mat.append(prop)
        self.rv_rel_freq_mat = np.transpose(np.array(self.rv_rel_freq_mat))
    def fw_fake_relative_freq(self):
        """
        Takes added scoring matrix from forward strand and creates a relative frequency matrix which is used as the forward classifier scoring matrix
        """
        fw_total_sm = self.fw_fake_add_scoring_matrix()
        self.fw_fake_rel_freq_mat = []
        for x in range(0, self.ksl):
            prop = []
            max = np.amax(fw_total_sm[:,x])
            for value in np.nditer(fw_total_sm[:,x]):
                if value < 0 :
                    prop.append(-1)
                else:
                    rel_freq = float(value)/max
                    prop.append(rel_freq) # convert the frequencies in each bin into raltive frequencies
            self.fw_fake_rel_freq_mat.append(prop)
        self.fw_fake_rel_freq_mat = np.transpose(np.array(self.fw_fake_rel_freq_mat))
        #self.rel_freq_mat = np.column_stack((self.bins, transposed_rel_matrix))
    def rv_fake_relative_freq(self):
        """
        Takes added scoring matrix from reverse complemented strand and creates a relative frequency matrix which is used as the reverse classifier scoring matrix
        """
        rv_total_sm = self.rv_fake_add_scoring_matrix()
        self.rv_fake_rel_freq_mat = []
        for x in range(0,self.ksl):
            prop = []
            for value in np.nditer(rv_total_sm[:,x]):
                if value < 0 :
                    prop.append(-1)
                else:
                    rel_freq = float(value)/np.amax(rv_total_sm[:,x])
                    prop.append(rel_freq) # convert the frequencies in each bin into raltive frequencies
            self.rv_fake_rel_freq_mat.append(prop)
        self.rv_fake_rel_freq_mat = np.transpose(np.array(self.rv_fake_rel_freq_mat))
    def add_rel_freq(self):
        """
        Takes the putative rad relative frequency matrix and the putative non-rad (fake) realtix frequency matrix and adds them together
        """
        self.fw_rel_freq_mat = self.fw_rel_freq_mat  - self.fw_fake_rel_freq_mat
        self.rv_rel_freq_mat = self.rv_rel_freq_mat - self.rv_fake_rel_freq_mat
    def scoring(self):
        """
        Should take an array of sites each row is a site, length is kmer_seq_len + 2 where the two columns have ID ingo
        """
        scores = []
        for x in range(0,len(self.putative_sites)):
            fwd_score = self.putative_sites[x].fw_score_site()
            rv_score = self.putative_sites[x].rv_score_site()
            location = self.putative_sites[x].first
            total = fwd_score + rv_score
            scores.append([location, fwd_score, rv_score, total])
        return(scores)
    def fake_scoring(self):
        """
        Should take an array of sites each row is a site, length is kmer_seq_len + 2 where the two columns have ID ingo
        """
        scores = []
        for x in range(0,len(self.fake_sites)):
            fwd_score = self.fake_sites[x].fw_score_site()
            rv_score = self.fake_sites[x].rv_score_site()
            location = self.fake_sites[x].first
            total = fwd_score + rv_score
            scores.append([location, fwd_score, rv_score, total])
        return(scores)
    def fw_create_heat_maps(self):
        ylabels = []
        for x in self.bins:
          ylabels.append(format(x, '.4f'))
          ylabels.reverse()
        df = np.flipud(self.fw_rel_freq_mat) # currently hard coded for rel_freq_matrix but will alter this
        ax = sns.heatmap(df, cmap="RdYlBu_r", cbar = False)
        ax.set_yticklabels(ylabels)
        plt.savefig("fw_scoring_matrix_heatmap")
    def rv_create_heat_maps(self):
        ylabels = []
        for x in self.bins:
          ylabels.append(format(x, '.4f'))
          ylabels.reverse()
        df = np.flipud(self.rv_rel_freq_mat) # currently hard coded for rel_freq_matrix but will alter this
        ax = sns.heatmap(df, cmap="RdYlBu_r", cbar = False)
        ax.set_yticklabels(ylabels)
        plt.savefig("rv_scoring_matrix_heatmap")
    def save_scores(self):
        f = open("RAD_scores.csv", "w+")
        for x in self.global_scores:
            print(*x, file = f)
        f.close()
        f_fake = open("FAKE_scores.csv", "w+")
        for y in self.fake_scores:
            print(*y, file = f_fake)
        f_fake.close()
    def plot_rad_under_0(self):
        RAD_under0 = []
        for x in range(0, len(self.putative_sites)):
           score = self.putative_sites[x].fw_score + self.putative_sites[x].rv_score
           total = sum(score)
           if  total < 0:
              RAD_under0.append(total)
              dff = np.flipud(self.putative_sites[x].fw_scoring_matrix)
              axf = sns.heatmap(dff, cmap="RdYlBu_r", cbar = False)
              plt.savefig(str(self.putative_sites[x].first) + "_RAD_fw_scoring_matrix_heatmap")
              dfr = np.flipud(self.putative_sites[x].rv_scoring_matrix)
              axr = sns.heatmap(dff, cmap="RdYlBu_r", cbar = False)
              plt.savefig(str(self.putative_sites[x].first) + "_RAD_rv_scoring_matrix_heatmap")
        plt.close()
    def plot_fake_over_0(self):
        Fake_over0 = []
        for x in range(0, len(self.fake_sites)):
           score = self.fake_sites[x].fw_score + self.fake_sites[x].rv_score
           total = sum(score)
           if  total > 0:
               Fake_over0.append(total)
               dff = np.flipud(self.fake_sites[x].fw_scoring_matrix)
               axf = sns.heatmap(dff, cmap="RdYlBu_r", cbar = False)
               plt.savefig(str(self.fake_sites[x].first) + "_FAKE_fw_scoring_matrix_heatmap")
               dfr = np.flipud(self.fake_sites[x].rv_scoring_matrix)
               axr = sns.heatmap(dff, cmap="RdYlBu_r", cbar = False)
               plt.savefig(str(self.fake_sites[x].first) + "_FAKE_rv_scoring_matrix_heatmap")
        plt.close()
    def create_classifier(self, cutter, genome, reads):
        """
        Work horse creates relative frequency scoring matrices for forward and reverse strands by running methods from the classifier class
        """
        self.set_cutter_motif(cutter)
        self.load_reference_genome(genome)
        self.prep_rad_read_ones(reads)
        self.save_workspace()
        print("Loaded genome, reads and workspace saved!")
        self.find_putative_windows()
        print("found", len(self.putative_sites),"potential RAD sites")
        self.create_empty_scoring_matrix()
        self.create_bins_for_scoring_matrix()
        print("Finished prepping, let's do this")
        for x in range(0, len(self.putative_sites)):
            self.putative_sites[x].fw_string()
            self.putative_sites[x].rv_string()
            self.putative_sites[x].project_coverage_of_fw_window()
            self.putative_sites[x].project_coverage_of_rv_window()
        self.check_sites_1()
        print("check 1 complete there are ", len(self.putative_sites), "remaining")
        for x in range(0,len(self.putative_sites)):
          self.putative_sites[x].normalise_fw_kmer_counts()
          self.putative_sites[x].normalise_rv_kmer_counts()
        self.check_sites_2()
        print("check 2 complete there are ", len(self.putative_sites), "remaining")
        for x in range(0,len(self.putative_sites)):
          self.putative_sites[x].create_fw_scoring_matrix()
          self.putative_sites[x].create_rv_scoring_matrix()
        self.find_fake_windows(len(self.putative_sites))
        for x in range(0, len(self.fake_sites)):
          self.fake_sites[x].fw_string()
          self.fake_sites[x].rv_string()
          self.fake_sites[x].project_coverage_of_fw_window()
          self.fake_sites[x].project_coverage_of_rv_window()
          self.fake_sites[x].normalise_fw_kmer_counts()
          self.fake_sites[x].normalise_rv_kmer_counts()
          self.fake_sites[x].create_fw_scoring_matrix()
          self.fake_sites[x].create_rv_scoring_matrix()
        self.fw_relative_freq()
        self.rv_relative_freq()
        self.fw_fake_relative_freq()
        self.rv_fake_relative_freq()
        self.add_rel_freq()
        self.global_scores = self.scoring()
        self.fake_scores = self.fake_scoring()
        self.save_scores()
        self.plot_rad_under_0()
        self.plot_fake_over_0()
    def score_individuals(self, fastq, ind):
        self.prep_rad_read_ones(fastq, main = ind)
        for x in range(0,len(self.putative_sites)):
          self.putative_sites[x].project_coverage_of_fw_window(main = ind)
          self.putative_sites[x].project_coverage_of_rv_window(main = ind)
          self.putative_sites[x].normalise_fw_kmer_counts(All_ind = False)
          self.putative_sites[x].normalise_rv_kmer_counts(All_ind = False)
        c = self.scoring()
        f = open(str(ind).strip()+"_scores.csv", "w+")
        for x in c:
          print(*x, file = f)
        f.close()
class PutativeSite:
    """
    This class is deisgned to take a single putative site and get its seqeucne string, project kmer coverage and create a normalised scoring matrix for both the forward and reverse genome
    """
    def __init__(self, classifier, nodeid, first, last, kmer_end): #creates a reference to memory where 'thing' stored that created when init triggered - used as a reference to link PutativeSite and Classifier due to special thing of init internal to python
        self.classifier = weakref.ref(classifier)
        self.nodeid = nodeid
        self.first = first
        self.last = last
        self.kmer_end = kmer_end
        self.fw_counts = []
        self.rv_counts = []
        self.fw_normalised_count = []
        self.rv_normalised_count = []
        self.fw_scoring_matrix = np.zeros([self.classifier().num_bins, self.classifier().ksl])
        self.rv_scoring_matrix = np.zeros([self.classifier().num_bins, self.classifier().ksl])
        self.fw_score = []
        self.rv_score = []
        self.total_score = self.fw_score + self.rv_score
    def fw_string(self):
        """
        Uses the coordinates passed to the class by the Classifier to retrieve a string of sequence for the putative window on the forward genome
        """
        fw_first = self.first - self.classifier().bp_before
        fw_last = self.last + self.classifier().bp_after
        fw_nodeseq = self.classifier().fw_ws.sdg.nodes[self.nodeid].sequence
        return fw_nodeseq[fw_first:fw_last]
    def rv_string(self):
        """
        Uses the coordinates passed to the class by the Classifier to retrieve a string of sequence for the putative window on the reverse genome in order to put the window in the same space in realtion to the forward window the corrdinates are alterered.
        """
        rv_nodeseq = self.classifier().rv_ws.sdg.nodes[self.nodeid].sequence
        rv_coord_f = len(rv_nodeseq) - self.last
        rv_coord_l = len(rv_nodeseq) - self.first
        rv_first = rv_coord_f - self.classifier().bp_before
        rv_last = rv_coord_l + self.classifier().bp_after
        return rv_nodeseq[rv_first:rv_last] # need to swap as last now smaller than first. If start was 3 and end was 6 and total 10 - 10 -3 = 7 and 10 -6 4 therefore first is now larger than last
    def project_coverage_of_fw_window(self, main='main'):
        """
        Projects kmers created from reads to the putative site's forward sequence string
        """
        fw_seq = self.fw_string()
        self.fw_counts = self.classifier().fw_ws.get_kmer_counter(str(main)+"fw_pooled_rads_noncanon").project_count(str(main)+"fw_coverage", fw_seq)
    def project_coverage_of_rv_window(self, main='main'):
        """
        Projects kmers created from reads to the putative site's reverse sequence string
        """
        rv_seq = self.rv_string()
        self.rv_counts = self.classifier().rv_ws.get_kmer_counter(str(main)+"rv_pooled_rads_noncanon").project_count(str(main)+"rv_coverage", rv_seq)
    def normalise_fw_kmer_counts(self, All_ind = True):
        """
        Normalises the kmer counts for the forward direction putative site
        """
        total = sum(self.fw_counts)# used to be all rows but from column 2. Column 1 = ID, column 2 = posnumber - is kmer_counter just like .cvg?
        if not All_ind:
            self.check_total_fw(total)
        else:
            self.fw_normalised_count = [x / total for x in self.fw_counts]
    def normalise_rv_kmer_counts(self, All_ind = True):
        """
        Normalises the kmer counts for the reverse direction putative site
        """
        total = sum(self.rv_counts)# used to be all rows but from column 2. Column 1 = ID, column 2 = posnumber - is kmer_counter just like .cvg?
        if not All_ind:
            self.check_total_rv(total)
        else:
            self.rv_normalised_count = [x / total for x in self.rv_counts]
    def check_total_rv(self, total):
        if total == 0:
            self.rv_normalised_count = [0 for x in self.rv_counts]
        else:
            self.rv_normalised_count = [x / total for x in self.rv_counts]
    def check_total_fw(self, total):
        if total == 0:
            self.fw_normalised_count = [0 for x in self.fw_counts]
        else:
            self.fw_normalised_count = [x / total for x in self.fw_counts]
    def create_fw_scoring_matrix(self): # make argument window id and put functions inside so you can do straight away just singel site node id 2 etc.
        """
        Uses the normalised values and the different bins created in Classifier to create a scoring matrix where the bin where the value lands get +1 - this is done for the forward direction.
        """
        site = self.fw_normalised_count
        window_scoring_matrix = self.fw_scoring_matrix # these dimensions will not be hardcoded eventually
        counter = -1
        for position in site:
            counter += 1
            value = float(position)
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        window_scoring_matrix[i, counter] += 1 # row = bin but column = position) # 0 in table is row names so everything starts 1 up
        self.fw_scoring_matrix = window_scoring_matrix # This is the thing we need!
    def create_rv_scoring_matrix(self): # make argument window id and put functions inside so you can do straight away just singel site node id 2 etc.
        """
        Uses the normalised values and the different bins created in Classifier to create a scoring matrix where the bin where the value lands get +1 - this is done for the reverse direction.
        """
        site = self.rv_normalised_count
        window_scoring_matrix = self.rv_scoring_matrix # these dimensions will not be hardcoded eventually
        counter = -1
        for position in site:
            counter += 1
            value = float(position)
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        window_scoring_matrix[i, counter] += 1 # row = bin but column = position) # 0 in table is row names so everything starts 1 up
        self.rv_scoring_matrix = window_scoring_matrix # This  This is the thing we need!
    def fw_create_heat_maps(self):
        ylabels = []
        name = self.first
        for x in self.classifier().bins:
          ylabels.append(format(x, '.4f'))
          ylabels.reverse()
        df = np.flipud(self.fw_scoring_matrix) # currently hard coded for rel_freq_matrix but will alter this
        ax = sns.heatmap(df, cmap="RdYlBu_r", cbar = False)
        ax.set_yticklabels(ylabels)
        plt.savefig("fw_" + str(name) + "_heatmap")
    def rv_create_heat_maps(self):
        ylabels = []
        name = self.first
        for x in self.classifier().bins:
          ylabels.append(format(x, '.4f'))
          ylabels.reverse()
        df = np.flipud(self.rv_scoring_matrix) # currently hard coded for rel_freq_matrix but will alter this
        ax = sns.heatmap(df, cmap="RdYlBu_r", cbar = False)
        ax.set_yticklabels(ylabels)
        plt.savefig("rv_" + str(name) + "_heatmap")
    def fw_score_site(self):
        """
        Takes a single site and gives score for that site by comparing
        site to scoring matrix
        """
        score = 0
        for position in range(0, (len(self.fw_normalised_count)-1)): # range works lower bound to 1- upper bound (up to not including upper bound)
            value = float(self.fw_normalised_count[position])
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        score += self.classifier().fw_rel_freq_mat[i, position ] # row = bin but column = position) # 0 in table is row names so everything starts 1 up
                    else:
                        pass
        self.fw_score.append(score)
        return(score)
    def rv_score_site(self):
        """
        Takes a single site and gives score for that site by comparing
        site to scoring matrix
        """
        score = 0
        for position in range(0, (len(self.rv_normalised_count)-1)): # range works lower bound to 1- upper bound (up to not including upper bound)
            value = float(self.rv_normalised_count[position])
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        score += self.classifier().rv_rel_freq_mat[i, position ] # row = bin but column = position) # 0 in table is row names so everything starts 1 up
                    else:
                        pass
        self.rv_score.append(score)
        return(score)
class FakeSite:
    """
    This class is deisgned to take a single putative site and get its seqeucne string, project kmer coverage and create a normalised scoring matrix for both the forward and reverse genome
    """
    def __init__(self, classifier, nodeid, first, last, kmer_end): #creates a reference to memory where 'thing' stored that created when init triggered - used as a reference to link PutativeSite and Classifier due to special thing of init internal to python
        self.classifier = weakref.ref(classifier)
        self.nodeid = nodeid
        self.first = first
        self.last = last
        self.kmer_end = kmer_end
        self.fw_counts = []
        self.rv_counts = []
        self.fw_normalised_count = []
        self.rv_normalised_count = []
        self.fw_scoring_matrix = np.zeros([self.classifier().num_bins, self.classifier().ksl])
        self.rv_scoring_matrix = np.zeros([self.classifier().num_bins, self.classifier().ksl])
        self.fw_score = []
        self.rv_score = []
    def fw_string(self):
        """
        Uses the coordinates passed to the class by the Classifier to retrieve a string of sequence for the putative window on the forward genome
        """
        fw_nodeseq = self.classifier().fw_ws.sdg.nodes[self.nodeid].sequence
        return fw_nodeseq[self.first:self.last]
#Need to reverse the coordinates so that the same poistion on fwd and rev genome being targeted
    def rv_string(self):
        """
        Uses the coordinates passed to the class by the Classifier to retrieve a string of sequence for the putative window on the reverse genome in order to put the window in the same space in realtion to the forward window the corrdinates are alterered.
        """
        rv_nodeseq = self.classifier().rv_ws.sdg.nodes[self.nodeid].sequence
        rv_first = len(rv_nodeseq) - self.first
        rv_last = len(rv_nodeseq) - self.last
        return rv_nodeseq[rv_last:rv_first] # need to swap as last now smaller than first. If start was 3 and end was 6 and total 10 - 10 -3 = 7 and 10 -6 4 therefore first is now larger than last
###get and normalise kmer coverage
    def project_coverage_of_fw_window(self, main = 'main'):
        """
        Projects kmers created from reads to the putative site's forward sequence string
        """
        fw_seq = self.fw_string()
        self.fw_counts = self.classifier().fw_ws.get_kmer_counter(str(main)+"fw_pooled_rads_noncanon").project_count(str(main)+"fw_coverage", fw_seq)
    def project_coverage_of_rv_window(self, main = 'main'):
        """
        Projects kmers created from reads to the putative site's reverse sequence string
        """
        rv_seq = self.rv_string()
        self.rv_counts = self.classifier().rv_ws.get_kmer_counter(str(main)+"rv_pooled_rads_noncanon").project_count(str(main)+"rv_coverage", rv_seq)
    def normalise_fw_kmer_counts(self):
        """
        Normalises the kmer counts for the forward direction putative site
        """
        total = sum(self.fw_counts)# used to be all rows but from column 2. Column 1 = ID, column 2 = posnumber - is kmer_counter just like .cvg?
        if total == 0:
            self.fw_normalised_count = [0 for x in self.fw_counts]
        else:
            self.fw_normalised_count = [x / total for x in self.fw_counts]
    def normalise_rv_kmer_counts(self):
        """
        Normalises the kmer counts for the reverse direction putative site
        """
        total = sum(self.rv_counts)# used to be all rows but from column 2. Column 1 = ID, column 2 = posnumber - is kmer_counter just like .cvg?
        if total == 0:
            self.rv_normalised_count = [0 for x in self.rv_counts]
        else:
            self.rv_normalised_count = [x / total for x in self.rv_counts]
    def create_fw_scoring_matrix(self): # make argument window id and put functions inside so you can do straight away just singel site node id 2 etc.
        """
        Uses the normalised values and the different bins created in Classifier to create a scoring matrix where the bin where the value lands get +1 - this is done for the forward direction.
        """
        site = self.fw_normalised_count
        window_scoring_matrix = self.fw_scoring_matrix # these dimensions will not be hardcoded eventually
        counter = -1
        for position in site:
            counter += 1
            value = float(position)
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        window_scoring_matrix[i, counter] += 1 # row = bin but column = position) # 0 in table is row names so everything starts 1 up
        self.fw_scoring_matrix = window_scoring_matrix # This is the thing we need!
    def create_rv_scoring_matrix(self): # make argument window id and put functions inside so you can do straight away just singel site node id 2 etc.
        """
        Uses the normalised values and the different bins created in Classifier to create a scoring matrix where the bin where the value lands get +1 - this is done for the reverse direction.
        """
        site = self.rv_normalised_count
        window_scoring_matrix = self.rv_scoring_matrix # these dimensions will not be hardcoded eventually
        counter = -1
        for position in site:
            counter += 1
            value = float(position)
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        window_scoring_matrix[i, counter] += 1 # row = bin but column = position) # 0 in table is row names so everything starts 1 up
        self.rv_scoring_matrix = window_scoring_matrix # This  This is the thing we need!
    def fw_score_site(self):
        """
        Takes a single site and gives score for that site by comparing
        site to scoring matrix
        """
        score = 0
        for position in range(0, (len(self.fw_normalised_count)-1)): # range works lower bound to 1- upper bound (up to not including upper bound)
            value = float(self.fw_normalised_count[position])
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        score += self.classifier().fw_rel_freq_mat[i, position ] # row = bin but column = position) # 0 in table is row names so everything starts 1 up
                    else:
                        pass
        self.fw_score.append(score)
        return(score)
    def rv_score_site(self):
        """
        Takes a single site and gives score for that site by comparing
        site to scoring matrix
        """
        score = 0
        for position in range(0, (len(self.rv_normalised_count)-1)): # range works lower bound to 1- upper bound (up to not including upper bound)
            value = float(self.rv_normalised_count[position])
            for i in range(0,len(self.classifier().bins)):
                if float(self.classifier().bins[i]) < float(self.classifier().bins[-1]): # if the bin is less than the last bin
                    start = float(self.classifier().bins[i])
                    end = float(self.classifier().bins[i+1])
                    if start <= value < end:
                        score += self.classifier().rv_rel_freq_mat[i, position ] # row = bin but column = position) # 0 in table is row names so everything starts 1 up
                    else:
                        pass
        self.rv_score.append(score)
        return(score)
