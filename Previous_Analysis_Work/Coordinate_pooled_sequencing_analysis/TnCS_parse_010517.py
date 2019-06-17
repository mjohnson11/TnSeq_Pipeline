#A program to demultiplex reads and count tn insertions for the Tn-seq project
#Milo Johnson
#Modified from SLAMseq_parse 11/04/2016
#Modified 01/05/2017

import time
import gzip
import csv
import subprocess
import os
import numpy as np
import Levenshtein

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('primer_map_file', help='file containing primer demultiplex info')
parser.add_argument('R2_file_base', help='the file base before the multiplex info and before the R2.fastq.gz')
parser.add_argument('output_file_base', help='output file base directory')
parser.add_argument('run_name', help='run name for file naming')
parser.add_argument('-sample_id_key', type=int, help='index of the sample to work on (from the list of samples in the demultiplex info file).'
                                                     'Input 1000 to run all samples sequentially.')
parser.add_argument("-b", help="run the count combining step (should be called only after all counting is complete)", action="store_true")
parser.add_argument("-run_combine_when_done", default="no", help="provide the .sh file to run when all libraries have been counted")
parser.add_argument("-bases_to_use", type=int, default=50, help='number of bases to use to define the insertion, default=50')
parser.add_argument("-qual_cutoff", type=int, default=30, help='average quality score cutoff for the barcode region, default=30')
parser.add_argument("-editdist_thresh", type=int, default=4, help='edit distance threshold for clustering within UMI families, default=4')

args = parser.parse_args()

###POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
PRIMER_MAP_FILE = args.primer_map_file
R2_FILE_BASE = args.R2_file_base
OUTPUT_FILE_BASE = args.output_file_base
RUN_NAME = args.run_name

###OPTIONAL ARGS
SAMPLE_ID_KEY = args.sample_id_key
DO_COMBINE = args.b
RUN_COMBINE_WHEN_DONE = args.run_combine_when_done
BASES_TO_USE = args.bases_to_use
QUALITY_CUTOFF = args.qual_cutoff
EDITDIST_THRESH = args.editdist_thresh

###HARDCODED OPTIONS
#expected position of the first genomic dna base after the transposon edge
EXPECTED_PREVIOUS_8 = 'CGCCCACA' #just to double check it is working
OPPOSITE_ADAPTER_SEQ = 'CTGTCTCTTA'
TN_EDGE_RC = 'TGTGGGCGGA'

TITLE_OUT = ['Edge.Bases', 'Count', 'Example.Full.Read']
UMI_OUT_TITLE = ['UMI', 'Edge.Bases', 'Counts']


def get_demultiplex_info(primer_map_file):
    # THE DEMULTIPLEX FILE SHOULD BE FORMATTED LIKE:
    # title_row
    # library_name, f_index, r_index, combining_group
    demult_info = []
    with open(primer_map_file, 'r') as infile:
        reader = csv.reader(infile)
        rc = 0
        for row in reader:
            if rc > 0:
                demult_info.append(row)
            rc += 1
    return demult_info

class DNACounter:

    def __init__(self, library_name):
        self.library_name = library_name
        #counts for various outcomes for each read
        self.full_count = 0 #any read
        self.not_tn_edge_count = 0 #did not identify the read as an actual tn7 edge
        self.quality_count = 0 #quality score in bc region was too low
        self.exactly_right_edge = 0 #exactly correct tn edge - indicates a definitely good read
        self.tag_too_close = 0 #tagmentation was too close, not enough edge sequence
        self.umi_count = 0 #counts excluded for

        self.dna_count_dict = dict()
        self.UMI_dict = dict()

    def add_entry(self, r2_tmp, r1_tmp, start_of_edge):
        #input is an entry (4 lines from 2 paired-end files), time to start checking things...
        #start of edge and the read with the edge will be different for my barcoded and unbarcoded libraries
        #this makes sure I'm looking in the right place
        if start_of_edge == 110:
            f_read = r1_tmp
            r_read = r2_tmp
        else:
            f_read = r2_tmp
            r_read = r1_tmp

        if TN_EDGE_RC in r_read: #looking for the tn edge in the reverse read
            if r_read.index(TN_EDGE_RC) < BASES_TO_USE: #tagmentation was too close to get a complete edge sequence
                self.tag_too_close += 1
        else:
            the_read = f_read[1]
            #first check: does this look like the tn7 edge: (exact match of the 8 bases at the edge)
            if not EXPECTED_PREVIOUS_8 in the_read[start_of_edge-16:start_of_edge+8]:
                self.not_tn_edge_count += 1
            else:
                found_edge = the_read.index(EXPECTED_PREVIOUS_8) + 8
                #check for quality of the region we'll be using
                if np.mean([ord(c)-33 for c in f_read[3][found_edge:found_edge+BASES_TO_USE]]) < QUALITY_CUTOFF:
                    self.quality_count += 1
                else:
                    #we're all good, count it
                    dna_region = the_read[found_edge:found_edge+BASES_TO_USE]
                    if len(dna_region) < BASES_TO_USE:
                        #this check is in place because there are some rare shortened reads (<150 bp)
                        self.quality_count += 1
                    elif OPPOSITE_ADAPTER_SEQ in the_read[:found_edge+BASES_TO_USE]:
                        self.tag_too_close += 1 #alternated way to catch the tagmentation being too close
                    else:
                        #the umi or unique molecular index is the 12 bases that start the opposite read
                        #since these bases come from tagmentation, seeing the same exact bases likely means
                        #it is a pcr duplicate, so I will not count it
                        tmp_UMI = r_read[1][:12]
                        if tmp_UMI in self.UMI_dict.keys(): #seen UMI
                            tmp_dict = self.UMI_dict[tmp_UMI]
                            if dna_region in tmp_dict.keys():
                                tmp_dict[dna_region][2] += 1
                            else:
                                #Entries like [Edge.Bases, read example, counts]
                                tmp_dict[dna_region] = [dna_region, the_read[:len(the_read)-1], 1]
                        else: #new UMI
                            self.UMI_dict[tmp_UMI] = {dna_region: [dna_region, the_read[:len(the_read)-1], 1]}

    def UMI_dict_parse(self):
        # This function parses the UMI dicts, and has 2 goals:
        # 1) create a list of BC counts where each count is from a unique UMI (UMI duplicates get removed)
        # 2) output statistics on UMI family sizes

        #this has entries like: [UMI R1, UMI R2, bc div, bc env, count]
        self.UMI_family_size_list = []

        #['Edge.Bases', 'Count', 'Example.Full.Read']

        for UMI in self.UMI_dict.keys():
            tmp_dict = self.UMI_dict[UMI]
            tmp_keys = [i for i in tmp_dict.keys()]
            if len(tmp_keys) == 1:
                #only one dna region corresponding to this UMI, simply add it into the bc dict and record the UMI fam
                use_region = tmp_keys[0]
                use_rec = tmp_dict[use_region]
                self.umi_count += use_rec[2] - 1
                self.UMI_family_size_list.append([UMI, use_rec[0], use_rec[2]])
                if use_region in self.dna_count_dict.keys():
                    self.dna_count_dict[use_region][1] += 1 #add 1 to corrected counts
                else:
                    self.dna_count_dict[use_region] = [use_rec[0], 1, use_rec[1]]
            else:
                #multiple barcodes for this UMI, error correct so we only count each pcr template once
                regions_in_umi_family = []
                #sort so that errors are clustered to more common bcs
                sorted_keys = sorted(tmp_keys, key=lambda x: -1*tmp_dict[x][2])
                for region in sorted_keys:
                    region_rec = tmp_dict[region]
                    matched = False
                    for region_target in regions_in_umi_family:
                        #check to see if this region is really an error off of another one within this UMI family
                        if Levenshtein.distance(region_rec[0], region_target[0]) <= EDITDIST_THRESH:
                            matched = True
                            #add counts to UMI family
                            region_target[2] += region_rec[2]
                            break
                    if not matched:
                        #didn't match, new region centroid
                        regions_in_umi_family.append(region_rec)
                #add counts for these
                for region_rec in regions_in_umi_family:
                    self.UMI_family_size_list.append([UMI, region_rec[0], region_rec[2]])
                    self.umi_count += region_rec[2] - 1
                    if region_rec[0] in self.dna_count_dict.keys():
                        self.dna_count_dict[region_rec[0]][1] += 1         #add 1 to corrected counts
                    else:
                        self.dna_count_dict[region_rec[0]] = [region_rec[0], 1, region_rec[1]]     #new barcode, add full rec + initialize corrected counts

    def read_file(self, f_in, start_of_edge):
        f_in_R1 = f_in.replace('R2.fastq', 'R1.fastq')
        with gzip.open(f_in_R1, 'rt') as r1in:
            with gzip.open(f_in, 'rt') as r2in:
                rc = 0
                r2_lines_tmp = []
                r1_lines_tmp = []
                for r2line in r2in:
                    r1_lines_tmp.append(r1in.readline())
                    r2_lines_tmp.append(r2line)
                    rc += 1
                    if (rc % 4) == 0:
                        self.add_entry(r2_lines_tmp, r1_lines_tmp, start_of_edge)
                        r2_lines_tmp = []
                        r1_lines_tmp = []
        print('Looked at', rc/4, 'reads.')
        self.full_count = rc/4

    def output_counts(self, fout, fout_UMI_fam):
        #outputting the file
        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(TITLE_OUT)
            for dna_region in self.dna_count_dict.keys():
                writer.writerow(self.dna_count_dict[dna_region])

        #output UMI family sizes
        with open(fout_UMI_fam, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(UMI_OUT_TITLE)
            for UMI_rec in self.UMI_family_size_list:
                writer.writerow(UMI_rec)

        with open(masterout, 'a') as outfile:
            writer = csv.writer(outfile)
            writer.writerow([self.library_name, self.full_count, self.not_tn_edge_count, self.quality_count,
                             self.tag_too_close, self.umi_count,
                             self.full_count - self.quality_count - self.not_tn_edge_count - self.tag_too_close - self.umi_count])

class UMIFamSizeCombiner:

    def __init__(self):
        #the way this will be organized is bc_dict[bc][column name] = value
        self.biggest_fam = 1
        self.umi_size_dict = dict()
        self.libs_in_order = []

    def add_sample(self, fin, sample_name):
        self.libs_in_order.append(sample_name)
        tmp_dict = self.umi_size_dict[sample_name] = dict()
        print('Looking at file:', fin)
        with open(fin, 'r') as infile:
            reader = csv.reader(infile)
            rc = 0
            for row in reader:
                #skip title row
                if rc > 0:
                    fam_size = int(row[2])
                    if fam_size > self.biggest_fam:
                        self.biggest_fam = fam_size
                    if fam_size in tmp_dict.keys():
                        tmp_dict[fam_size] += 1
                    else:
                        tmp_dict[fam_size] = 1
                rc += 1

    def write_data(self, fout):

        top_val = self.biggest_fam

        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            fam_sizes = [i+1 for i in range(top_val)]
            writer.writerow(['Library'] + fam_sizes)
            for lib in self.libs_in_order:
                tmp_row = [lib]
                for fs in fam_sizes:
                    if fs in self.umi_size_dict[lib].keys():
                        tmp_row.append(self.umi_size_dict[lib][fs])
                    else:
                        tmp_row.append(0)
                writer.writerow(tmp_row)

class FileCombiner:

    def __init__(self):
        #the way this will be organized is bc_dict[bc][column name] = value
        self.column_names = ['Edge.Bases', 'Total.Count', 'Example.Full.Read']
        self.count_names = []
        self.data = dict()

    def add_sample(self, fin, sample_name):
        self.count_names.append(sample_name)
        print('Looking at file:', fin)
        with open(fin, 'r') as infile:
            reader = csv.reader(infile)
            rc = 0
            for row in reader:
                #skip title row
                if rc > 0:
                    dna_region = row[0]
                    counts = int(row[1])
                    if not (dna_region in self.data.keys()):
                        self.data[dna_region] = {'Edge.Bases': dna_region, 'Total.Count': counts,
                                                 'Example.Full.Read': row[2], sample_name: counts}
                    else:
                        self.data[dna_region]['Total.Count'] += counts
                        self.data[dna_region][sample_name] = counts
                rc += 1

    def write_data(self, fout):
        #sort by total counts
        use_data = self.data
        sorted_keys = sorted(use_data.keys(), key=lambda x: -1*use_data[x]['Total.Count'])

        with open(fout, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(self.column_names + self.count_names)
            for dna_region in sorted_keys:
                entry = use_data[dna_region]
                writer_row = []
                for colname in (self.column_names + self.count_names):
                    writer_row.append(entry.setdefault(colname, 0))
                writer.writerow(writer_row)


def run_one_sample():

    #for timing data
    orig_time = time.time()

    #sample id key input is 1 indexed
    dc = 1
    for dem in demult_info:
        #sample id of 1000 means run them all (no jabbaraying)
        if dc == SAMPLE_ID_KEY or SAMPLE_ID_KEY == 1000:
            print('Counting reads for:', dem[0])
            print(R2_FILE_BASE + dem[0].replace('_', '-') + '_S' + dem[3] + '.R2.fastq.gz')
            already_demult_file = R2_FILE_BASE + dem[0].replace('_', '-') + '_S' + dem[3] + '.R2.fastq.gz'
            if os.path.exists(already_demult_file):
                print('Demultiplexed file found, reading from:', already_demult_file)
                counter = DNACounter(dem[0])
                counter.read_file(already_demult_file, int(dem[4]))
                counter.UMI_dict_parse()
                counter.output_counts(count_file_dir + '/' + dem[0] + '_counts.csv', count_file_dir + '/' + dem[0] + '_UMI_families.csv')
            else:
                print('Demultiplexed file not found! exiting...')
                exit()
        dc += 1

    print('Done!  Time elapsed:', time.time()-orig_time)

def get_all_libraries_done():
    libraries_done = dict()
    with open(masterout, 'r') as infile:
        reader = csv.reader(infile)
        for row in reader:
            if row[0] in libraries_done.keys():
                print('\n\n\nWARNING: IT LOOKS LIKE', row[0], 'HAS BEEN COUNTED TWICE')
            libraries_done[row[0]] = row
    return libraries_done

def finish_and_get_statistics(combine_name):
    #If all libraries are done, format the statistics file
    libraries_done = get_all_libraries_done()
    if len(libraries_done) == len(demult_info):
        print('All libraries are completed:\nReformatting library stats file...')
        with open(masterout[:masterout.index('.csv')] + '_formatted.csv', 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Library', 'Total.Reads', 'Did.Not.Find.Tn.Edge', 'Quality.Failed',
                             'Opposite.Adapter.Found', 'UMI.counts', 'Usable.Reads'])
            for dem in demult_info:
                writer.writerow(libraries_done[dem[0]])
        print('Combining count files...')
        combinerator = FileCombiner()
        UMIcomb = UMIFamSizeCombiner()
        for dem in demult_info:
            combinerator.add_sample(count_file_dir + '/' + dem[0] + '_counts.csv', dem[0])
            UMIcomb.add_sample(count_file_dir + '/' + dem[0] + '_UMI_families.csv', dem[0])
        combinerator.write_data(combine_file_dir + '/' + combine_name + '_combined_counts.csv')
        UMIcomb.write_data(combine_file_dir + '/' + combine_name + '_UMI_familysize_distribs.csv')
    else:
        print('Not combining statistics because some libraries from the demultiplex file are not done.')


#
#MAKING DIRECTORIES
#


if not os.path.isdir(OUTPUT_FILE_BASE):
    print('Making main output directory:', OUTPUT_FILE_BASE)
    subprocess.call(['mkdir', OUTPUT_FILE_BASE])
else:
    print('Making main output directory was already made:', OUTPUT_FILE_BASE)


masterout = OUTPUT_FILE_BASE + RUN_NAME + '_library_statistics.csv'
count_file_dir = OUTPUT_FILE_BASE + 'counts'
combine_file_dir = OUTPUT_FILE_BASE + 'combined_counts'

if not os.path.isdir(count_file_dir):
    print('Making counts output directory:', count_file_dir)
    subprocess.call(['mkdir', count_file_dir])
else:
    print('BC counts output directory was already made:', count_file_dir)

if DO_COMBINE:
    if not os.path.isdir(combine_file_dir):
        print('Making combined counts output directory:', combine_file_dir)
        subprocess.call(['mkdir', combine_file_dir])
    else:
        print('Combined BC counts output directory was already made:', combine_file_dir)


#
#MAIN CALLS
#

demult_info = get_demultiplex_info(PRIMER_MAP_FILE)
if not DO_COMBINE:
    run_one_sample()

    if len(get_all_libraries_done()) == len(demult_info):
        if not (RUN_COMBINE_WHEN_DONE == 'no'):
            print('Calling', RUN_COMBINE_WHEN_DONE, 'to combine files.')
            subprocess.call(['sbatch', RUN_COMBINE_WHEN_DONE])

#Just a combining run, check all groups for combining
if DO_COMBINE:
    finish_and_get_statistics(RUN_NAME)


