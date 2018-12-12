"""
A program to parse single-end BT - (Barcoded TnSeq) reads, counting barcodes

Read the file, for each read

    1) check that the read is good:
        - inline index is correct
        - quality of bc region is above some threshold
        - can extract barcode by regex
    2) record the barcode and unique molecular index (umi, just the first 7 bp), keeping track of counts of each

"""

import gzip
import numpy as np
import csv
import pandas as pd
import os
import subprocess
from collections import defaultdict, Counter
import regex
from milo_tools import FourLineFastq

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('demult_file', help='demultiplex info file')
parser.add_argument('output_base', help='basic output directory')
parser.add_argument('reads_base', help='directory with read files')
parser.add_argument('demult_id_key', help='key to identify which file to parse (which row of the demult file) - for jabba-rays (job arrays)')
parser.add_argument('run_name', help='run name for stats file naming')
parser.add_argument("-quality_cutoff", type=int, default=25, help='quality threshold for bc region')

args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
DEMULT_FILE = args.demult_file
OUTPUT_BASE = args.output_base
if OUTPUT_BASE[-1] != '/':
    OUTPUT_BASE += '/'
READS_BASE = args.reads_base
DEMULT_ID_KEY = int(args.demult_id_key)
QUALITY_CUTOFF = int(args.quality_cutoff)
RUN_NAME = args.run_name

def make_dir(path, name=''):
    if not os.path.isdir(path):
        print('Making', name, ':', path)
        subprocess.call(['mkdir', path])
    else:
        print(name, 'directory was already made:', path)

make_dir(OUTPUT_BASE, 'output')
make_dir(OUTPUT_BASE + 'counts', 'count output')
make_dir(OUTPUT_BASE + 'offtarget_indices', 'offtarget_indices output')
make_dir(OUTPUT_BASE + 'umis', 'UMI output')

BC_REGEX_LIST = [regex.compile('(TGACGC)(\D{28})(TTCTAC)'),
                 regex.compile('(TGACGC)(\D{26,30})(TTCTAC)'),
                 regex.compile('(TGACGC){e<=1}(\D{28})(TTCTAC){e<=1}'),
                 regex.compile('(TGACGC){e<=1}(\D{26,30})(TTCTAC){e<=1}')]

class BcCounter:

    def __init__(self, library_name):
        self.library_name = library_name
        # counts for various outcomes for each read
        self.full_count = 0     # any read
        self.quality_count = 0  # quality score in bc region was too low
        self.inline_count = 0   # inline indices were not right
        self.regex_count = 0    # could not find bc with regex

        # dictionaries with the barcodes as keys like:
        # {bc: [total_counts (int), umis (dict like {umi: count})]
        self.bc_dict = defaultdict(lambda: [0, Counter()])

        # this keeps track of the different inline indices seen, which can provide information about
        # primer cross-contamination
        self.inline_indices_dict = Counter()

    def add_entry(self, r1_seq, r1_qual, bc_start, inline_expected):
        # input is an entry (4 lines from 1 single-end file), time to start checking things...
        # first check: is the inline index correct:
        # note that I am allowing multiple inline indices to be input (separated by ;)
        # (this is because of one library where I screwed up and added two TnF primers... nice.)
        inline = r1_seq[7:15]
        if inline not in inline_expected.split(';'):
            self.inline_count += 1
            self.inline_indices_dict[inline] += 1
        else:
            # Check quality of bc region, continue if above threshold
            if np.mean([ord(c)-33 for c in r1_qual[bc_start:bc_start+28]]) < QUALITY_CUTOFF:
                self.quality_count += 1
            else:
                bc = None
                # regex check, looking for barcode by trying an increasingly lenient series of regular expressions
                for reg_search in BC_REGEX_LIST:
                    reg_hit = reg_search.search(r1_seq[bc_start-8:bc_start+36])
                    if reg_hit:
                        bc = reg_hit.group(2)
                        break
                if not bc:
                    self.regex_count += 1
                else:
                    # extract UMI (first 7 bp)
                    tmp_UMI = r1_seq[:7]
                    self.bc_dict[bc][0] += 1
                    self.bc_dict[bc][1][tmp_UMI] += 1
        # increment total real count
        self.full_count += 1

    def output_counts(self, fout, fout_UMI_fam, masterout):
        # Variable to keep track of the number of UMIs.  Note that since I just have 7 bp UMIs,
        # I expect a certain amount of UMI repeats just by chance (depending the number of counts) 
        # The reason I am counting UMIs anyways is to be able to detect 1st-round pcr bottleneck issues
        UMI_repeat_count = 0
        # output results like BC, Count, UMI.Count (number of unique UMIs seen with this bc)
        outfile1 = open(fout, 'w')
        writer1 = csv.writer(outfile1)
        writer1.writerow(['BC', 'Reads', 'UMI.Count'])
        # output UMI family sizes (number of reads with the same barcode and UMI)
        outfile2 = open(fout_UMI_fam, 'w')
        writer2 = csv.writer(outfile2)
        writer2.writerow(['BC', 'Reads', 'UMI.fam.size.1', 'UMI.fam.size.2', 'etc'])
        for bc_tmp in self.bc_dict.keys():
            tmp_rec = self.bc_dict[bc_tmp]
            # the number of UMI repeats is the total count minus the number of unique UMIs
            UMI_repeat_count += (tmp_rec[0] - len(tmp_rec[1])) 
            writer1.writerow([bc_tmp, tmp_rec[0], len(tmp_rec[1])]) # Simple output
            # getting UMI family sizes
            tmp_dict = dict() # this dict is like {umi family size: number of fams with that size}
            for umi in tmp_rec[1]:
                UMI_fam_size = tmp_rec[1][umi]
                if UMI_fam_size in tmp_dict:
                    tmp_dict[UMI_fam_size] += 1
                else:
                    tmp_dict[UMI_fam_size] = 1
            max_fam_size = max([i for i in tmp_dict.keys()])
            output_row = [bc_tmp, tmp_rec[0]]
            writer2.writerow(output_row + [tmp_dict.setdefault(i, 0) for i in range(1, max_fam_size+1)])
        
        outfile1.close()
        outfile2.close()
        #Print stats
        print([self.library_name, self.full_count, self.inline_count, self.quality_count, self.regex_count,
               UMI_repeat_count])

        # output statistics to a single file (not in a nice order though, will need to sort)
        with open(masterout, 'a') as outfile:
            writer = csv.writer(outfile)
            # there will be many title rows, I will delete all but one when cleaning this up
            writer.writerow(['Library', 'Reads', 'Inline.Index.Excluded', 'Quality.Excluded', 'Regex.Excluded',
                             'UMI.Repeats'])
            writer.writerow([self.library_name, self.full_count, self.inline_count, self.quality_count,
                             self.regex_count, UMI_repeat_count])

    def output_offtarget_inlines(self, inline_out):
        # outputs a csv with the counts for all inline indices seen that did not match the expectation
        with open(inline_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['R1_inline_index', 'Reads'])
            sorted_indices = sorted([i for i in self.inline_indices_dict], key=lambda x: -1*self.inline_indices_dict[x])
            for index in sorted_indices:
                writer.writerow([index, self.inline_indices_dict[index]])


def count_file(dem_info, sample_id, count_output_base, umi_output_base, stats_output, offind_base):
    library_name = dem_info[0]
    inline_index = dem_info[1]
    filenames = dem_info[3]
    # the bc_start region depends on the TnF primer offset:
    bc_start_point = int(dem_info[2]) + 7 + 23
    # initialize bc counting class
    bc_counterator = BcCounter(library_name)
    
    for filename in filenames.split(';'):
        R1 = READS_BASE + filename
        print('Reading file:', R1)
        with gzip.open(R1, 'rt') as r1in:
            for R1_title, R1_seq, R1_qual in FourLineFastq(r1in):
                bc_counterator.add_entry(R1_seq, R1_qual, bc_start_point, inline_index)

    bc_counterator.output_counts(count_output_base + library_name + '_counts.csv',
                                 umi_output_base + library_name + '_umi_fam_sizes.csv',
                                 stats_output)

    bc_counterator.output_offtarget_inlines(offind_base + library_name + '_offtarget_indices.csv')


### BEGIN MAIN CALLS###

# Reading demult file, parsing one file based on demult id key
dem_dat = pd.read_csv(DEMULT_FILE)
rc = 1
for row in dem_dat.as_matrix(['Library', 'TnF_inline', 'TnF_offset', 'Filenames']):
    if rc == DEMULT_ID_KEY:
        count_file(row, rc, OUTPUT_BASE + 'counts/', OUTPUT_BASE + 'umis/', 
                   OUTPUT_BASE + RUN_NAME + '_library_stats.csv', OUTPUT_BASE + 'offtarget_indices/')
    rc += 1 
