"""
A program to combine and cluster BT barcodes

1) combines barcodes

2) error-corrects the barcodes to known barcodes, splits into multiple files based on plasmid libraries

3) adds edges based on a BC-edge association file

"""

import csv
import pandas as pd
import os
import subprocess
from milo_tools import reverse_transcribe
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('assay_file', help='file containing assay library info')
parser.add_argument('output_base', help='basic output directory')
parser.add_argument('bc_assoc_file', help='bc-edge association file')
parser.add_argument('assay_id_key', help='assay id key to identify which assay to work on (which row in assay_file), for jabba-rays (job arrays)')
parser.add_argument('-excluded_bcs_file', type=str, default='none', help='file with bcs to be excluded from analysis, 1 per line')

args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
ASSAY_FILE = args.assay_file
OUTPUT_BASE = args.output_base
if OUTPUT_BASE[-1] != '/':
    OUTPUT_BASE += '/'
BC_ASSOC_FILE = args.bc_assoc_file
ASSAY_ID_KEY = int(args.assay_id_key)

# BC exclustion is included because in the BT1 BFA I sequenced the barcodes of a few clones in separate libraries in the same lane
# and because of some minor sequencing prep contamination or illumina index swapping, the frequency of these barcodes are very 
# jumpy in the data, which doesn't represent the reality.  Excluding these barcodes doesn't affect estimates of s in any way
# (in the same way that missing the bc->plasmid-library association for some barcodes is ok)
excluded_bcs_file = args.excluded_bcs_file
excluded_bcs = set()
if excluded_bcs_file != 'none':
    with open(excluded_bcs_file, 'r') as infile:
        for line in infile:
            excluded_bcs.add(line.rstrip())
    
def get_deletion_neighborhood(stringer):
    # returns a set of all single character deletions of a string (includes the string)
    return set([stringer] + [stringer[:x] + stringer[x+1:] for x in range(len(stringer))])
    
def find_bc_match(del_dict, b):
    # finds the correct barcode by looking for overlap with the deletion network
    # this accounts for all single errors - substitutions, insertions, and deletions
    if b in del_dict:
        return del_dict[b]
    else:
        hits = set()
        for b_edit in get_deletion_neighborhood(b):
            if b_edit in del_dict:
                hits.add(del_dict[b_edit])
        if len(hits) == 1:
            return hits.pop()
        elif len(hits) == 0:
            return 'None found'
        else:
            return 'Multiple hits'
        
    
## Reading information about barcodes and their associated edges and libraries ##
bc_dat = pd.read_csv(BC_ASSOC_FILE)
plasmid_libs = [i for i in set(bc_dat['plasmid.lib']).union({'Unclustered'}) if 'none' not in i]
bc_dat['bc.rt'] = bc_dat.apply(lambda row: reverse_transcribe(row['BC']), axis=1)

bc_del_dict = dict() # dict points from all single bp deletions of each barcode to the barcode 
bc_info_dict = dict() # points from the barcode to its edge and library
for row in bc_dat.as_matrix(['bc.rt', 'Edge', 'plasmid.lib']):
    if row[0] not in excluded_bcs:
        if 'none' not in row[2]:
            assert row[2] in plasmid_libs
            bc_del_dict.update({d: row[0] for d in get_deletion_neighborhood(row[0])})
            bc_info_dict[row[0]] = [row[1], row[2]]
        else:
            bc_info_dict[row[0]] = [row[1], 'Unclustered']


## GETTING LIBRARY INFORMATION FOR COMBINING LIBRARIES ##
with open(ASSAY_FILE, 'r') as infile:
    reader = csv.reader(infile)
    rc = 1
    for row in reader:
        if rc == ASSAY_ID_KEY:
            run = row[0]
            libs = row[1:] # the "libs" are the files for timepoints in the assay (run)
        rc += 1
        
if not os.path.isdir(OUTPUT_BASE + run):
    print('Making output directory:', OUTPUT_BASE + run)
    subprocess.call(['mkdir', OUTPUT_BASE + run])
else:
    print('Output directory was already made:', OUTPUT_BASE + run)


## CLUSTERING AND COMBINING READS ETC. ##        
        
# Main data structure - bc_dict - nested dict like:
# {Plasmid Library: 
# {BC: {'BC': bc, 'Edge': edge, 'Total.Reads’: total_reads, Lib1': reads, 'Lib1.UMI.Count': UMI_counts...}
# }
bc_dict = {pl: dict() for pl in plasmid_libs} 
column_names = ['BC', 'Edge', 'Total.Reads'] + libs + [i + '.UMI.Count' for i in libs] 
# Stats by library: {Library: {‘Total.Reads’: total_reads, ‘Total.BCs’: total_bcs, 
# ‘plasmid_lib1.Reads’, p_lib1_reads, 'plasmid_lib1.BCs'; p_lib1_bcs, etc.}
stats_d = dict()

# the "libs" are the files for timepoints in the assay
for lib in libs:
    stats_d[lib] = {'Total.Reads': 0, 'BCs.Clustering.Twice': 0}
    stats_d[lib].update({pl + '.Reads': 0 for pl in plasmid_libs})
    stats_d[lib].update({pl + '.BCs': set() for pl in plasmid_libs})
    infile = open(OUTPUT_BASE + 'counts/' + lib + '_counts.csv', 'r')
    reader = csv.reader(infile)
    rc = 0
    for row in reader:
        rc += 1
        if rc > 1:
            bc = row[0]
            reads = int(row[1])
            umis = int(row[2])
            stats_d[lib]['Total.Reads'] += reads
            ## BARCODE ERROR CORRECTION ##
            if bc in bc_info_dict:  # Check if the barcode is one of the expected barcodes (no errors)
                correct_bc = bc
                edge, plasmid_lib = bc_info_dict[correct_bc]
            else:
                bc_look = find_bc_match(bc_del_dict, bc)  # using deletion neighborhoods of known barcodes
                if bc_look == 'Multiple hits':  
                    # This barcode had a deletion neighborhood that matched with 2 real bcs
                    # It will be assigned Unclustered
                    stats_d[lib]['BCs.Clustering.Twice'] += 1
                    correct_bc = bc
                    edge, plasmid_lib = 'NA', 'Unclustered'
                elif bc_look == 'None found':
                    # No real bc found, assigned Unclustered
                    correct_bc = bc
                    edge, plasmid_lib = 'NA', 'Unclustered'
                else:
                    # BC successfully corrected 
                    correct_bc = bc_look
                    edge, plasmid_lib = bc_info_dict[correct_bc]
            ## ADDING TO DATA STRUCTURES ##
            stats_d[lib][plasmid_lib + '.Reads'] += reads
            stats_d[lib][plasmid_lib + '.BCs'].add(correct_bc)
            if correct_bc not in bc_dict[plasmid_lib]:
                bc_dict[plasmid_lib][correct_bc] = {'BC': correct_bc, 'Edge': edge, 'Total.Reads': reads}
            else:
                bc_dict[plasmid_lib][correct_bc]['Total.Reads'] += reads
            if lib in bc_dict[plasmid_lib][correct_bc]:
                bc_dict[plasmid_lib][correct_bc][lib] += reads
                bc_dict[plasmid_lib][correct_bc][lib + '.UMI.Count'] += umis
            else:
                bc_dict[plasmid_lib][correct_bc][lib] = reads
                bc_dict[plasmid_lib][correct_bc][lib + '.UMI.Count'] = umis
                
    infile.close()              
    stats_d[lib]['Total.BCs'] = rc - 1
    stats_d[lib].update({pl + '.BCs': len(stats_d[lib][pl + '.BCs']) for pl in plasmid_libs})


## OUTPUTTING COUNT FILES FOR EACH PLASMID LIBRARY ##   
for pl in plasmid_libs: 
    d = bc_dict[pl] # dict of barcode counts for this plasmid library
    sorted_bcs = sorted(d.keys(), key=lambda x: -1*d[x]['Total.Reads']) # sort by total reads
    with open(OUTPUT_BASE + run + '/' + pl + '_counts.csv', 'w') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(column_names)
        for bc in sorted_bcs:
            writer_row = []
            for colname in column_names:
                writer_row.append(d[bc].setdefault(colname, 0))
            writer.writerow(writer_row)
            

## OUTPUTTING STATS ##
with open(OUTPUT_BASE + run + '/' + run + '_error_correction_statistics.csv', 'w') as outfile:
    writer = csv.writer(outfile)
    titles = (['Total.Reads', 'Total.BCs', 'BCs.Clustering.Twice'] + 
              [pl + '.Reads' for pl in plasmid_libs] + [pl + '.BCs' for pl in plasmid_libs])
    writer.writerow(['Library'] + titles)
    for lib in libs:
        writer.writerow([lib] + [stats_d[lib].setdefault(c, 0) for c in titles])