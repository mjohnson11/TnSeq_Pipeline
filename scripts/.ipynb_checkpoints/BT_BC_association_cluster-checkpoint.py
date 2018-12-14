"""
A program to error-correct barcodes and edges in barcoded TnSeq bc association sequencing
Uses a single-bp-deletion-overlaps error correction algorithm (described in "DeletionClusterator")
Outputs both an unfiltered file and a filtered one, where bcs that are associated with multiple edges are excluded
"""

import csv
import pandas as pd
import numpy as np
import Levenshtein
from DeletionErrorCorrector import DeletionClusterator
from milo_tools import reverse_transcribe
from collections import defaultdict

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='input file')
parser.add_argument('unfiltered_output_csv', help='output csv - all bc edge counts')
parser.add_argument('filtered_output_csv', help='output csv - filtered for conflicting edges etc.')
parser.add_argument('excluded_bc_output', help='file with list of excluded barcodes (due to being associated with >1 edge)')
parser.add_argument('-row_column_libraries', action='store_true', help='call if libraries are rows and columns and bcs should be called to wells')

args = parser.parse_args()

# POSITIONAL ARGS
input_file = args.input_file
# output files
unfiltered_output_csv = args.unfiltered_output_csv
filtered_output_csv = args.filtered_output_csv
excluded_bc_output = args.excluded_bc_output
row_column_libraries = args.row_column_libraries

# max edits for errors to cluster
EDITDIST_THRESH_EDGE = 3
EDITDIST_THRESH_BC = 3

# potential centroids with less than or equal to this number of total reads are not included
CENTROID_COUNT_THRESH_EDGE = 20
CENTROID_COUNT_THRESH_BC = 3

def decide_lib(row, which_ones, which_sum):
    # function for calling which library a bc is in - must have >3 reads and >95% of all reads in one column or row
    if row[which_sum] > 0:
        for lib in which_ones:
            if row[lib]/row[which_sum] > 0.95 and row[lib] > 3:
                return lib
    return 'none'

#
# MAIN CALLS
#
dat = pd.read_csv(input_file)
lib_cols = [i for i in dat.keys() if i not in ['BC', 'Edge', 'Total.Counts']]
dat.sort_values(by='Total.Counts', ascending=False, inplace=True)

edge_counts = dat.as_matrix(['Edge', 'Total.Counts'])
use_edges = [i for i in edge_counts if 'CTGTCTCT' not in i[0]]  # getting rid of edges with tagmentation remnants

# Doing error correction on edges and barcodes separately
edge_corrector = DeletionClusterator(use_edges)
edge_corrector.get_deletion_clusters(False)
edge_corrector.cluster_within_delnets(EDITDIST_THRESH_EDGE, CENTROID_COUNT_THRESH_EDGE, False)

bc_corrector = DeletionClusterator(dat.as_matrix(['BC', 'Total.Counts']))
bc_corrector.get_deletion_clusters(False)
bc_corrector.cluster_within_delnets(EDITDIST_THRESH_BC, CENTROID_COUNT_THRESH_BC, False)

rows_total = 0
rows_used = 0
counts_total = 0
counts_used = 0
bc_not_clustered = 0
edge_not_clustered = 0

bc_edge_counts = defaultdict(lambda: ['', '', 0] + [0] * len(lib_cols))
conflicting_bc_to_edge = dict()
conflicting_bcs = set()

# Adding counts based on error correction
for row in dat.as_matrix(['BC', 'Edge', 'Total.Counts'] + lib_cols):
    rows_total += 1
    counts_total += row[2]
    raw_bc = row[0]
    edge = row[1]
    if edge in edge_corrector.corrector:
        real_edge = edge_corrector.corrector[edge]
        if raw_bc in bc_corrector.corrector:
            real_bc = bc_corrector.corrector[raw_bc]
            corrected_row = [real_bc, real_edge] + list(row[2:])
            rows_used += 1
            counts_used += row[2]
            bc_edge = real_bc + '_' + real_edge
            record = bc_edge_counts[bc_edge]  # don't need to check membership b/c of defaultdict
            if record[0] == '':
                record[:2] = corrected_row[:2]
                for i in range(2, len(bc_edge_counts[bc_edge])):
                    bc_edge_counts[bc_edge][i] += corrected_row[i]
        else:
            bc_not_clustered += row[2]
    else:
        edge_not_clustered += row[2]

print('Used', rows_used, 'rows out of', rows_total, 'and', counts_used, 'counts out of', counts_total)
print('Lost', edge_not_clustered, 'counts b/c the edge did not cluster, and', bc_not_clustered, 'b/c bc did not cluster')

with open(unfiltered_output_csv, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['BC', 'Edge', 'Total.Counts'] + lib_cols)
    for bc_edge in sorted([i for i in bc_edge_counts.keys()], key=lambda x: bc_edge_counts[x][2], reverse=True):
        writer.writerow(bc_edge_counts[bc_edge])


# Filtering output so that bcs associated with more than one edge are excluded
excluded_bcs = set()
data_by_bc = dict()
short_to_long_edge = dict()
for bc_edge in sorted([i for i in bc_edge_counts.keys()], key=lambda x: bc_edge_counts[x][2], reverse=True):
    result_row = bc_edge_counts[bc_edge]
    bc, edge, counts = result_row[:3]
    #### Checking for barcodes that are associated with more than one short edge
    if bc in data_by_bc: # this is at least the second bc_edge combo for this bc
        # if there is a second short edge that has >10% the reads of the most common short edge, exclude the 
        # original (top) barcode since we can't be confident about bc association in this case
        top_edge_counts = data_by_bc[bc][2]
        if counts/top_edge_counts > 0.1:
            excluded_bcs.add(bc)
        else:
            # otherwise this bc_edge combo is just excluded from the data, since it is likely from a chimera
            pass
    elif bc not in excluded_bcs:
        # this is the bc edge combination for this bc with the most reads ("top" short edge associated with this bc)
        data_by_bc[bc] = result_row

with open(excluded_bc_output, 'w') as outfile:
    for bc in excluded_bcs:
        outfile.write(bc + '\n')
        
excluded_reads = np.sum([bc_edge_counts[bc_edge][2] for bc_edge in bc_edge_counts if bc_edge.split('_')[0] in excluded_bcs])
print('After filtering for multiple assocations, excluded', len(excluded_bcs), 'bcs, representing', excluded_reads, 'reads.')

# Calling which libraries each BC belongs to based on the filtered data
filtered_data = []
for bc in sorted([i for i in data_by_bc.keys() if i not in excluded_bcs], key=lambda x: data_by_bc[x][2], reverse=True):
   filtered_data.append(data_by_bc[bc])
bc_dat = pd.DataFrame(filtered_data, columns=['BC', 'Edge', 'Total.Counts'] + lib_cols)
bc_dat['bc.rt'] = bc_dat.apply(lambda row: reverse_transcribe(row['BC']), axis=1)

if row_column_libraries:
    # This is more complicated specific case because I sequenced rows and columns and need to identify which well each bc is in
    rows = ['TP-A', 'TP-B', 'TP-C', 'TP-D', 'TP-E',
           'TP-F', 'TP-G', 'TP-H']
    cols = ['TP-1', 'TP-2', 'TP-3', 'TP-4', 'TP-5', 'TP-6',
           'TP-7', 'TP-8', 'TP-9', 'TP-10', 'TP-11', 'TP-12']
    bc_dat['row.sum'] = np.sum(bc_dat[rows], axis=1)
    bc_dat['col.sum'] = np.sum(bc_dat[cols], axis=1)
    bc_dat['row.call'] = bc_dat.apply(lambda row: decide_lib(row, rows, 'row.sum').split('-')[-1], axis=1)
    bc_dat['col.call'] = bc_dat.apply(lambda row: decide_lib(row, cols, 'col.sum').split('-')[-1], axis=1)
    bc_dat['plasmid.lib'] = 'TP-' + bc_dat['row.call'] + bc_dat['col.call']  # specifies a well
else:
    bc_dat['plasmid.lib'] = bc_dat.apply(lambda row: decide_lib(row, lib_cols, 'Total.Counts'), axis=1)

# Exporting filtered output with plasmid library calls, this is the file we will use to parse fitness assay data
bc_dat[['BC', 'Edge', 'Total.Counts', 'plasmid.lib'] + lib_cols].to_csv(filtered_output_csv, index=False)

