import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from glob import glob
import numpy as np
from calculate_edge_stats import add_analysis
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('num_subsamples', help='number of subsamplings to get confidence intervals on r^2')
parser.add_argument('output_file', help='output_file')
parser.add_argument('-only_use_two_seg_reps', action='store_true')
args = parser.parse_args()

NUM_SUBSAMPLES = int(args.num_subsamples)
output_file = args.output_file
only_use_two_seg_reps = args.only_use_two_seg_reps


edge_info_file = '../accessory_files/Tn96_edges_chosen_final.csv'
ann_info_file = '../../Mutation_Annotation/Edges_annotated.csv'

# Reading extra info
edge_info = pd.read_csv(edge_info_file)[['Edge', 'Is.Reference']]
edge_info['Long.Edge'] = edge_info['Edge']
edge_info['Type'] = edge_info['Is.Reference'].apply(lambda r: {True: 'Reference', False: 'Experiment'}[r])
edge_info['Edge'] = edge_info['Long.Edge'].str[:15]

# these were edges added to the libraries for comparison to a previous experiment (Kryazhimskiy et al. 2014),
# but are not included in the final analysis for this paper
control_edges = ['CTGATTTGTGCTGTC', 'GCTGCTTATGAGGAT', 'CTAAGTGTGAAGGAG'] 
accidentally_included = ['TTAATTGCACTTATG'] # this edge is in the experiment due to a mistake during re-arraying or multiple clones in a single well

mat = [[e, 'Extra_control'] for e in control_edges] + [[e, 'Extra_mistake'] for e in accidentally_included]

td = pd.concat([edge_info[['Edge', 'Type']], pd.DataFrame(mat, columns=['Edge', 'Type'])], ignore_index=True)

# ADDING ANNOTATION INFO
ann = pd.read_csv(ann_info_file)
ann['Edge'] = ann['Edge'].str[:15]
ann_use = ann[['Edge.ID', 'Edge', 'Gene.Use', 'Gene', 'Gene.nearby', 'ORF', 'ORF.nearby',
               'briefDescription', 'briefDescription.nearby', 'chromosome',
               'description', 'description.nearby', 'end', 'end.nearby',
               'functionSummary', 'functionSummary.nearby', 'insertion_edge',
               'insertion_strand', 'orf_strand', 'orf_strand.nearby',
               'phenotypeSummary', 'phenotypeSummary.nearby', 'start', 'start.nearby']]
edge_dat = td.merge(ann_use, on='Edge', how='left')

# ADDING FITNESS INFO
cols = ['mean.s', 'stderr.s', 'total.cbcs', 'pval', 'rep1.s', 'rep1.cbcs', 'rep1.stderr.s', 'rep2.s', 'rep2.cbcs', 'rep2.stderr.s']
TP_segs = [i.split('/')[-1] for i in glob('../../S_Estimation/TP_output/*')]
for s in TP_segs:
    ed = {i[0]: i[1:] for i in pd.read_csv('../../S_Estimation/TP_output/' + s + '/' + s + '_edge_s.csv').as_matrix(['Edge'] + cols)}
    for c in range(len(cols)):
        edge_dat[s + '.' + cols[c]] = edge_dat['Edge'].apply(lambda e: ed.setdefault(e, [np.nan for i in range(len(cols))])[c])

# ADDING GENETIC DETERMINANT ANALYSIS
add_analysis('TP', edge_dat, output_file, NUM_SUBSAMPLES, use_only_two_rep_segs=only_use_two_seg_reps)