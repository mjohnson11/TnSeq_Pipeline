import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from glob import glob
import numpy as np
from calculate_edge_stats import add_analysis

NUM_SUBSAMPLES = 1 # no subsampling for confidence intervals r^2 of each model, bc tehre are many mutations and few segregants

ann_info_file = '../../Mutation_Annotation/Edges_annotated.csv'

# ADDING ANNOTATION INFO
ann = pd.read_csv(ann_info_file)
edge_dat = ann[['Edge.ID', 'Edge', 'Gene.Use', 'Gene_ORF', 'Gene_ORF.nearby',
               'briefDescription', 'briefDescription.nearby', 'chromosome',
               'description', 'description.nearby', 'end', 'end.nearby',
               'functionSummary', 'functionSummary.nearby', 'insertion_edge',
               'insertion_strand', 'orf_strand', 'orf_strand.nearby',
               'phenotypeSummary', 'phenotypeSummary.nearby', 'start', 'start.nearby']]

# ADDING FITNESS INFO
cols = ['mean.s', 'stderr.s', 'total.cbcs', 'pval', 'rep1.s', 'rep1.cbcs', 'rep1.stderr.s', 'rep2.s', 'rep2.cbcs', 'rep2.stderr.s']
BT_segs = [i.split('/')[-1] for i in glob('../../S_Estimation/BT_output/*')]
for s in BT_segs:
    ed = {i[0]: i[1:] for i in pd.read_csv('../../S_Estimation/BT_output/' + s + '/' + s + '_edge_s.csv').as_matrix(['Edge'] + cols)}
    for c in range(len(cols)):
        edge_dat[s + '.' + cols[c]] = edge_dat['Edge'].apply(lambda e: ed.setdefault(e, [np.nan for i in range(len(cols))])[c])

# ADDING GENETIC DETERMINANT ANALYSIS        
add_analysis('BT', edge_dat, '../../Analysis/BT_data_by_edge.csv', NUM_SUBSAMPLES)
