"""
OK so the plan:
1) calculate log slopes for all time intervals (0->1, 0->2, 0->3, 0->4, 1->2, 1->3, 1->4, 2->3, 2->4, 3->4)
    only include it if both timepoints have >= some # reads
2) Use the median log slope of all barcodes as the reference - subtract it from all logslopes
3) take the mean and std of the intervals that have enough counts for each barcode
    so I'm just looking at intervals between timepoints with >= some # reads, so for example, I will use
    0->1, 1->2, 2->3, 3->4 if all the timepoints have enough reads, but 1->3, 3->4 if T0 and T2 don't have enough reads
"""

import pandas as pd
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_base', help='input file directory')
parser.add_argument('output_base', help='basic output directory')
parser.add_argument('assays', help='assay names split by commas ')
parser.add_argument('-excluded_bcs_file', type=str, default='none', help='file with bcs to be excluded from analysis, 1 per line')
args = parser.parse_args()

# POSITIONAL ARGS - THESE SHOULD ALWAYS BE THE SAME FOR A RUN
input_base = args.input_base
output_base = args.output_base
assays = args.assays.split(',')
excluded_bcs_file = args.excluded_bcs_file

excluded_bcs = set()
if excluded_bcs_file != 'none':
    with open(excluded_bcs_file, 'r') as infile:
        for line in infile:
            excluded_bcs.add(line.rstrip())

read_cutoff = 10
libs = ['BT-N' + str(i) for i in range(1,11)] + ['BT-H' + str(i) for i in range(1,11)]
tps = [0, 1, 2, 3, 4]

plas_lib_info = pd.read_csv('Plasmid_Lib_Demult.csv')
pl3 = {i[0]: i[1] for i in plas_lib_info.as_matrix(['Lib 3', 'quick number'])}
pl1 = {i[0]: i[1] for i in plas_lib_info.as_matrix(['Lib 1', 'quick number'])}
pl_dict = {'A1': pl1, 'A3': pl3}

def calc_interval_logslope(row, assay, t1, t2, sum1, sum2, cutoff):
    if row['BT-' + assay + '-T' + str(t1)] < cutoff or row['BT-' + assay + '-T' + str(t2)] < cutoff:
        return np.nan
    else:
        logdif = np.log(row['BT-' + assay + '-T' + str(t2)]/sum2) - np.log(row['BT-' + assay + '-T' + str(t1)]/sum1)
        return logdif / (10*(t2-t1))

def get_s_rec(row, assay, cutoff):
    tmp_rec = []
    tps = [i for i in range(5) if row['BT-' + assay + '-T' + str(i)] >= cutoff]
    for i in range(len(tps)-1):
        t1 = tps[i]
        t2 = tps[i+1]
        tmp_rec.append(row['scaled.s_' + str(t1) + '_' + str(t2)])
    return tmp_rec

cols2 = ['mean.s', 'std.s', 'len.s']
for a in assays:
    columns = ['BC', 'Edge', 'Total.Reads'] + ['BT-' + a + '-T' + str(i) for i in range(5)]
    for lib in libs:
        clone = pl_dict[a][lib]
        d_tmp = pd.read_csv(input_base + a + '/' + lib + '_counts.csv')
        d = d_tmp.loc[~d_tmp['BC'].isin(excluded_bcs)]
        for i in range(4):
            for j in range(i+1, 5):
                s1 = sum(d['BT-' + a + '-T' + str(i)])
                s2 = sum(d['BT-' + a + '-T' + str(j)])
                d['s_' + str(i) + '_' + str(j)] = d.apply(lambda r: calc_interval_logslope(r, a, i, j, s1, s2, read_cutoff), axis=1)

        for i in range(4):
            for j in range(i+1, 5):
                median_fit = np.nanmedian(d['s_' + str(i) + '_' + str(j)])
                d['scaled.s_' + str(i) + '_' + str(j)] = d['s_' + str(i) + '_' + str(j)] - median_fit

        d['mean.s'] = d.apply(lambda r: np.nanmean(get_s_rec(r, a, read_cutoff)), axis=1)
        d['std.s'] = d.apply(lambda r: np.nanstd(get_s_rec(r, a, read_cutoff)), axis=1)
        d['len.s'] = d.apply(lambda r: len(get_s_rec(r, a, read_cutoff)), axis=1)

        d_good = d.loc[~np.isnan(d['mean.s'])]
        print(str(clone) + '_' + a + '_' + lib[3:], len(d_tmp), len(d), len(d_good), len(set(d_good['Edge'])))
        d_good[columns + cols2].to_csv(output_base + str(clone) + '_' + a + '_' +
                                       lib[3:] + '_bc_s.csv', index=False)
