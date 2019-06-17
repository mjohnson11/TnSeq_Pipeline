import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from collections import Counter
from joint_s_inference import BfaParamEstimator
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('experiment', help='experiment name')
parser.add_argument('demult_id_key', help='key to identify which file to parse (which row of the demult file) - for jabba-rays (job arrays)')
args = parser.parse_args()

experiment = args.experiment
DEMULT_ID_KEY = int(args.demult_id_key)

#parameters
tp_read_count_cutoff = 5000 # cutoff # reads in a timepoint
read_cutoff = 10 # for a set of timepoints to be used in the simple s calculation, they both have to have >= this # reads 
neut_reference_read_cutoff = 30 # to be used to calculate median fitness, barcodes from reference mutations need to have more than this # total reads
near_median_s_buffer = 0.01 # range on either side of the median bc s to use as the reference for finding outliers

# Simple s measurement functions
def calc_interval_logslope(row, assay, t1, t2, sum1, sum2):
    if row[assay + '-T' + str(t1)] < read_cutoff or row[assay + '-T' + str(t2)] < read_cutoff:
        return np.nan
    else:
        logdif = np.log(row[assay + '-T' + str(t2)]/sum2) - np.log(row[assay + '-T' + str(t1)]/sum1)
        return logdif / (10*(t2-t1))

def get_s_rec(row, assay, tps):
    tmp_rec = []
    tmp_noise = []
    use_tps = [i for i in tps if row[assay + '-T' + str(i)] >= read_cutoff]
    for i in range(len(use_tps)-1):
        tmp_rec.append(row['scaled.s_' + str(use_tps[i]) + '_' + str(use_tps[i+1])])
        tmp_noise.append(row['interval.noise.term_' + str(use_tps[i]) + '_' + str(use_tps[i+1])])
    weighted_mean = np.sum(np.array(tmp_rec)*np.power(np.array(tmp_noise), -2)) / np.sum(np.power(np.array(tmp_noise), -2))
    return tmp_rec, weighted_mean


def simple_s_measuring(td, assay, tps):
    td['Total.Reads'] = np.sum(td[[assay + '-T' + str(t) for t in tps]], axis=1)
    for i in range(len(tps)-1):
        for j in range(i+1, len(tps)):
            s1 = sum(td[assay + '-T' + str(tps[i])])
            s2 = sum(td[assay + '-T' + str(tps[j])])
            td['s_' + str(tps[i]) + '_' + str(tps[j])] = td.apply(lambda r: calc_interval_logslope(r, assay, tps[i], tps[j], s1, s2), axis=1)
            if experiment == 'TP':
                neut_data = td.loc[td['Edge'] < 5].loc[td['Total.Reads'] > neut_reference_read_cutoff]
            else:
                neut_data = td.loc[td['Total.Reads'] > neut_reference_read_cutoff]            
            median_neut_fit = np.nanmedian(neut_data['s_' + str(tps[i]) + '_' + str(tps[j])])
            td['scaled.s_' + str(tps[i]) + '_' + str(tps[j])] = td['s_' + str(tps[i]) + '_' + str(tps[j])] - median_neut_fit
            # noise term is 1 / (sqrt(geom mean of reads at t1 and t2) * delta t)
            td['interval.noise.term_' + str(tps[i]) + '_' + str(tps[j])] = np.power(td[assay + '-T' + str(tps[i])] * td[assay + '-T' + str(tps[j])], -0.25) / (tps[j]-tps[i])
    td['simple.s'] = td.apply(lambda r: get_s_rec(r, assay, tps)[1], axis=1)
    td['simple.len.s'] = td.apply(lambda r: len([s for s in get_s_rec(r, assay, tps)[0] if not np.isnan(s)]), axis=1)
    td['simple.stderr.s'] = td.apply(lambda r: np.nanstd(get_s_rec(r, assay, tps)[0])/np.sqrt(r['simple.len.s']), axis=1)

# Helper functions
def get_median_s_vals(td, s_column):
    edge_grouped = td[['Edge', s_column]].groupby('Edge', as_index=False)
    edge_bc_stats = {i[0]: [np.nanmedian(i[1][s_column]), # median s for each edge
                            len(i[1].loc[~np.isnan(i[1][s_column])]),  # number of bcs for each edge
                            np.nanstd(i[1][s_column])/np.sqrt(len(i[1].loc[~np.isnan(i[1][s_column])]))] for i in edge_grouped} # stderr s for each edge
    return edge_bc_stats

def outlier_ll_test(tmp_R, gen_times):
    #td should be a 2 bc (2 row) dataframe
    assert len(tmp_R) == 2
    if min(np.sum(tmp_R, axis=1)) < 40:
        return 0
    else:
        inf_neut = BfaParamEstimator(tmp_R.T, [], [], gen_times)
        inf_neut.run_ll_max(max_iters=10, neutral_assumption=True)
        inf = BfaParamEstimator(tmp_R.T, [], [], gen_times)
        inf.run_ll_max(max_iters=1000)
        return inf_neut.current_ll - inf.current_ll
    
def add_ll_column(td, tp_names, s_col, ll_col_name):
    td_edge_grouped = td.groupby('Edge', as_index=False)
    bc_to_ll_dif = Counter()
    for e in td_edge_grouped:
        num_bcs = len(e[1])
        if num_bcs > 2:
            median_s = np.nanmedian(e[1][s_col])
            near_median = e[1].loc[e[1][s_col] > (median_s - near_median_s_buffer)].loc[e[1][s_col] < (median_s + near_median_s_buffer)]
            if len(near_median) > 0:
                e_outside = e[1].loc[~e[1]['BC'].isin(near_median['BC'])]
                median_rec = near_median[tp_names].sum().as_matrix(tp_names)
                for entry in e_outside.as_matrix(['BC'] + tp_names):
                    bc_to_ll_dif[entry[0]] = outlier_ll_test(np.array([list(median_rec), list(entry[1:])]), [i*10 for i in tps])
    td[ll_col_name] = td.apply(lambda row: bc_to_ll_dif.setdefault(row['BC'], 0), axis=1)


bc1_out_cols = ['BC', 'Edge', 'simple.s', 'LL.ratio.simple', 'Edge.s', 'Used.s']
# Doing it now
for sn in range(DEMULT_ID_KEY*5 + 1, (DEMULT_ID_KEY+1)*5 + 1):
    for rep in ['r1', 'r2']:
        ftp = 0
        d_tmp = pd.read_csv('sim/' + experiment + '-sim' + str(sn) + '_' + rep + '_counts.csv')
        tp_names = [tp for tp in ['sim-T' + str(t) for t in range(ftp, 5)] if np.sum(d_tmp[tp]) > tp_read_count_cutoff]
        d_tmp['Total.Reads'] = np.sum(d_tmp[tp_names], axis=1)
        d = d_tmp.loc[d_tmp['Total.Reads'] > 0]  # Taking this step because some bcs have zero reads after tps are excluded
        tps = [int(tp.split('-T')[-1]) for tp in tp_names]
        if experiment == 'TP':
            neut_data = d.loc[d['Edge'] < 5].loc[d['Total.Reads'] > neut_reference_read_cutoff]
        else:
            neut_data = d.loc[d['Total.Reads'] > neut_reference_read_cutoff]
        if len(neut_data) > 4 and len(tps) > 2:
            simple_s_measuring(d, 'sim', tps)
            add_ll_column(d, tp_names, 'simple.s', 'LL.ratio.simple')
            d[bc1_out_cols + tp_names].to_csv('sim_s_v5/' + experiment + '_' + str(sn) + '_' + rep + '_bc_s_initial.csv', index=False)
        elif len(neut_data) < 5:
            print(sn, rep, 'excluded by neut data')
        else:
            print(sn, rep, 'excluded by low read counts')




