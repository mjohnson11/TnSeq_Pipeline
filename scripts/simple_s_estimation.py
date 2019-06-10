import csv
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
np.warnings.filterwarnings('ignore')
from joint_s_inference import BfaParamEstimator
from collections import Counter, defaultdict
import os
import subprocess
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

def make_dir(path, name=''):
    if not os.path.isdir(path):
        print('Making', name, ':', path)
        subprocess.call(['mkdir', path])
    else:
        print(name, 'directory was already made:', path)

# Simple s measurement functions
def calc_interval_logslope(row, assay, t1, t2, sum1, sum2, read_cutoff):
    # simply calculates the slope of log(freq) if boths tps have >=read_cutoff reads (else returns nan)
    if row[assay + '-T' + str(t1)] < read_cutoff or row[assay + '-T' + str(t2)] < read_cutoff:
        return np.nan
    else:
        logdif = np.log(row[assay + '-T' + str(t2)]/sum2) - np.log(row[assay + '-T' + str(t1)]/sum1)
        return logdif / (10*(t2-t1))

def get_s_rec(row, assay, tps, read_cutoff):
    # returns a list of s values from different time intervals for this barcode
    use_tps = [i for i in tps if row[assay + '-T' + str(i)] >= read_cutoff]
    return [row['scaled.s_' + str(use_tps[i]) + '_' + str(use_tps[i+1])] for i in range(len(use_tps)-1)]

def simple_s_measuring(td, assay, tps, neut_edges, neut_reference_read_cutoff, read_cutoff=10):
    # calculates log-slopes for all lineages for all possible time intervals, scales by the median log-slope of the neutral class and then for
    # each lineages, using all timepoints with at least read_cutoff reads, averages the consecutive time interval scaled s values to get a lineage s ('s'). 
    # Also records the number of intervals used ('len.s'), and the standard error of the different measures ('stderr.s')
    td['Total.Reads'] = np.sum(td[[assay + '-T' + str(t) for t in tps]], axis=1)
    for i in range(len(tps)-1):
        for j in range(i+1, len(tps)):
            s1 = sum(td[assay + '-T' + str(tps[i])])
            s2 = sum(td[assay + '-T' + str(tps[j])])
            td['s_' + str(tps[i]) + '_' + str(tps[j])] = td.apply(lambda r: calc_interval_logslope(r, assay, tps[i], tps[j], s1, s2, read_cutoff), axis=1)
            neut_data = td.loc[td['Edge'].isin(neut_edges)].loc[td['Total.Reads'] > neut_reference_read_cutoff]
            median_neut_fit = np.nanmedian(neut_data['s_' + str(tps[i]) + '_' + str(tps[j])])
            td['scaled.s_' + str(tps[i]) + '_' + str(tps[j])] = td['s_' + str(tps[i]) + '_' + str(tps[j])] - median_neut_fit
    # calling get_s_rec 3 times is silly but hey I'm doing it...
    td['s'] = td.apply(lambda r: np.nanmean(get_s_rec(r, assay, tps, read_cutoff)), axis=1) # mean s across timepoint intervals
    td['len.s'] = td.apply(lambda r: len([s for s in get_s_rec(r, assay, tps, read_cutoff) if not np.isnan(s)]), axis=1) # this is number of timepoint intervals
    td['stderr.s'] = td.apply(lambda r: np.nanstd(get_s_rec(r, assay, tps, read_cutoff))/np.sqrt(r['len.s']), axis=1) # this is stderr across timepoint intervals

# Helper functions
def get_bc_s_vals(td, s_column):
    # returns a dictionary like {Edge: [median s, mean s, num bcs, stderr of mean s]} (all in reference to the set of bcs corresponding to this edge)
    edge_grouped = td[['Edge', s_column]].groupby('Edge', as_index=False)
    edge_bc_stats = {i[0]: [np.nanmedian(i[1][s_column]), # median s for each edge
                            np.nanmean(i[1][s_column]),  # mean s for each edge
                            len(i[1].loc[~np.isnan(i[1][s_column])]),  # number of bcs for each edge
                            np.nanstd(i[1][s_column])/np.sqrt(len(i[1].loc[~np.isnan(i[1][s_column])]))] for i in edge_grouped} # stderr s for each edge
    return edge_bc_stats

def reduce_bcs(td, tp_names, bc_target):
    """
    combines barcode counts randomly (really alphabetically) to have some maximum number of barcodes. Columns should be ['BC', 'Edge'] and then read 
    count columns (tp_names) (bc column is actually not necessary). Returns a dataframe with columns like ['cBC', 'Edge'] + tp_names where 
    cBC is like edge_C where C ranges from 0 to bc_target (or from 0 to the number of original bcs-1, in which case they are just renamed)
    """
    td_shuf = td.sort_values(by='BC')  # sorts rows alphabetically, so they are random with respect to total counts
    edge_bc_counter = Counter()
    dm = td_shuf.as_matrix(['Edge'] + tp_names)
    # adding an extra columns that just counts up for bcs of the same edge
    new_rows = []
    for row in dm:
        edge_bc_counter[row[0]] += 1
        new_rows.append(list(row) + [edge_bc_counter[row[0]]])

    new_df = pd.DataFrame(new_rows, columns=(['Edge'] + tp_names + ['edge_count']))  # back to pandas world
    new_df['cBC'] = new_df.apply(lambda row: str(row['Edge']) + '_' + str(row['edge_count'] % bc_target), axis=1)  # make cBC column
    combined = new_df[['cBC'] + tp_names].groupby('cBC', as_index=False).sum()  # combine bcs with same cBC
    combined['Edge'] = combined['cBC'].str.split('_').str[0]  # remake Edge column
    return combined

# Log-likelihood ratio exclusion functions
def outlier_ll_test(tmp_R, gen_times):
    # Usese a maximum likelihood model based on reads as multinomial draws (see BfaParamEstimator class) to compare two hypotheses:
    # 1) Two lineages have the same fitness, 2) They have different fitnesses
    # tmp_R should be a numpy array of read counts from two lineages
    assert len(tmp_R) == 2
    if min(np.sum(tmp_R, axis=1)) < 40:  # if either lineage has low reads the test doesn't work well, so we ignore them when one lineage has total reads < 40
        return 0
    else:
        inf_neut = BfaParamEstimator(tmp_R.T, [], [], gen_times)
        inf_neut.run_ll_max(max_iters=10, neutral_assumption=True) # This should actually one take one iteration, since there is only one free parameter
        inf = BfaParamEstimator(tmp_R.T, [], [], gen_times)
        inf.run_ll_max(max_iters=1000)  # 1000 iterations should be more than enough for 2 params
        return inf_neut.current_ll - inf.current_ll # log-likelihood ratio - higher = more likely that the two lineages have different s
    
def add_ll_column(td, tps, tp_names, s_col, ll_col_name, near_median_s_buffer=0.01):
    # Adding a log-likelihood ratio to bcs that quantifies whether they look like outliers from other barcodes corresponding to the same mutation/edge
    # We define the neutral class of bcs for each edge as all the lineages with s measured within near_median_s_buffer of the median s for the edge
    # we combined the counts from these bcs to make a reference lineage.  Then we do the log-likelihood ratio test on all bcs not in this class,
    # testing how likely it is that they have an s different from the reference class (just to be clear this is only a reference within the edge)
    td_edge_grouped = td.groupby('Edge', as_index=False)
    bc_to_ll_dif = dict()
    for e in td_edge_grouped:
        num_bcs = len(e[1])
        if num_bcs > 2:  # There must be at least 3 bcs to even try this exclusion mehod
            median_s = np.nanmedian(e[1][s_col])
            near_median = e[1].loc[e[1][s_col] > (median_s - near_median_s_buffer)].loc[e[1][s_col] < (median_s + near_median_s_buffer)]
            if len(near_median) > 0:
                e_outside = e[1].loc[~e[1]['BC'].isin(near_median['BC'])]
                median_rec = near_median[tp_names].sum().as_matrix(tp_names)
                for entry in e_outside.as_matrix(['BC'] + tp_names):
                    bc_to_ll_dif[entry[0]] = outlier_ll_test(np.array([list(median_rec), list(entry[1:])]), [i*10 for i in tps])
    td[ll_col_name] = td.apply(lambda row: bc_to_ll_dif.setdefault(row['BC'], 0), axis=1)


# Main running function
def s_estimation(segs, rep_info, output_base, input_base, experiment, ll_cutoff, bc_target_num, neut_edges, exclusion_file=None, neut_ref_read_cutoff=30, tp_read_count_cutoff=5000, assay_to_first_tp=None, bc_out_extra_cols=[], consider_all_edges_neutral=False):
    if exclusion_file:
        ex_out = open(exclusion_file, 'w')
        ex_writer = csv.writer(ex_out)
    
    for seg in segs:
        print('Working on', seg)
        make_dir(output_base + seg, seg)
        rep_count = 1
        bc1_out_cols = ['BC', 'Edge', 's', 'LL.ratio'] + bc_out_extra_cols
        edge_out_cols = ['Edge', 's', 'stderr.s', 'len.s']
        edge_to_bc_s = defaultdict(lambda: defaultdict(list))
        for rep in rep_info[seg]:
            assay, plasmid_lib = rep.split('.')
            if assay_to_first_tp is None:
                ftp = 0
            else:
                ftp = assay_to_first_tp[assay]  # for an experiment where I am excluding timepoint zero
            d_tmp = pd.read_csv(input_base + assay + '/' + experiment + '-' + plasmid_lib + '_counts.csv')
            d_tmp['Edge'] = d_tmp.Edge.astype(str)  # casting to string because interpretting simulation results as ints causes errors
            tp_names = [tp for tp in [assay + '-T' + str(t) for t in range(ftp, 5)] if np.sum(d_tmp[tp]) > tp_read_count_cutoff]
            d_tmp['Total.Reads'] = np.sum(d_tmp[tp_names], axis=1)
            d = d_tmp.loc[d_tmp['Total.Reads'] > 0]  # Taking this step because some bcs have zero reads after tps are excluded
            tps = [int(tp.split('-T')[-1]) for tp in tp_names]
            if consider_all_edges_neutral:
                neut_edges = list(d['Edge'])
            neut_data = d.loc[d['Edge'].isin(neut_edges)].loc[d['Total.Reads'] > neut_ref_read_cutoff]
            if len(neut_data) > 4 and len(tps) > 2:
                simple_s_measuring(d, assay, tps, neut_edges, neut_ref_read_cutoff)
                add_ll_column(d, tps, tp_names, 's', 'LL.ratio')
                d[bc1_out_cols + tp_names].to_csv(output_base + seg + '/' + seg + '_rep' + str(rep_count) + '_bc_s_initial.csv', index=False)
                d_cut_pre = d.loc[d['LL.ratio'] < ll_cutoff][['BC', 'Edge'] + tp_names]
                d_cut = reduce_bcs(d_cut_pre, tp_names, bc_target_num)
                simple_s_measuring(d_cut, assay, tps, neut_edges, neut_ref_read_cutoff)
                for entry in d_cut.as_matrix(['Edge', 's']):
                    if not np.isnan(entry[1]):
                        edge_to_bc_s[entry[0]][rep].append(entry[1])
                d_cut[['cBC', 'Edge', 's', 'stderr.s', 'len.s'] + tp_names].to_csv(output_base + seg + '/' + seg + '_rep' + str(rep_count) + '_bc_s.csv', index=False)
                d_edge = d_cut[['Edge'] + tp_names].groupby("Edge", as_index=False).sum()
                simple_s_measuring(d_edge, assay, tps, neut_edges, neut_ref_read_cutoff)
                bc_stats = get_bc_s_vals(d_cut, 's')
                d_edge['bc_stats'] = d_edge.apply(lambda row: bc_stats.setdefault(row['Edge'], [np.nan, np.nan, np.nan, np.nan]), axis=1)
                d_edge['median.bc.s'], d_edge['mean.bc.s'], d_edge['num.cbcs'], d_edge['bc.stderr.s'] = d_edge['bc_stats'].str  # .str splits lists somehow...
                d_edge[edge_out_cols + ['median.bc.s', 'mean.bc.s', 'num.cbcs', 'bc.stderr.s'] + tp_names].to_csv(output_base + seg + '/' + seg + '_rep' + str(rep_count) + '_edge_s.csv', index=False)
            elif len(neut_data) < 5:
                if exclusion_file:
                    ex_writer.writerow([seg, rep, 'Insufficient reference/neutral barcodes'])
                print(seg, rep, 'excluded by neut data')
            else:
                if exclusion_file:
                    ex_writer.writerow([seg, rep, 'Low sequencing coverage'])
                print(seg, rep, 'excluded by low read counts')
            rep_count += 1

        with open(output_base + seg + '/' + seg + '_edge_s.csv', 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Edge', 'mean.s', 'stderr.s', 'total.cbcs', 'pval', 'rep1.cbcs', 'rep1.s', 'rep1.stderr.s', 
                             'rep2.cbcs', 'rep2.s', 'rep2.stderr.s'])
            neut_bc_s = defaultdict(list)
            for edge in neut_edges:
                if edge in edge_to_bc_s:
                    for rep in edge_to_bc_s[edge]:
                        neut_bc_s[rep] += edge_to_bc_s[edge][rep]
            for edge in edge_to_bc_s:
                tmp = edge_to_bc_s[edge]
                usable_reps = [r for r in tmp if len(tmp[r]) >= 1 and len(neut_bc_s[r]) >= 2]  # must have at least 1 cbc for the replicate
                num_cbcs = [len(tmp[r]) for r in usable_reps]
                if len(usable_reps) > 0 and np.sum(num_cbcs) > 3: # must have at least 3 cbcs overall
                    # I get the mean and standard error from inverse variance weighting
                    means = np.array([np.mean(tmp[r]) for r in usable_reps])
                    all_s = []
                    # What we really want for standard errors is the standard error of the difference between s for the edge and the neutral edges
                    # The standard error of a difference is the square root of the sum of the two errors
                    std_err_d = dict()
                    for r in usable_reps:
                        all_s += tmp[r]
                        if len(tmp[r]) > 1:
                            std_err_d[r] = np.sqrt((np.std(tmp[r], ddof=1)/np.sqrt(len(tmp[r])))**2 + (np.std(neut_bc_s[r], ddof=1)/np.sqrt(len(neut_bc_s[r])))**2)
                        else: # if there is only one barcode, we use the standard deviation of bc s in the other replicate as an estimate of the standard error
                            other_rep = [rep for rep in usable_reps if rep != r][0] # since we require 3 cbcs, the other rep must have at least 2
                            std_err_d[r] = np.sqrt(np.std(tmp[other_rep], ddof=1)**2 + (np.std(neut_bc_s[r], ddof=1)/np.sqrt(len(neut_bc_s[r])))**2)
                    std_errs = np.array([std_err_d[r] for r in usable_reps])
                    # inverse variance averaging code modified from Venkataram et al. 2016 code
                    use_mean = np.sum(means*np.power(std_errs, -2))/np.sum(np.power(std_errs,-2))
                    # standard error is based on deviations of each bc s from this mean
                    use_std_err = np.sqrt(np.sum([(s - use_mean)**2 for s in all_s])/(len(all_s)-1))/np.sqrt(len(all_s))
                    # for significance testing, I fit the cbc s data by ordinary least squares with replicate as a fixed effect predictor
                    # and ask whether the fitness effect is different from zero (for each mutation)
                    mat = []
                    for rep in usable_reps:
                        for val in tmp[rep]:
                            mat.append([rep, val]) 
                    td = pd.DataFrame(mat, columns=['replicate', 's'])
                    # sig testing
                    ps = list(ols('s ~ replicate', data=td).fit().pvalues)
                    tmp_row = [edge, use_mean, use_std_err, np.sum(num_cbcs), ps[0]] # first p val is for teh intercept - is the mutation's effect non-zero
                    for r in rep_info[seg]:
                        if r in usable_reps:
                            tmp_row += [len(tmp[r]), np.mean(tmp[r]), std_err_d[r]]
                        else:
                            tmp_row += [np.nan, np.nan, np.nan]

                    writer.writerow(tmp_row)