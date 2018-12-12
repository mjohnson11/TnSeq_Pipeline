import pandas as pd
import numpy as np
import csv
import os
from scipy import stats as sci_stats
from glob import glob
from collections import defaultdict
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg
import statsmodels.api as sm

NUM_PERMUTATIONS = 10000
NUM_SUBSAMPLES = 10000

# Reading in segregant fitness information
seg_to_fit = {i[0]: i[1] for i in pd.read_csv('Clones_For_Tn96_Experiment.csv').as_matrix(['segregant', 'initial fitness, YPD 30C'])}

# Chromosome lengths (S288C) for seeing how big intervals are 
chromo_lens = {
    'chr01': 230218,
    'chr02': 813184,
    'chr03': 316620,
    'chr04': 1531933,
    'chr05': 576874,
    'chr06': 270161,
    'chr07': 1090940,
    'chr08': 562643,
    'chr09': 439888,
    'chr10': 745751,
    'chr11': 666816,
    'chr12': 1078177,
    'chr13': 924431,
    'chr14': 784333,
    'chr15': 1091291,
    'chr16': 948066
}

def change_well_format(w):
    if '_' in w:
        plate = int(w[1:3])
        t = 'LK' + str(plate) + '-'
        n = int(w.split('_')[1])
        lets = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        l = lets[int(np.floor((n-1)/12))]
        return t + l + str(((n-1) % 12) + 1).zfill(2)
    else:
        return w

    
def get_genotype_dataframe(seg_names):
    # gets genotype matrix, excludes alleles that are at less than 10% either allele
    print('Reading in segregant genotype data')
    d = pd.read_csv('BYxRM_GenoData.csv')
    map_genos = {'B': 0, 'R': 1}
    for w in d.keys():
        if change_well_format(w) in seg_names:
            d[change_well_format(w)] = d[w].map(map_genos)
    d['p'] = np.mean(d[[column for column in d.columns if column != 'marker']], axis=1)
    d_use = d.loc[d['p']>0.1].loc[d['p']<0.9]
    print(len(d) - len(d_use), 'alleles excluded because they were nearly fixed.')
    assert len([s for s in seg_names if s in d.columns]) == len(seg_names)
    return d_use[['marker'] + seg_names]


def calculate_lods_no_weighting(genotype_matrix_centered, geno_mat_std, Ycen, Ystd, Ylen):
    # Modified from Elizabeth Jerison's eLife 2017 scripts
    rs = np.dot(genotype_matrix_centered.T, Ycen)
    rscaled = rs/(geno_mat_std * Ystd * Ylen)
    return -0.5*Ylen*np.log10(1-np.square(rscaled)), np.square(rs)


def refitting_ols(param_names, param_values, phenotypes):
    # ordinary least squares that iteratively excludes non-significant parameters at the 0.05 level
    names = param_names
    all_params = param_values
    model_fit = sm.OLS(phenotypes, sm.add_constant(all_params)).fit()
    sig_indices = [i for i in range(np.shape(all_params)[1]) if model_fit.pvalues[i+1] < 0.05]
    sig_names = [names[i] for i in sig_indices]
    sig_params = all_params[:, sig_indices]
    while np.shape(sig_params)[1] < np.shape(all_params)[1]:
        all_params = sig_params
        names = sig_names
        model_fit = sm.OLS(phenotypes, sm.add_constant(all_params)).fit()
        sig_indices = [i for i in range(np.shape(all_params)[1]) if model_fit.pvalues[i+1] < 0.05]
        sig_names = [names[i] for i in sig_indices]
        sig_params = all_params[:, sig_indices]
    return names, model_fit.rsquared, model_fit.resid, model_fit.params


def detect_qtls(seg_names, geno_mat_df, phenotypes):
    # Modified from Elizabeth Jerison's eLife 2017 scripts
    # Input: list of segregant names, pandas dataframe with segregant names as columns (BY allele = 0, RM allele=1), list of phenotype values
    geno_mat_raw = geno_mat_df.as_matrix(seg_names).T
    p = np.mean(geno_mat_raw, axis=0)
    genotype_mat_centered = (geno_mat_raw - p) / np.sqrt(len(p)*p*(1-p))  # rescaled genotype values for LOD calculation
    # the genotype values are structured so they all have the same std dev : 1/sqrt(num loci)
    genotype_value_std = 1/np.sqrt(len(p))
    r_squared = 0
    one_qtl_r_squared = 0
    n_segs = len(phenotypes)
    QTLs = []
    intervals = []
    all_QTLs_found = False
    while all_QTLs_found ==False:
        #print(QTLs)
        # Fit a linear model using the current QTL list
        if len(QTLs) > 0:
            qtl_matrix = genotype_mat_centered[:,QTLs]
            model_fit = sm.OLS(phenotypes, sm.add_constant(qtl_matrix)).fit()
            residuals = model_fit.resid
            r_squared = model_fit.rsquared
            if len(QTLs) == 1:
                one_qtl_r_squared = r_squared
        else:
            residuals = phenotypes
        resid_cen = residuals - np.mean(residuals)
        resid_std = np.std(resid_cen)
        lods, rvals = calculate_lods_no_weighting(genotype_mat_centered, genotype_value_std, resid_cen, resid_std, n_segs)
        top_rval = np.nanmax(rvals)
        top_lod = np.nanmax(lods)
        top_lod_idx = np.nanargmax(lods)
        # To define the range, we look for a drop in LOD of 2, and make sure it is a real drop by waiting until we see 20 consecutive low LODs
        lb_index = top_lod_idx
        first_consecutive_low_idx = 0  # if it fails, full range
        consecutive_low_scores = 0
        while consecutive_low_scores < 20:
            lb_index -= 1
            if lb_index < 0:
                break
            if lods[lb_index] < top_lod - 2:
                consecutive_low_scores += 1
            else:
                consecutive_low_scores = 0
            if consecutive_low_scores == 1:
                first_consecutive_low_idx = lb_index
        ub_index = top_lod_idx
        first_consecutive_high_idx = len(lods)  # if it fails, full range
        consecutive_low_scores = 0
        while consecutive_low_scores < 20:
            ub_index += 1
            if ub_index >= len(lods):
                break
            if lods[ub_index] < top_lod - 2:
                consecutive_low_scores += 1
            else:
                consecutive_low_scores = 0
            if consecutive_low_scores == 1:
                first_consecutive_high_idx = ub_index
        # Permutations
        permutations = np.array([np.random.permutation(resid_cen) for i in range(NUM_PERMUTATIONS)])
        # calculating dot-products in matrix multiplication form
        permuted_rs = np.matmul(genotype_mat_centered.T, permutations.T)
        permuted_rvals = np.square(permuted_rs)
        bootstrapped_rvals = np.nanmax(permuted_rvals, axis=0)
        sig_threshold = np.sort(bootstrapped_rvals)[int(-0.01 * NUM_PERMUTATIONS)]  # p < .01
        if top_rval > sig_threshold:
            QTLs.append(top_lod_idx)
            intervals.append([first_consecutive_low_idx, first_consecutive_high_idx])
        else:
            all_QTLs_found = True
    # Now fitting by ordinary least squares and dropping qtls that are not significant at the 0.05 level
    if len(QTLs) > 0:
        sig_indices, r_squared, residuals, coeffs = refitting_ols([i for i in range(len(QTLs))], genotype_mat_centered[:,QTLs], phenotypes)
        QTLs = [QTLs[i] for i in sig_indices]
        intervals = [intervals[i] for i in sig_indices]
    if len(QTLs) > 0: 
        # subsampling for r^2 confidence interval
        sub_r2 = []
        for b in range(NUM_SUBSAMPLES):
            indices_chosen = np.random.choice([i for i in range(len(seg_names))], size=int(np.floor(len(seg_names)/2)), replace=False)
            qtl_matrix = genotype_mat_centered[:,QTLs][indices_chosen,:]
            sub_fit = sm.OLS(phenotypes[indices_chosen], sm.add_constant(qtl_matrix)).fit()
            sub_r2.append(sub_fit.rsquared)
        conf_low = np.percentile(sub_r2, 2.5)
        conf_high = np.percentile(sub_r2, 97.5)
    else:
        conf_low, conf_high, r_squared = 0, 0, 0
    return QTLs, intervals, r_squared, conf_low, conf_high, one_qtl_r_squared, residuals 


def ols_fit_and_qtls(seg_names, geno_mat_df, phenotypes, qtl_indices, seg_fitnesses):
    # Modified from Elizabeth Jerison's eLife 2017 scripts
    # Input: list of segregant names, pandas dataframe with segregant names as columns (BY allele = 0, RM allele=1), list of phenotype values
    geno_mat_raw = geno_mat_df.as_matrix(seg_names).T
    p = np.mean(geno_mat_raw, axis=0)
    genotype_mat_centered = (geno_mat_raw - p) / np.sqrt(len(p)*p*(1-p))  # see note on rescaled genotype values 
    # just calculating this because its easy and keeps the lod calculation general, but all these stds are 1/sqrt(num loci) (probably should change)
    geno_mat_stds = np.std(genotype_mat_centered, axis=0) 
    n_segs = len(phenotypes)
    qtl_matrix = genotype_mat_centered[:,qtl_indices]
    all_params = np.column_stack([qtl_matrix, seg_fitnesses])
    sig_indices, r_squared, sig_residuals, coeffs = refitting_ols([i for i in range(len(qtl_indices)+1)], all_params, phenotypes)
    sig_qtl_indices = [qtl_indices[i] for i in sig_indices if i != len(qtl_indices)]
    # I am always including background fitness in the full model - even if it is not significant in the ols
    qtl_matrix = genotype_mat_centered[:,sig_qtl_indices]
    all_params = np.column_stack([qtl_matrix, seg_fitnesses])
    final_fit = sm.OLS(phenotypes, sm.add_constant(all_params)).fit()
    # subsampling for r^2 confidence interval
    sub_r2 = []
    for b in range(NUM_SUBSAMPLES):
        indices_chosen = np.random.choice([i for i in range(len(seg_names))], size=int(np.floor(len(seg_names)/2)), replace=False)
        qtl_matrix = genotype_mat_centered[:,sig_qtl_indices][indices_chosen,:]
        all_params = np.column_stack([qtl_matrix, seg_fitnesses[indices_chosen]])
        sub_fit = sm.OLS(phenotypes[indices_chosen], sm.add_constant(all_params)).fit()
        sub_r2.append(sub_fit.rsquared)
    conf_low = np.percentile(sub_r2, 2.5)
    conf_high = np.percentile(sub_r2, 97.5)
    return final_fit.rsquared, final_fit.resid, conf_low, conf_high, sig_qtl_indices, final_fit.params

def get_interval_size(interval, qtl, markers):
    i1, i2 = interval
    jnk, chromo_q, loc_q, junk, junkie = markers[qtl].split('_')
    jnk, chromo1, loc1, junk, junkie = markers[i1].split('_')
    jnk, chromo2, loc2, junk, junkie = markers[i2].split('_')
    if chromo1 != chromo2: # interval stretches over the chromosome edge, call the interval size
        if chromo_q == chromo1:
            return chromo_lens[chromo_q] - int(loc1)
        else:
            return int(loc2)
    else:
        return int(loc2) - int(loc1)

def qtl_dedup(qtl_infos_list, markers):
    # deduplicating qtls found before and after regressing out fitness
    # if there is a duplicate, we use the qtl with the tighter confidence interval
    overlapping_indices = [[], []]
    qtls1, intervals1 = qtl_infos_list[0]
    qtls2, intervals2 = qtl_infos_list[1]
    for i in range(len(intervals1)):
        for j in range(len(intervals2)):
            if intervals1[i][0] < intervals2[j][0]:
                if intervals2[j][0] < intervals1[i][1]:
                    if get_interval_size(intervals1[i], qtls1[i], markers) < get_interval_size(intervals2[j], qtls2[j], markers):
                        overlapping_indices[0].append([i, j])
                    else:
                        overlapping_indices[1].append([i, j])
            else:
                if intervals1[i][0] < intervals2[j][1]:
                    if get_interval_size(intervals1[i], qtls1[i], markers) < get_interval_size(intervals2[j], qtls2[j], markers):
                        overlapping_indices[0].append([i, j])
                    else:
                        overlapping_indices[1].append([i, j])
    qtls_keep = [qtls1[q] for q in range(len(qtls1)) if q not in [i[0] for i in overlapping_indices[1]]]
    qtls_keep += [qtls2[q] for q in range(len(qtls2)) if q not in [i[1] for i in overlapping_indices[0]]]
    qint_keep = [intervals1[q] for q in range(len(qtls1)) if q not in [i[0] for i in overlapping_indices[1]]]
    qint_keep += [intervals2[q] for q in range(len(qtls2)) if q not in [i[1] for i in overlapping_indices[0]]]
    
    return qtls_keep, qint_keep


def simple_x_regression(background_fits, fit_effects):
    reg_result = sci_stats.linregress(background_fits, fit_effects)
    subsamples = []
    for b in range(NUM_SUBSAMPLES):
        inds_chosen = np.random.choice([i for i in range(len(fit_effects))], size=int(np.floor(len(fit_effects)/2)), replace=False)
        sub_regress = sci_stats.linregress([background_fits[i] for i in inds_chosen], [fit_effects[i] for i in inds_chosen])
        subsamples.append(sub_regress[2]**2)
    return reg_result, subsamples


def analyze_determinants(seg_list, phenos, gm_df):
    sd = dict()
    x_list = np.array([seg_to_fit[seg] for seg in seg_list])
    markers = list(gm_df['marker'])
    regress, sub_r2 = simple_x_regression(x_list, phenos)
    sd['x.slope'] = regress[0]
    sd['x.r2'] = regress[2]**2
    sd['x.p.value'] = regress[3]
    sd['x.r2.95.conf.low'] = np.percentile(sub_r2, 2.5)
    sd['x.r2.95.conf.high'] = np.percentile(sub_r2, 97.5)
    x_residuals = [phenos[i] - (x_list[i]*regress[0]+regress[1]) for i in range(len(phenos))]
    # finding qtls
    qtl_info = detect_qtls(seg_list, gm_df, np.array(phenos))
    sd['qtls.r2'] = qtl_info[2]
    sd['qtls'] = [(markers[qtl_info[0][i]], markers[qtl_info[1][i][0]], markers[qtl_info[1][i][1]]) for i in range(len(qtl_info[0]))]
    sd['qtls.r2.95.conf.low'] = qtl_info[3]
    sd['qtls.r2.95.conf.high'] = qtl_info[4]
    sd['qtls.one.r2'] = qtl_info[5]
    # regressing with background fitness over and above qtl effects (y values are qtl model residuals)
    res_regress, res_sub_r2 = simple_x_regression(x_list, qtl_info[6])
    sd['resid.x.slope'] = res_regress[0]
    sd['resid.x.r2'] = res_regress[2]**2
    sd['resid.x.p.value'] = res_regress[3]
    sd['resid.x.r2.95.conf.low'] = np.percentile(res_sub_r2, 2.5)
    sd['resid.x.r2.95.conf.high'] = np.percentile(res_sub_r2, 97.5)
    # finding qtls over and above qtl effects
    qtl_resid_info = detect_qtls(seg_list, gm_df, np.array(x_residuals))
    sd['resid.qtls.r2'] = qtl_resid_info[2]
    sd['resid.qtls'] = [(markers[qtl_resid_info[0][i]], markers[qtl_resid_info[1][i][0]], markers[qtl_resid_info[1][i][1]]) for i in range(len(qtl_resid_info[0]))]
    sd['resid.qtls.r2.95.conf.low'] = qtl_resid_info[3]
    sd['resid.qtls.r2.95.conf.high'] = qtl_resid_info[4]
    dedup_qtl_info = qtl_dedup([qtl_info[:2], qtl_resid_info[:2]], markers)
    sd['dedup.qtls'] = [(markers[dedup_qtl_info[0][i]], markers[dedup_qtl_info[1][i][0]], markers[dedup_qtl_info[1][i][1]]) for i in range(len(dedup_qtl_info[0]))]
    # fitting full model
    r2, resid, clow, chigh, sig_qtls, coeffs = ols_fit_and_qtls(seg_list, gm_df, np.array(phenos), 
                                                                dedup_qtl_info[0], x_list)
    sig_qtl_tmp_indices = [i for i in range(len(dedup_qtl_info[0])) if dedup_qtl_info[0][i] in sig_qtls]
    sd['full.model.qtls'] = [(markers[dedup_qtl_info[0][i]], markers[dedup_qtl_info[1][i][0]], markers[dedup_qtl_info[1][i][1]]) for i in sig_qtl_tmp_indices]
    sd['full.model.r2'] = r2
    sd['full.model.r2.95.conf.low'] = clow
    sd['full.model.r2.95.conf.high'] = chigh
    sd['full.model.coeffs'] = ';'.join([str(i) for i in coeffs])
    return sd

    
def get_stats_for_one_edge(row, segs, gm_df):
    measured = [seg for seg in segs if row[seg + '.total.cbcs'] >= 4]
    stat_dict = {'num.measured': len(measured)}
    if len(measured) >= 10:
        pvals = [row[seg + '.pval'] for seg in measured]
        pval_sig_boolean = benjamini_hochberg(pvals)[0]  # B/H with alpha=0.05 by default
        sig = [measured[s] for s in range(len(measured)) if pval_sig_boolean[s]]
        stat_dict['num.sig'] = len(sig)
        means = [row[seg + '.mean.s'] for seg in measured]
        variances = [row[seg + '.stderr.s']**2 for seg in measured]
        Vg = np.var(means)
        Ve = np.mean(variances)
        stat_dict['H2'] = (Vg - Ve) / Vg
        # doing n/2 sub-samplings to get error on that
        sub_H2 = []
        seg_indices = [i for i in range(len(measured))]
        for b in range(NUM_SUBSAMPLES):
            segs_chosen = np.random.choice(seg_indices, size=int(np.floor(len(measured)/2)), replace=False)
            Vg = np.var([means[seg_ind] for seg_ind in segs_chosen])
            Ve = np.mean([variances[seg_ind] for seg_ind in segs_chosen])
            sub_H2.append((Vg - Ve) / Vg)
        stat_dict['H2.95.conf.low'] = np.percentile(sub_H2, 2.5)
        stat_dict['H2.95.conf.high'] = np.percentile(sub_H2, 97.5)
        stat_dict.update(analyze_determinants(measured, means, gm_df))
    return stat_dict

def add_analysis(exp, df, output_name):
    o_cols = ['num.measured', 'num.sig', 'x.slope', 'x.p.value', 'resid.x.slope', 'resid.x.p.value', 'full.model.coeffs', 'qtls.one.r2']
    cols_with_conf_ints = ['H2', 'qtls.r2', 'x.r2', 'full.model.r2', 'resid.qtls.r2', 'resid.x.r2']
    full_cols = o_cols + cols_with_conf_ints
    for c in cols_with_conf_ints:
        full_cols += [c + '.95.conf.low', c + '.95.conf.high']
    qtl_cols = ['qtls', 'resid.qtls', 'full.model.qtls']
    exp_segs = [i.split('.')[0] for i in df if '.mean.s' in i]
    geno_df = get_genotype_dataframe(exp_segs)
    edge_stats = {r['Edge']: get_stats_for_one_edge(r, exp_segs, geno_df) for index, r in df.iterrows()}
    for col in full_cols:
        df[col] = df['Edge'].apply(lambda edge: edge_stats[edge].setdefault(col, np.nan))
    for col in qtl_cols:
        df[col] = df['Edge'].apply(lambda edge: '|'.join([';'.join(['_'.join(q.split('_')[1:3]) for q in qtl]) for qtl in edge_stats[edge].setdefault(col, [])])) 
    df.to_csv(output_name, index=False)