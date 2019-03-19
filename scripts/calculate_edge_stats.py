import pandas as pd
import numpy as np
import csv
import os
from scipy import stats as sci_stats
from glob import glob
from collections import defaultdict
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg
import statsmodels.api as sm
import statsmodels.formula.api as smf
import time

otime = time.time()

NUM_PERMUTATIONS = 10000

# Reading in segregant fitness information
seg_to_fit = {i[0]: i[1] for i in pd.read_csv('../accessory_files/Clones_For_Tn96_Experiment.csv').as_matrix(['segregant', 'initial fitness, YPD 30C'])}

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
    d = pd.read_csv('../accessory_files/BYxRM_GenoData.csv')
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


def center_geno_mat(g, seg_names):
    # input pandas dataframe with segregant names as columns (BY allele = 0, RM allele=1), list of segregant names
    geno_mat_raw = g.as_matrix(seg_names).T  # this works even when seg_names has duplicates for replicates (thanks pandas!)
    p = np.mean(geno_mat_raw, axis=0)
    assert len(p) == np.shape(geno_mat_raw)[1]
    return (geno_mat_raw - p) / np.sqrt(len(p)*p*(1-p))  # rescaled genotype values for LOD calculation


def detect_qtls(seg_names, genotype_mat_centered, phenotypes):
    # Modified from Elizabeth Jerison's eLife 2017 scripts
    # Input: list of segregant names, centered genotype matrix, list of phenotype values
    # the genotype values are structured so they all have the same std dev : 1/sqrt(num loci)
    genotype_value_std = 1/np.sqrt(np.shape(genotype_mat_centered)[1])
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
    return QTLs, intervals


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


def fit_and_subsample(use_formula, tmp_df, num_subsamples):
    main_fit = smf.ols(formula = use_formula, data = tmp_df).fit()
    sub_r2 = []
    for b in range(num_subsamples):
        indices_chosen = np.random.choice([i for i in range(len(tmp_df))], size=int(np.floor(len(tmp_df)/2)), replace=False)
        sub_fit = smf.ols(formula = use_formula, data = tmp_df.iloc[indices_chosen]).fit()
        sub_r2.append(sub_fit.rsquared)
    return main_fit, np.percentile(sub_r2, 2.5), np.percentile(sub_r2, 97.5)


def analyze_determinants(seg_list, phenos, gm_df, num_subsamples):
    # In this function, I fit a series of nested models to explain fitness measurements
    # I both fit models to the residuals from restricted models (i.e. fit a qtl model to the residuals from a background fitness model)
    # And do F tests to compare models and restricted models.  The first is more conservative, but the p values are similar in both cases
    sd = {'var': np.var(phenos)}
    markers = list(gm_df['marker'])
    gm_centered = center_geno_mat(gm_df, seg_list)
    # Making a dataframe for modeling
    df = pd.DataFrame(seg_list, columns=['segregant'])
    df['s'] = phenos
    df['x'] = df['segregant'].apply(lambda s: seg_to_fit[s])
    qtl_info = detect_qtls(seg_list, gm_centered, np.array(phenos))
    if len(qtl_info[0]) > 0:
        df[['locus_' + markers[i] for i in qtl_info[0]]] = pd.DataFrame(gm_centered[:,qtl_info[0]], index=df.index) # adding qtls detected to the dataframe
        qtl_model_tmp = smf.ols(formula='s ~ ' + ' + '.join(['locus_' + markers[i] for i in qtl_info[0]]), data=df).fit()
        df['s_qtl_resids'] = qtl_model_tmp.resid
    else:
        df['s_qtl_resids'] = df['s']
    # this is a bit redundant, but I'm just quickly fitting these models so I can get the residuals to use in other models
    x_model_tmp = smf.ols(formula='s ~ x', data=df).fit()
    df['s_x_resids'] = x_model_tmp.resid
    # more qtl fitting
    qtl_resid_info = detect_qtls(seg_list, gm_centered, np.array(x_model_tmp.resid))
    dedup_qtl_info = qtl_dedup([qtl_info, qtl_resid_info], markers)
    new_qtls = [i for i in qtl_resid_info[0] if i not in qtl_info[0]]
    all_qtl_indices = list(set(qtl_info[0] + qtl_resid_info[0]))
    df[['locus_' + markers[i] for i in new_qtls]] = pd.DataFrame(gm_centered[:,new_qtls], index=df.index) # adding new qtls detected to the dataframe
    # fitting the full model (redundantly, to get residuals)
    full_model_tmp = smf.ols(formula='s ~ ' + ' + '.join(['x'] + ['locus_' + markers[i] for i in dedup_qtl_info[0]]), data=df).fit()
    df['s_full_resids'] = full_model_tmp.resid
    models = {
        'segregant': 's ~ segregant',
        'x': 's ~ x',
        'qtl': 's ~ ' + ' + '.join(['locus_' + markers[i] for i in qtl_info[0]]),
        'resid_qtl': 's_x_resids ~ ' + ' + '.join(['locus_' + markers[i] for i in qtl_resid_info[0]]),
        'resid_x': 's_qtl_resids ~ x',
        'full': 's ~ ' + ' + '.join(['x'] + ['locus_' + markers[i] for i in dedup_qtl_info[0]]),
        'full_plus_seg': 's ~ ' + ' + '.join(['segregant', 'x'] + ['locus_' + markers[i] for i in dedup_qtl_info[0]]),
        'full_resid_seg': 's_full_resids ~ segregant'
    }
    fits = dict()
    for model in models:
        if models[model][-2:] == '~ ':  # no qtls found, no point fitting a model
            sd[model + '_model_r2'], sd[model + '_model_p'], sd[model + '_model_r2_95_conf_low'], sd[model + '_model_r2_95_conf_high'] = np.nan, np.nan, np.nan, np.nan
            sd[model + '_model_p_values'], sd[model + '_model_params'], sd[model + '_model_coeffs'] = 'NA', 'NA', 'NA'
        else:
            fits[model], sd[model + '_model_r2_95_conf_low'], sd[model + '_model_r2_95_conf_high'] = fit_and_subsample(models[model], df, num_subsamples)
            sd[model + '_model_p'] = fits[model].f_pvalue
            sd[model + '_model_r2'] = fits[model].rsquared
            sd[model + '_model_p_values'] = ';'.join([str(i) for i in fits[model].pvalues])
            sd[model + '_model_params'] = models[model][models[model].index('~')+2:].replace(' + ', ';')
            sd[model + '_model_coeffs'] = ';'.join([str(i) for i in fits[model].params]) # the first parameter is the intercept

    # Adding some nice columns
    sd['x_slope'] = fits['x'].params[1]
    sd['full_model_x_slope'] = fits['full'].params[1]
    # This is an unscaled effect size that says how much the regression predicts the difference in s should be btwn the most and least fit segregant
    sd['full_model_x_effect_size_measure'] = sd['full_model_x_slope'] * (np.max(df['x'])-np.min(df['x']))  
    # Getting the qtl coefficients in the full model so they correspond to the difference in mean fitness effect between alleles
    qtl_effects = []
    for i in range(len(dedup_qtl_info[0])):
        high = np.max(df['locus_' + markers[dedup_qtl_info[0][i]]])
        low = np.min(df['locus_' + markers[dedup_qtl_info[0][i]]])
        qtl_effects.append((high-low)*fits['full'].params[i+2])
    sd['full_model_qtl_effect_sizes'] = ';'.join([str(ef) for ef in qtl_effects])

    # model comparison using F test
    model_comps = [('full', 'x'), ('full', 'qtl'), ('full_plus_seg', 'full')]
    for (m1, m2) in model_comps:
        if sd[m1 + '_model_p_values'] != 'NA' and sd[m2 + '_model_p_values'] != 'NA':
            f_jnk, sd['model_comp_p_' + m1 + '_vs_' + m2], df_diff_jnk = fits[m1].compare_f_test(fits[m2])
        else:
            if m2 == 'qtl' and len(qtl_info[0]) == 0:
                sd['model_comp_p_full_vs_qtl'] = sd['x_model_p']  # no qtls so no model comparison
            else:
                sd['model_comp_p_' + m1 + '_vs_' + m2] = np.nan

    sd['qtls'] = [(markers[qtl_info[0][i]], markers[qtl_info[1][i][0]], markers[qtl_info[1][i][1]]) for i in range(len(qtl_info[0]))]
    sd['resid.qtls'] = [(markers[qtl_resid_info[0][i]], markers[qtl_resid_info[1][i][0]], markers[qtl_resid_info[1][i][1]]) for i in range(len(qtl_resid_info[0]))]
    sd['full.model.qtls'] = [(markers[dedup_qtl_info[0][i]], markers[dedup_qtl_info[1][i][0]], markers[dedup_qtl_info[1][i][1]]) for i in range(len(dedup_qtl_info[0]))]
    return sd


def get_stats_for_one_edge(row, segs, gm_df, num_subsamples, use_only_two_rep_segs):
    # The difference here is I am including replicate s measurements in the model instead of using the mean s
    if use_only_two_rep_segs:  #only include segregants w s measured in both replicates
        measured = [seg for seg in segs if row[seg + '.rep1.cbcs'] >= 2 and row[seg + '.rep2.cbcs'] >= 2]
    else:
        measured = [seg for seg in segs if row[seg + '.rep1.cbcs'] >= 2 or row[seg + '.rep2.cbcs'] >= 2]
    reps = ['rep1', 'rep2']
    measured_by_rep = {r: [seg for seg in segs if row[seg + '.' + r + '.cbcs'] >= 2] for r in reps}
    stat_dict = {'num.measured': len(measured)}
    if len(measured) > 0:
        pvals = [row[seg + '.pval'] for seg in measured]
        pval_sig_boolean = benjamini_hochberg(pvals)[0]  # B/H with alpha=0.05 by default
        sig = [measured[s] for s in range(len(measured)) if pval_sig_boolean[s]]
        stat_dict['num.sig'] = len(sig)
        stat_dict['avg_s'] = np.nanmean([row[seg + '.mean.s'] for seg in measured])
    if len(measured) >= 10:
        means = [row[seg + '.mean.s'] for seg in measured]
        variances = [row[seg + '.stderr.s']**2 for seg in measured]
        Vg = np.var(means)
        Ve = np.mean(variances)
        stat_dict['H2'] = (Vg - Ve) / Vg
        # doing n/2 sub-samplings to get error on that
        sub_H2 = []
        seg_indices = [i for i in range(len(measured))]
        for b in range(num_subsamples):
            segs_chosen = np.random.choice(seg_indices, size=int(np.floor(len(measured)/2)), replace=False)
            Vg = np.var([means[seg_ind] for seg_ind in segs_chosen])
            Ve = np.mean([variances[seg_ind] for seg_ind in segs_chosen])
            sub_H2.append((Vg - Ve) / Vg)
        stat_dict['H2_95_conf_low'] = np.percentile(sub_H2, 2.5)
        stat_dict['H2_95_conf_high'] = np.percentile(sub_H2, 97.5)
        rep_means = [row[seg + '.rep1.s'] for seg in measured_by_rep['rep1']] + [row[seg + '.rep2.s'] for seg in measured_by_rep['rep2']]
        stat_dict.update(analyze_determinants(measured_by_rep['rep1']+measured_by_rep['rep2'], rep_means, gm_df, num_subsamples))
    return stat_dict

def add_analysis(exp, df, output_name, num_subsamples, use_only_two_rep_segs=False):
    full_cols = ['num.measured', 'num.sig', 'H2', 'H2_95_conf_low', 'H2_95_conf_high',
                 'model_comp_p_full_vs_qtl', 'model_comp_p_full_vs_x', 'model_comp_p_full_plus_seg_vs_full', 'var',
                 'avg_s', 'x_slope', 'full_model_x_slope', 'full_model_x_effect_size_measure', 'full_model_qtl_effect_sizes']
    mods = ['segregant', 'x', 'qtl', 'resid_qtl', 'resid_x', 'full', 'full_plus_seg', 'full_resid_seg']
    suffixes = ['_model_r2', '_model_p', '_model_r2_95_conf_low', '_model_r2_95_conf_high', '_model_p_values', '_model_params', '_model_coeffs']
    for c in mods:
        full_cols += [c + s for s in suffixes]
    qtl_cols = ['qtls', 'resid.qtls', 'full.model.qtls']
    exp_segs = [i.split('.')[0] for i in df if '.mean.s' in i]
    geno_df = get_genotype_dataframe(exp_segs)
    edge_stats = {r['Edge']: get_stats_for_one_edge(r, exp_segs, geno_df, num_subsamples, use_only_two_rep_segs) for index, r in df.iterrows()}
    for col in full_cols:
        df[col] = df['Edge'].apply(lambda edge: edge_stats[edge].setdefault(col, np.nan))
    for col in qtl_cols:
        df[col] = df['Edge'].apply(lambda edge: '|'.join([';'.join(['_'.join(q.split('_')[1:3]) for q in qtl]) for qtl in edge_stats[edge].setdefault(col, [])]))
    df.to_csv(output_name, index=False)
    print('Time:', time.time()-otime)
