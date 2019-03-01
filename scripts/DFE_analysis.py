import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from glob import glob
import numpy as np
from scipy import stats as sci_stats
from statsmodels.stats.multitest import fdrcorrection as benjamini_hochberg
from collections import defaultdict
from calculate_edge_stats import analyze_determinants, get_genotype_dataframe

NUM_SUBSAMPLES = 1000

tp_all = pd.read_csv('../../Analysis/TP_data_by_edge.csv')
tp = tp_all.loc[tp_all['Type']=='Experiment']
bt = pd.read_csv('../../Analysis/BT_data_by_edge.csv')
bt['Long.Edge'] = bt['Edge']
bt['Edge'] = bt['Long.Edge'].str[:15]
dats = {'BT': bt, 'TP': tp}

dfe_cols = ['background.fitness', 'mean', 'median', 'variance', 'skew', 'kurtosis', 'significant.beneficial.mutations', 'significant.deleterious.mutations',
           'mean.sub.low', 'mean.sub.high', 'median.sub.low', 'median.sub.high', 'variance.sub.low', 'variance.sub.high', 'skew.sub.low',
           'skew.sub.high', 'kurtosis.sub.low', 'kurtosis.sub.high']

cols_to_analyze = ['background.fitness', 'mean', 'median', 'variance', 'skew', 'kurtosis', 'significant.beneficial.mutations', 'significant.deleterious.mutations']

seg_to_fit = {i[0]: i[1] for i in pd.read_csv('../accessory_files/Clones_For_Tn96_Experiment.csv').as_matrix(['segregant', 'initial fitness, YPD 30C'])}

dfes = defaultdict(lambda: defaultdict(dict))
dfes_sig = defaultdict(lambda: defaultdict(dict))
exps = {'BT': 'MM', 'TP': 'FM'}
exp_segs = {exp: [i.split('.')[0] for i in dats[exp].columns if '.mean.s' in i] for exp in dats}
for exp in exps:
    segs = exp_segs[exp]
    d = dats[exp]
    for seg in segs:
        measured = d.loc[d[seg + '.total.cbcs'] >= 4]
        dfes[exp][seg] = list(measured[seg + '.mean.s'])
        pvals = list(measured[seg + '.pval'])
        sig = measured.loc[benjamini_hochberg(pvals)[0]] # B/H with alpha=0.05 by default
        dfes_sig[exp][seg] = list(sig[seg + '.mean.s'])

for exp in exps:
    td = dfes[exp]
    td_sig = dfes_sig[exp]
    segs = [s for s in td if len(td[s]) > 50]
    tmp_dict = dict()
    for seg in segs:
        sub_means = [np.nanmean(np.random.choice(td[seg], size=int(len(td[seg])/2), replace=False)) for i in range(NUM_SUBSAMPLES)]
        sub_medians = [np.nanmedian(np.random.choice(td[seg], size=int(len(td[seg])/2), replace=False)) for i in range(NUM_SUBSAMPLES)]
        sub_variances = [np.nanvar(np.random.choice(td[seg], size=int(len(td[seg])/2), replace=False)) for i in range(NUM_SUBSAMPLES)]
        sub_skews = [sci_stats.skew(np.random.choice(td[seg], size=int(len(td[seg])/2), replace=False)) for i in range(NUM_SUBSAMPLES)]
        sub_kurtosis = [sci_stats.kurtosis(np.random.choice(td[seg], size=int(len(td[seg])/2), replace=False)) for i in range(NUM_SUBSAMPLES)]
        tmp_dict[seg] = {
            'background.fitness': seg_to_fit[seg],
            'mean': np.nanmean(td[seg]),
            'median': np.nanmedian(td[seg]),
            'variance': np.nanvar(td[seg]),
            'skew': sci_stats.skew(td[seg]),
            'kurtosis': sci_stats.kurtosis(td[seg]),
            'significant.beneficial.mutations': len([i for i in td_sig[seg] if i > 0]),
            'significant.deleterious.mutations': len([i for i in td_sig[seg] if i < 0]),
            'mean.sub.low': np.percentile(sub_means, 2.5),
            'mean.sub.high': np.percentile(sub_means, 97.5),
            'median.sub.low': np.percentile(sub_medians, 2.5),
            'median.sub.high': np.percentile(sub_medians, 97.5),
            'variance.sub.low': np.percentile(sub_variances, 2.5),
            'variance.sub.high': np.percentile(sub_variances, 97.5),
            'skew.sub.low': np.percentile(sub_skews, 2.5),
            'skew.sub.high': np.percentile(sub_skews, 97.5),
            'kurtosis.sub.low': np.percentile(sub_kurtosis, 2.5),
            'kurtosis.sub.high': np.percentile(sub_kurtosis, 97.5)
        }
    rows = [[c] + [tmp_dict[s][c] for s in segs] for c in dfe_cols]
    df = pd.DataFrame(rows, columns=['DFE.statistic'] + segs)
    full_cols = ['var']
    mods = ['segregant', 'x', 'qtl', 'resid_qtl', 'resid_x', 'full', 'full_plus_seg', 'full_resid_seg']
    suffixes = ['_model_r2', '_model_p', '_model_r2_95_conf_low', '_model_r2_95_conf_high', '_model_p_values', '_model_params', '_model_coeffs']
    for c in mods:
        full_cols += [c + s for s in suffixes]
    qtl_cols = ['qtls', 'resid.qtls', 'full.model.qtls']
    exp_segs = [i.split('.')[0] for i in df if '.mean.s' in i]
    geno_df = get_genotype_dataframe(segs)
    determinant_stats = defaultdict(dict)
    determinant_stats.update({r['DFE.statistic']: analyze_determinants(segs, list(r[segs]), geno_df, NUM_SUBSAMPLES) for index, r in df.loc[df['DFE.statistic'].isin(cols_to_analyze)].iterrows()})
    for col in full_cols:
        df[col] = df['DFE.statistic'].apply(lambda edge: determinant_stats[edge].setdefault(col, np.nan))
    for col in qtl_cols:
        df[col] = df['DFE.statistic'].apply(lambda edge: '|'.join([';'.join(['_'.join(q.split('_')[1:3]) for q in qtl]) for qtl in determinant_stats[edge].setdefault(col, [])])) 
    df.to_csv('../../Analysis/' + exp + '_DFE_statistics.csv', index=False)
