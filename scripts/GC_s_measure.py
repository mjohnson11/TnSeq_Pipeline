import pandas as pd
import numpy as np
import csv
from scipy import stats as sci_stats
pd.options.mode.chained_assignment = None  # default='warn'

# From looking at the growth curves, I can see that there is some issue with segregant LK3-D09 rep1, 
# which saturates way lower than the rest - excluded
segs = [
    'LK1-G09',
     'LK1-B04',
     #'LK3-D09',
     'LK5-B01',
     'LK4-B01',
     'LK2-A12',
     'LK2-D07',
     'LK6-A05',
     'LK1-G11',
     'LK1-F05',
     'LK2-B04',
     'LK1-A05'
]

# Functions for simple s estimation, very similar to in simple_s_estimation.py, just differences in column naming etc.
def calc_interval_logslope(row, tname1, tname2, t1, t2, sum1, sum2, read_cutoff):
    # simply calculates the slope of log(freq) if boths tps have >=read_cutoff reads (else returns nan)
    if row[tname1] < read_cutoff or row[tname2] < read_cutoff:
        return np.nan
    else:
        logdif = np.log(row[tname2]/sum2) - np.log(row[tname1]/sum1)
        return logdif / (t2-t1)
    
def get_s_rec(row, tp_names, read_cutoff):
    # returns a list of s values from different time intervals for this barcode
    use_tps = [i for i in tp_names if row[i] >= read_cutoff]
    return [row['scaled.s_' + use_tps[i] + '_' + use_tps[i+1]] for i in range(len(use_tps)-1)]

def simple_s_measuring_gc(td, tps, tp_names, neut_edges, neut_reference_read_cutoff, s_name, read_cutoff=10):
    td['Total.Reads'] = np.sum(td[tp_names], axis=1)
    for i in range(len(tps)-1):
        for j in range(i+1, len(tps)):
            s1 = sum(td[tp_names[i]])
            s2 = sum(td[tp_names[j]])
            td['s_' + tp_names[i] + '_' + tp_names[j]] = td.apply(lambda r: calc_interval_logslope(r, tp_names[i], tp_names[j], tps[i], tps[j], s1, s2, read_cutoff), axis=1)
            neut_data = td.loc[td['Edge'].isin(neut_edges)].loc[td['Total.Reads'] > neut_reference_read_cutoff]
            median_neut_fit = np.nanmedian(neut_data['s_' + tp_names[i] + '_' + tp_names[j]])
            td['scaled.s_' + tp_names[i] + '_' + tp_names[j]] = td['s_' + tp_names[i] + '_' + tp_names[j]] - median_neut_fit
    td[s_name] = td.apply(lambda r: np.nanmean(get_s_rec(r, tp_names, read_cutoff)), axis=1)


# Read general info
edge_info_file = '../accessory_files/Tn96_edges_chosen_final.csv'
edge_info = pd.read_csv(edge_info_file)
neut_edges = [i[:15] for i in edge_info.loc[edge_info['new.num.sig']==0]['Edge']]

# Reading in the clustering information, which will allow us to calculate the cell density for each segregant at each timepoint
gcf = pd.read_csv('../../BT_Bioinformatic_Work/GC_BFA_output/GCF/GCF_error_correction_statistics.csv')
gcp = pd.read_csv('../../BT_Bioinformatic_Work/GC_BFA_output/GCP/GCP_error_correction_statistics.csv')
segreps = pd.read_csv('../accessory_files/GC_segregant_replicate_info.csv')
r2s = list(segreps['Replicate_2'])
rs = list(segreps['Replicate_1']) + list(segreps['Replicate_2'])
seg_to_reps = {i[0]: i[1:] for i in segreps.as_matrix(['segregant', 'Replicate_1', 'Replicate_2'])}
gcp['Time'] = gcp['Library'].apply(lambda r: int(r.split('-')[1]))
gcf['Time'] = gcf['Library'].apply(lambda r: int(r.split('-')[1]))

cell_counts = pd.read_csv('../accessory_files/GC_cell_count_estimates.csv')
cell_counts['Library'] = cell_counts.apply(lambda r: r['Experiment'].split('-')[0] + '-' + str(r['Time']), axis=1)
gcf_new = gcf.merge(cell_counts[['Library', 'cells_ml_ave']], on='Library', how='left')
gcp_new = gcp.merge(cell_counts[['Library', 'cells_ml_ave']], on='Library', how='left')

for lib in [i for i in gcf if i in ['TP-'+j+'.Reads' for j in r2s]]:
    gcf_new[lib.split('.')[0][3:] + '.freq'] = gcf_new[lib] / (gcf_new['Total.Reads'] - gcf_new['Unclustered.Reads'])
    gcf_new[lib.split('.')[0][3:] + '.cells'] = gcf_new[lib.split('.')[0][3:] + '.freq'] * gcf_new['cells_ml_ave'] * 12
    
for lib in [i for i in gcp if i in ['TP-'+j+'.Reads' for j in rs]]:
    gcp_new[lib.split('.')[0][3:] + '.freq'] = gcp_new[lib] / (gcp_new['Total.Reads'] - gcp_new['Unclustered.Reads'])
    gcp_new[lib.split('.')[0][3:] + '.cells'] = gcp_new[lib.split('.')[0][3:] + '.freq'] * gcp_new['cells_ml_ave'] * 24

gcp_renamer = {seg_to_reps[seg][0] + '.cells': seg + '_r1.cells' for seg in segs}
gcp_renamer.update({seg_to_reps[seg][1] + '.cells': seg + '_r2.cells' for seg in segs})
gcp_new.rename(index=str, columns=gcp_renamer)[['Library', 'Time'] + [s+'_r1.cells' for s in segs] + [s+'_r2.cells' for s in segs]].to_csv('../../Analysis/GCP_cell_density_estimates.csv', index=False)

gcf_renamer = {seg_to_reps[seg][1] + '.cells': seg + '_r2.cells' for seg in segs}
gcf_new.rename(index=str, columns=gcf_renamer)[['Library', 'Time'] + [s+'_r2.cells' for s in segs]].to_csv('../../Analysis/GCF_cell_density_estimates.csv', index=False)

# Getting the exponential growth rate for each segregant in each assay (log(cell density) / time in hours)
s_rec = dict()
times = [4, 7, 10, 13]
with open('GC_seg_growth_rates.csv', 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['segregant', 'plate_expo_rate_r1', 'plate_expo_rate_r2', 'flask_expo_rate'])
    for seg in segs:
        r1, r2 = seg_to_reps[seg]
        lr = sci_stats.linregress(times, [np.log(list(gcp_new.loc[gcp_new['Time']==t][r1+'.cells'])[0]) for t in times])
        lr2 = sci_stats.linregress(times, [np.log(list(gcp_new.loc[gcp_new['Time']==t][r2+'.cells'])[0]) for t in times])
        lrf = sci_stats.linregress(times, [np.log(list(gcf_new.loc[gcf_new['Time']==t][r2+'.cells'])[0]) for t in times])
        s_rec[seg] = {'Pexpo1': lr[0], 'Pexpo2': lr2[0], 'Fexpo': lrf[0]}
        writer.writerow([seg, lr[0], lr2[0], lrf[0]])
        
# Only using barcodes that we saw in the TP experiment and were not excluded as outliers, I combine bc counts into edge counts:
gcf_tps = list(gcf_new['Library'])
gcp_tps = list(gcp_new['Library'])
dats = dict()
for seg in segs:
    r1_lib, r2_lib = seg_to_reps[seg]
    r1_tp = pd.read_csv('../../S_Estimation/TP_output/' + seg + '/' + seg + '_rep1_bc_s_initial.csv')
    r1_bcs_ok = list(r1_tp.loc[r1_tp['LL.ratio']<15]['BC'])
    r2_tp = pd.read_csv('../../S_Estimation/TP_output/' + seg + '/' + seg + '_rep2_bc_s_initial.csv')
    r2_bcs_ok = list(r2_tp.loc[r2_tp['LL.ratio']<15]['BC'])
    flask = pd.read_csv('../../BT_Bioinformatic_Work/GC_BFA_output/GCF/TP-' + r2_lib + '_counts.csv')[['BC', 'Edge'] + gcf_tps]
    p1 = pd.read_csv('../../BT_Bioinformatic_Work/GC_BFA_output/GCP/TP-' + r1_lib + '_counts.csv')[['BC', 'Edge'] + gcp_tps]
    p2 = pd.read_csv('../../BT_Bioinformatic_Work/GC_BFA_output/GCP/TP-' + r2_lib + '_counts.csv')[['BC', 'Edge'] + gcp_tps]
    fe = flask.loc[flask['BC'].isin(r2_bcs_ok)][['Edge'] + gcf_tps].groupby('Edge', as_index=False).sum()
    p1e = p1.loc[p1['BC'].isin(r1_bcs_ok)][['Edge'] + gcp_tps].groupby('Edge', as_index=False).sum()
    p2e = p2.loc[p2['BC'].isin(r2_bcs_ok)][['Edge'] + gcp_tps].groupby('Edge', as_index=False).sum()
    td = p1e.merge(p2e, on='Edge', how='outer', suffixes=('_r1', '_r2')).fillna(0)
    dats[seg] = td.merge(fe, on='Edge', how='outer').fillna(0)
    
# measuring s
# * During exponential phase, hours 4-13, using the estimated cell densities to define the number of generations 
# * During the lag transition, from 0 to 7, 24 to 31, 29 to 36, again using estimated cell densities to define the number of generations, for comparison to the exponential phase measurement
# * During saturation, from 20-29, not defining generations, just making the measurement over a series of hours

for seg in dats:
    td = dats[seg]
    lib1 = list(segreps.loc[segreps['segregant']==seg]['Replicate_1'])[0]
    lib2 = list(segreps.loc[segreps['segregant']==seg]['Replicate_2'])[0]
    times_to_cells1 = {t: list(gcp_new.loc[gcp_new['Library']=='GCP-'+str(t)][lib1 + '.cells'])[0] for t in 
                       [0, 4, 7, 10, 13, 16, 20, 24, 29, 31, 36]}
    times_to_cells2 = {t: list(gcp_new.loc[gcp_new['Library']=='GCP-'+str(t)][lib2 + '.cells'])[0] for t in 
                       [0, 4, 7, 10, 13, 16, 20, 24, 29, 31, 36]}
    times_to_cellsf = {t: list(gcf_new.loc[gcf_new['Library']=='GCF-'+str(t)][lib2 + '.cells'])[0] for t in 
                       [4, 7, 10, 13, 16, 20, 24, 29, 36]}
    # calculting s during exponential phase
    expo_times = [4, 7, 10, 13]
    gens1 = np.log2(np.array([times_to_cells1[t] for t in expo_times]))
    gens2 = np.log2(np.array([times_to_cells2[t] for t in expo_times]))
    gensf = np.log2(np.array([times_to_cellsf[t] for t in expo_times]))
    simple_s_measuring_gc(td, gens1, ['GCP-' + str(t) + '_r1' for t in expo_times], neut_edges, 30, 'pr1_expo_s')
    simple_s_measuring_gc(td, gens2, ['GCP-' + str(t) + '_r2' for t in expo_times], neut_edges, 30, 'pr2_expo_s')
    simple_s_measuring_gc(td, gensf, ['GCF-' + str(t) for t in expo_times], neut_edges, 30, 'f_expo_s')
    # calculating s during the extended saturation phase
    saturation_times = [20, 24, 29]
    simple_s_measuring_gc(td, saturation_times, ['GCP-' + str(t) + '_r1' for t in saturation_times], neut_edges, 30, 'pr1_sat_s')
    simple_s_measuring_gc(td, saturation_times, ['GCP-' + str(t) + '_r2' for t in saturation_times], neut_edges, 30, 'pr2_sat_s')
    simple_s_measuring_gc(td, saturation_times, ['GCF-' + str(t) for t in saturation_times], neut_edges, 30, 'f_sat_s')
    # calculating "lag phase s", which is at least partially confounded by exponential s because it involves the first few divisions.
    simple_s_measuring_gc(td, np.log2(np.array([times_to_cells2[0]/1024, times_to_cellsf[7]])), 
                       ['GCP-0_r2', 'GCF-7'], neut_edges, 30, 'f_lag_s')
    # There are 3 7-hours-after-dilution intervals which we can use to measure lag s in the plates
    # I will not use the 29-36 hr one in my final measure of lag s because it is different from the normal growth conditions (over saturated)
    # but I calculate it anyways, to compare to the other measures later.
    for time_combo in [[0, 7], [24, 31], [29, 36]]:
        simple_s_measuring_gc(td, np.log2(np.array([times_to_cells1[time_combo[0]]/1024, times_to_cells1[time_combo[1]]])), 
                           ['GCP-' + str(t) + '_r1' for t in time_combo], neut_edges, 30, 'pr1_lag_s_' + '_'.join([str(i) for i in time_combo]))
        simple_s_measuring_gc(td, np.log2(np.array([times_to_cells2[time_combo[0]]/1024, times_to_cells2[time_combo[1]]])), 
                           ['GCP-' + str(t) + '_r2' for t in time_combo], neut_edges, 30, 'pr2_lag_s_' + '_'.join([str(i) for i in time_combo]))
        
    td['p_expo_mean'] = np.nanmean(td[['pr1_expo_s', 'pr2_expo_s']], axis=1)
    td['p_sat_mean'] = np.nanmean(td[['pr1_sat_s', 'pr2_sat_s']], axis=1)
    td['pr1_lag_s'] = np.nanmean(td[['pr1_lag_s_0_7', 'pr1_lag_s_24_31']], axis=1)
    td['pr2_lag_s'] = np.nanmean(td[['pr2_lag_s_0_7', 'pr2_lag_s_24_31']], axis=1)
    td['p_lag_mean'] = np.nanmean(td[['pr1_lag_s', 'pr2_lag_s']], axis=1)
    
    td[[i for i in td if i[0] != 's']].to_csv('../../S_Estimation/GC_seg_s_files/GC_' + seg + '_counts_and_s.csv', index=False)
    
# Combining files together
ann_info_file = '../../Mutation_Annotation/Edges_annotated.csv'

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
ann_use = ann[['Edge.ID', 'Edge', 'Gene.Use']]
edge_dat = td.merge(ann_use, on='Edge', how='left')

# ADDING FITNESS INFO
measures = ['p_expo_mean', 'p_lag_mean', 'p_sat_mean', 'f_expo_s', 'f_lag_s', 'f_sat_s']
rep_measures = ['expo_s', 'lag_s', 'sat_s']
names = ['plate exponential s', 'plate lag s', 'plate saturation s', 'flask exponential s', 'flask lag s', 'flask saturation s']
cols = measures + ['pr1_' + c for c in rep_measures] + ['pr2_' + c for c in rep_measures]
for s in segs:
    ed = {i[0]: i[1:] for i in pd.read_csv('../../S_Estimation/GC_seg_s_files/GC_' + s + '_counts_and_s.csv').as_matrix(['Edge'] + cols)}
    for c in range(len(cols)):
        edge_dat[s + '.' + cols[c]] = edge_dat['Edge'].apply(lambda e: ed.setdefault(e, [np.nan for i in range(len(cols))])[c])

edge_dat.to_csv('../../Analysis/GC_data_by_edge.csv', index=False)
