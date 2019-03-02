import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as pl
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
sns.set_style("white")

colors = ['#FFB000', '#648FFF']

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
    'chr16': 948066,
}

segs_w_data_in_both_exps = ['LK4-A04',
 'LK1-E09',
 'LK4-H11',
 'LK2-A12',
 'LK4-A12',
 'LK4-D01',
 'LK3-G02',
 'LK3-D09',
 'LK2-E07',
 'LK2-B04',
 'LK1-C11',
 'LK4-B12',
 'LK1-C09',
 'LK3-D08',
 'LK2-A01',
 'LK1-D06']

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

def get_geno_matrix(seg_names):
    d = pd.read_csv('../accessory_files/BYxRM_GenoData.csv')
    map_genos = {'B': 0, 'R': 1}
    for w in d.keys():
        if change_well_format(w) in seg_names:
            d[change_well_format(w)] = d[w].map(map_genos)
    assert len([s for s in seg_names if s in d.columns]) == len(seg_names)
    return d[['marker'] + seg_names]

def get_dfe(df, segname, show_only='nope'):
    use_df = df.loc[df[segname + '.total.cbcs']>=4]
    if show_only == 'sig_bh':
        pvals = list(use_df[segname + '.pval'])
        pval_sig_boolean = benjamini_hochberg(pvals)[0]
        use_df = use_df.loc[pval_sig_boolean]
    elif show_only == 'sig_uncorrected':
        use_df = use_df.loc[use_df[segname + '.pval']<0.05]
    elif show_only == 'not_near_zero':
        use_df = use_df.loc[(use_df[segname + '.mean.s'] < -0.015) | (use_df[segname + '.mean.s'] > 0.015)]
    return list(use_df[segname + '.mean.s'])

def plot_dfe(sub, df, segname, inset=None, n_bins=20, n_inset_bins=20, single_plot_output_name=None):
    if single_plot_output_name:
        f, sub = pl.subplots(1, 1, figsize=(5, 4))
    sub.hist(np.clip(get_dfe(df, segname), -0.2, 0.1), bins=n_bins, facecolor='#333333')
    sub.set_xlim([-0.2, 0.1])
    if inset:
        ax = inset_axes(sub, height="40%", width="40%", loc=2, bbox_to_anchor=(0.1, 0, 1, 0.9), bbox_transform=sub.transAxes)
        ax.hist(np.clip(get_dfe(df, segname, show_only=inset), -0.2, 0.1), bins=n_inset_bins, facecolor='#333333')
        ax.set_xlim([-0.2, 0.1])
        ax.set_xticks([])
        ax.tick_params(axis='y', which='major', labelsize=14)
    if single_plot_output_name:
        sns.despine()
        sub.tick_params(axis='both', which='major', labelsize=14)
        sub.set_xlabel('Fitness Effect', fontsize=16)
        pl.tight_layout()        
        f.savefig(single_plot_output_name, background='transparent')
        pl.close('all')

def make_jointplots(df, cols, output_file, col_names=['match_em'], limit_size=True, sig_column=None):
    if col_names[0] == 'match_em':
        col_names = cols
    if len(cols) == 2:
        f = pl.figure(figsize=(6, 5))
        gs0 = gridspec.GridSpec(5, 5)
        scatters = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1:,:4])
        y_hists = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[1:,4:])
        x_hists = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[:1,:4])
        cols_x = [cols[0]]
        cols_y = [cols[1]]
        cn_x = [col_names[0]]
        cn_y = [col_names[1]]
    else:
        f = pl.figure(figsize=(4*len(cols) + 2, 3*len(cols) + 2))
        gs0 = gridspec.GridSpec(4*len(cols) + 2, 3*len(cols) + 2)
        scatters = gridspec.GridSpecFromSubplotSpec(len(cols), len(cols), subplot_spec=gs0[2:,:3*len(cols)])
        y_hists = gridspec.GridSpecFromSubplotSpec(len(cols), 1, subplot_spec=gs0[2:,3*len(cols):])
        x_hists = gridspec.GridSpecFromSubplotSpec(1, len(cols), subplot_spec=gs0[:2,:3*len(cols)])
        cols_x, cols_y, cn_x, cn_y = cols, cols, col_names, col_names

    for i in range(len(cols_x)):
        for j in range(len(cols_y)):
            sub_s = pl.Subplot(f, scatters[j, i])
            f.add_subplot(sub_s)
            td = df.loc[pd.notnull(df[cols_x[i]])].loc[pd.notnull(df[cols_y[j]])]
            sub_s.axhline(y=0, xmin=0, xmax=1, color='#333333', linestyle='dashed', alpha=0.5)
            sub_s.axvline(x=0, ymin=0, ymax=1, color='#333333', linestyle='dashed', alpha=0.5)
            if sig_column:
                td1 = td.loc[td[sig_column]<0.05]
                td2 = td.loc[~td['Edge'].isin(td1['Edge'])] 
                sub_s.scatter(td1[cols_x[i]], td1[cols_y[j]], color='#333333')
                sub_s.scatter(td2[cols_x[i]], td2[cols_y[j]], color='#333333', alpha=0.4)
            else:
                sub_s.scatter(td[cols_x[i]], td[cols_y[j]], color='#333333')
            if j == (len(cols_y)-1): 
                sub_s.set_xlabel(cn_x[i], fontsize=20)
                sub_xh = pl.Subplot(f, x_hists[i])
                f.add_subplot(sub_xh)
                sub_xh.hist(td[cols_x[i]], bins=20, color='#333333')
                sub_xh.set_xticks([])
                sub_xh.set_yticks([])
            if i == 0: 
                sub_s.set_ylabel(cn_y[j], fontsize=20)
                sub_yh = pl.Subplot(f, y_hists[j])
                f.add_subplot(sub_yh)
                sub_yh.hist(td[cols_y[j]], bins=20, orientation="horizontal", color='#333333')
                sub_yh.set_xticks([])
                sub_yh.set_yticks([])
            
            sub_s.tick_params(axis='both', which='major', labelsize=16)
            
            sns.despine(bottom=True, left=True)

    #pl.gcf().subplots_adjust(bottom=0.2)
    pl.tight_layout()
    f.savefig(output_file, background='transparent')
    pl.close('all')
    
def plot_one(sub, df_row, segs, x_list, gm, qtl_color=False, show_pvals=False):
    y_list = list(df_row[segs])
    pval = df_row['x_model_p']
    if show_pvals:
        if pval < 0.001:
            sub.annotate('p=' + "{:.2E}".format(pval), xy=(0.9, 0.9), xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')
        else:
            sub.annotate('p=' + str(pval)[:4], xy=(0.9, 0.9), xycoords='axes fraction', fontsize=16, horizontalalignment='right', verticalalignment='bottom')

    sub.tick_params(axis='both', which='major', labelsize=12)
    if not qtl_color:
        sub.scatter(x_list, y_list, c='#333333')
    else:
        geno_dat = gm.as_matrix(['marker'] + segs)
        qtl_str = str(df_row['full.model.qtls'])
        if qtl_str != 'nan':
            edge_qtls = [i.split(';')[0] for i in qtl_str.split('|')]
        else:
            edge_qtls = []
        if len(edge_qtls) == 1:
            for eq in edge_qtls:
                geno_row = [row[1:] for row in geno_dat if '_'.join(row[0].split('_')[1:3]) == eq][0]
                for i in range(2):
                    xs = [x_list[s] for s in range(len(segs)) if geno_row[s]==i]
                    ys = [y_list[s] for s in range(len(segs)) if geno_row[s]==i]
                    sub.scatter(x=xs, y=ys, marker='o', c=colors[i], s=25)
        elif len(edge_qtls) > 1: # only distinguish maximum 2 QTLs
            geno_row1 = [row[1:] for row in geno_dat if '_'.join(row[0].split('_')[1:3]) == edge_qtls[0]][0]
            geno_row2 = [row[1:] for row in geno_dat if '_'.join(row[0].split('_')[1:3]) == edge_qtls[1]][0]
            for i in range(2):
                xs1 = [x_list[s] for s in range(len(segs)) if geno_row1[s]==i and geno_row2[s]==0]
                ys1 = [y_list[s] for s in range(len(segs)) if geno_row1[s]==i and geno_row2[s]==0]
                xs2 = [x_list[s] for s in range(len(segs)) if geno_row1[s]==i and geno_row2[s]==1]
                ys2 = [y_list[s] for s in range(len(segs)) if geno_row1[s]==i and geno_row2[s]==1]
                sub.scatter(x=xs1, y=ys1, marker='o', c=colors[i], s=25)
                sub.scatter(x=xs2, y=ys2, marker='v', linewidth=1, facecolors='none', edgecolors=colors[i], s=25)
        else:
            xs = [x_list[s] for s in range(len(segs))]
            ys = [y_list[s] for s in range(len(segs))]
            sub.scatter(x=xs, y=ys, marker='o', c='k', s=25)
    sub.set_xlim([-0.16, 0.12])
    y_range = max(y_list)-min(y_list)
    sub.set_ylim([min(y_list)-y_range*0.25, max(y_list)+y_range*0.25]) # I don't know why pyplot is not doing this well automatically here
    
def make_dfe_plot(exp, outname):
    hm_x, hm_y = data_by_fit_ranks[exp]
    dfe_range = [-0.15, 0.07]
    f = pl.figure(figsize=(26, 10))
    gs0 = gridspec.GridSpec(30, 48)
    top_fit_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[:6,1:12])
    fit_heatmaps_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[7:26,1:12])
    stats_gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[:26,30:38])
    stats_gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[:26,40:])
    dfe_combined_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs0[:26,17:25])
    top_sub = pl.Subplot(f, top_fit_gs[0])
    hm_sub = pl.Subplot(f, fit_heatmaps_gs[0])
    stats_subs = [[pl.Subplot(f, stats_gs1[j]) for j in range(2)], [pl.Subplot(f, stats_gs2[j]) for j in range(2)]]
    dfe_combined_subs = [pl.Subplot(f, dfe_combined_gs[j]) for j in range(2)]
    if len(segs_use[exp]) > 40:
        combine_bin = len(segs_use[exp])/3
    else:
        combine_bin = len(segs_use[exp])
    hm_sub.hist2d(hm_x, hm_y, bins=[combine_bin, 30], cmin=1, cmap=pl.cm.viridis)
    hm_sub.tick_params(axis='both', which='major', labelsize=16)

    hm_sub.set_ylim(dfe_range)
    hm_sub.annotate('Background Fitness Rank', fontsize=20, xy=(0.5, -0.15), xycoords="axes fraction", horizontalalignment="center")
    hm_sub.set_ylabel('Fitness Effect', fontsize=20)
    #hm_sub.set_yticks([])
    f.add_subplot(hm_sub)

    top_sub.plot([i for i in range(len(sorted_segs[exp]))], [seg_to_fit[s] for s in sorted_segs[exp]])
    top_sub.tick_params(axis='both', which='major', labelsize=16)
    top_sub.set_xticks([])
    top_sub.set_ylim([-0.12, 0.12])
    top_sub.set_ylabel('Background\nFitness', fontsize=20)
    f.add_subplot(top_sub)

    num_to_combine = int(len(segs_use[exp])/4) #plotting quartiles
    stats = ['mean', 'skew', 'significant.beneficial.mutations', 'significant.deleterious.mutations']
    stat_names = {'mean': 'DFE Mean', 'skew': 'DFE skew', 'significant.beneficial.mutations': '# of Bene. Mutations', 
                  'significant.deleterious.mutations': '# of Del. Mutations'}
    i = 0
    for j in range(2):
        for k in range(2):
            if exp == 'TP':
                plot_one(stats_subs[j][k], dats['TP.DFE'].loc[dats['TP.DFE']['DFE.statistic']==stats[i]].iloc[0], segs_use['TP'], 
                         [seg_to_fit[s] for s in segs_use['TP']], gm, qtl_color=True)
            else:
                plot_one(stats_subs[j][k], dats['BT.DFE'].loc[dats['BT.DFE']['DFE.statistic']==stats[i]].iloc[0], segs_use['BT'], 
                         [seg_to_fit[s] for s in segs_use['BT']], gm, show_pvals=True)
            stats_subs[j][k].set_ylabel(stat_names[stats[i]], fontsize=20)
            if k == 0:
                stats_subs[j][k].set_xticks([])
            f.add_subplot(stats_subs[j][k])
            i += 1
    stats_subs[0][1].annotate('Background Fitness', fontsize=20, xy=(1.2, -0.24), xycoords="axes fraction", horizontalalignment="center")

    top_dfe, bottom_dfe = [], [],    
    for i in range(num_to_combine):
        bottom_dfe += get_dfe(dats[exp], sorted_segs[exp][i])
        top_dfe += get_dfe(dats[exp], sorted_segs[exp][len(sorted_segs[exp])-1-i])
        
    dfe_combined_subs[0].hist(top_dfe, bins=50, color="#333333")
    dfe_combined_subs[1].hist(bottom_dfe, bins=50, color="#333333")

    dfe_combined_subs[1].annotate('Fitness Effect', fontsize=20, xy=(0.5, -0.24), xycoords="axes fraction", horizontalalignment="center")

    top_sub.annotate('A', fontsize=40, xy=(-0.2, 1.2), xycoords="axes fraction", horizontalalignment="center")
    stats_subs[0][0].annotate('C', fontsize=40, xy=(-0.2, 1.1), xycoords="axes fraction", horizontalalignment="center")
    dfe_combined_subs[0].annotate('B', fontsize=40, xy=(-0.2, 1.1), xycoords="axes fraction", horizontalalignment="center")

    dfe_combined_subs[0].annotate('Most Fit Quartile\nof Segregants', fontsize=14, xy=(0.32, 0.7), xycoords="axes fraction", horizontalalignment="center")
    dfe_combined_subs[1].annotate('Least Fit Quartile\nof Segregants', fontsize=14, xy=(0.32, 0.7), xycoords="axes fraction", horizontalalignment="center")

    jnk = [f.add_subplot(dfe_combined_subs[j]) for j in range(2)]
    jnk = [dfe_combined_subs[x].set_xlim(dfe_range) for x in range(2)]
    jnk = [dfe_combined_subs[x].tick_params(axis='both', which='major', labelsize=16) for x in range(2)]
    jnk = [stats_subs[y][x].tick_params(axis='both', which='major', labelsize=16) for x in range(2) for y in range(2)]

    sns.despine()
    sns.despine(ax=hm_sub, bottom=True)
    sns.despine(ax=top_sub, bottom=True, left=True)
    f.savefig(outname, background='transparent')
    pl.close('all')

def make_correlation_plot(segs, td, xvar, yvar, xerr_var, yerr_var, xlabel, ylabel, criteria_one, criteria_two, output_name):
    nrows = int(np.ceil(len(segs)/4))
    fig, subps = pl.subplots(nrows, 4, figsize=(16, nrows*(16/5)), sharex=True, sharey=True)
    if nrows == 1:
        subps = [subps]
    for i in range(nrows):
        for j in range(4):
            if i*4 + j < len(segs):
                seg = segs[i*4 + j]
                d = td.loc[td[seg+criteria_one[0]] >= criteria_one[1]].loc[td[seg+criteria_two[0]] >= criteria_two[1]]
                subps[i][j].plot([-0.15, 0.1], [-0.15, 0.1], linestyle='dashed', c='k', alpha=0.5)
                subps[i][j].errorbar(data=d, x=seg+xvar, y=seg+yvar, xerr=seg+xerr_var, yerr=seg+yerr_var, c='k', 
                                     marker='.', linestyle='', alpha=0.5)
                subps[i][j].set_xlim([-0.15, 0.1])
                subps[i][j].set_ylim([-0.15, 0.1])
                subps[i][j].tick_params(axis='both', which='major', labelsize=12)
                if j == 0:
                    subps[i][j].set_ylabel(ylabel, fontsize=14)
                if i == nrows-1:
                    subps[i][j].set_xlabel(xlabel, fontsize=14)
                subps[i][j].set_title(seg, y=0.9, fontsize=16)
            else:
                subps[i][j].axis('off')
    sns.despine()
    pl.tight_layout()
    fig.savefig(output_name, background='transparent')
    pl.close('all')

def make_single_determinant_plot(sub, df_row, segs, gm, plot_errors, show_title, single_plot_output_name=None):
    if single_plot_output_name:
        f, sub = pl.subplots(1, 1, figsize=(5, 4))
    measured = [seg for seg in segs if df_row[seg + '.total.cbcs'] >= 4]
    geno_dat = gm.as_matrix(['marker'] + measured)
    sub.axhline(y=0, xmin=0, xmax=1, color='#333333', linestyle='dashed', alpha=0.5)
    if show_title:
        sub.set_title(str(df_row['Gene.Use']), fontsize=14)
    sub.tick_params(axis='both', which='major', labelsize=12)
    if str(df_row['full.model.qtls']) != 'nan':
        edge_qtls_unsorted = [i.split(';')[0] for i in str(df_row['full.model.qtls']).split('|')]
        qtl_effects = [float(i) for i in str(df_row['full_model_coeffs']).split(';')[2:]]
        edge_qtls = [i for jnk,i in sorted(zip(qtl_effects, edge_qtls_unsorted), reverse=True)]
    else:
        edge_qtls = []
    if len(edge_qtls) == 1:
        for eq in edge_qtls:
            geno_row = [row[1:] for row in geno_dat if '_'.join(row[0].split('_')[1:3]) == eq][0]
            for i in range(2):
                xs = [seg_to_fit[measured[s]] for s in range(len(measured)) if geno_row[s]==i]
                ys = [df_row[measured[s] + '.mean.s'] for s in range(len(measured)) if geno_row[s]==i]
                ye =[df_row[measured[s] + '.stderr.s'] for s in range(len(measured)) if geno_row[s]==i]
                sub.scatter(x=xs, y=ys, marker='o', c=colors[i], s=25)
                if plot_errors:
                    sub.errorbar(x=xs, y=ys, yerr=ye, marker='', c=colors[i], linestyle='')
    elif len(edge_qtls) > 1:
        geno_row1 = [row[1:] for row in geno_dat if '_'.join(row[0].split('_')[1:3]) == edge_qtls[0]][0]
        geno_row2 = [row[1:] for row in geno_dat if '_'.join(row[0].split('_')[1:3]) == edge_qtls[1]][0]
        for i in range(2):
            xs1 = [seg_to_fit[measured[s]] for s in range(len(measured)) if geno_row1[s]==i and geno_row2[s]==0]
            ys1 = [df_row[measured[s] + '.mean.s'] for s in range(len(measured)) if geno_row1[s]==i and geno_row2[s]==0]
            ye1 = [df_row[measured[s] + '.stderr.s'] for s in range(len(measured)) if geno_row1[s]==i and geno_row2[s]==0]
            xs2 = [seg_to_fit[measured[s]] for s in range(len(measured)) if geno_row1[s]==i and geno_row2[s]==1]
            ys2 = [df_row[measured[s] + '.mean.s'] for s in range(len(measured)) if geno_row1[s]==i and geno_row2[s]==1]
            ye2 = [df_row[measured[s] + '.stderr.s'] for s in range(len(measured)) if geno_row1[s]==i and geno_row2[s]==1]
            sub.scatter(x=xs1, y=ys1, marker='o', c=colors[i], s=25)
            if plot_errors:
                sub.errorbar(x=xs1, y=ys1, yerr=ye1, marker='', c=colors[i], linestyle='')
            if plot_errors:
                sub.errorbar(x=xs2, y=ys2, yerr=ye2, marker='', linestyle='', c=colors[i], alpha=0.5)
            sub.scatter(x=xs2, y=ys2, marker='o', linewidth=1, facecolors='none', edgecolors=colors[i], s=25, alpha=0.5)
    else:
        xs = [seg_to_fit[measured[s]] for s in range(len(measured))]
        ys = [df_row[measured[s] + '.mean.s'] for s in range(len(measured))]
        ye =[df_row[measured[s] + '.stderr.s'] for s in range(len(measured))]
        sub.scatter(x=xs, y=ys, marker='o', c='k', s=25)
        if plot_errors:
            sub.errorbar(x=xs, y=ys, yerr=ye, marker='', c='k', linestyle='') 
    sub.set_xlim([-0.16, 0.12])
    sub.set_ylim([-0.2, 0.08])
    if single_plot_output_name:
        sns.despine()
        sub.tick_params(axis='both', which='major', labelsize=14)
        sub.set_xlabel('Background Fitness', fontsize=16)
        sub.set_ylabel('Fitness Effect', fontsize=16)
        pl.tight_layout()
        f.savefig(single_plot_output_name, background='transparent')
        pl.close('all')
    
def plot_20_determinants(df_rows, segs, output_name, plot_errors=False, show_title=True, big_title=False):
    nrows = int(np.ceil(len(df_rows)/4))
    fig, subps = pl.subplots(nrows, 4, figsize=(16, nrows*(16/5)), sharex=True, sharey=True)
    if nrows == 1:
        subps = [subps]
    for i in range(nrows):
        for j in range(4):
            if i*4 + j < len(df_rows):
                make_single_determinant_plot(subps[i][j], df_rows[i*4 + j], segs, gm, plot_errors, show_title)
                subps[i][j].tick_params(axis='both', which='major', labelsize=12)
                if j == 0:
                    subps[i][j].set_ylabel('Fitness Effect', fontsize=14)
                if i == nrows-1:
                    subps[i][j].set_xlabel('Background Fitness', fontsize=14)
            else:
                subps[i][j].axis('off')
    if big_title:
        fig.suptitle(big_title)
    sns.despine()
    pl.tight_layout()
    if output_name:
        fig.savefig(output_name, background='transparent')
        pl.close('all')
        
        
def make_determinants_figure(outname, gm, plot_errors=False, plot_std_dev=True):
    
    df = tp

    f = pl.figure(figsize=(24, 28))
    pl.subplots_adjust(wspace=0.4)
    gs0 = gridspec.GridSpec(10, 52)

    top_sub = pl.subplot(gs0[:2,2:50])
    gs01 = gridspec.GridSpecFromSubplotSpec(5, 3, subplot_spec=gs0[3:,2:45])
    subps = [[pl.Subplot(f, gs01[i, j]) for j in range(3)] for i in range(5)]
    
    example_loc = len(df)+7
    xbars = np.array([i for i in range(len(df))] + [example_loc])
    gene_descrips = list(df['Gene.Use'])
    edges_in_order = list(df['Edge'])
    pl.xticks(xbars, [i.split(' ')[1] for i in gene_descrips], rotation='vertical', fontsize=12)
    ticks = top_sub.get_xticklabels()
    cd = {'in': 'red', 'nearby': 'black'}
    jnk = [ticks[i].set_color(cd[gene_descrips[i].split(' ')[0]]) for i in range(len(gene_descrips))]
    for i, label in enumerate(ticks):
        label.set_y(label.get_position()[1] + 0.025)
    top_sub.tick_params(axis='y', which='major', labelsize=14)
    top_sub.set_xlim([-1, len(gene_descrips)+16])
    top_sub.set_ylim([0, 0.055])
    ex_dev, ex_full, ex_qtl, ex_x = 0.05, 0.03, 0.022, 0.007
    top_sub.bar(xbars, list(np.sqrt(df['var'])) + [ex_dev], color='none', edgecolor='black', width=0.7)
    top_sub.bar(xbars, list(np.sqrt(df['full_sig_only_r2']*df['var'])) + [ex_full], color='#BBBBBB', width=0.7)
    top_sub.scatter(xbars-0.05, list(np.sqrt(df['x_sig_only_r2']*df['var'])) + [ex_x], color='#222222', marker='x', zorder=3)
    top_sub.scatter(xbars-0.05, list(np.sqrt(df['qtl_sig_only_r2']*df['var'])) + [ex_qtl], color=colors[1], marker='o', zorder=4)

    top_sub.annotate('$\sigma_{XM}$\n(variance\nexplained by\nX model)$^{1/2}$', xy=(example_loc+8, ex_x), xycoords='data', fontsize=20, ha='center', va='center')
    top_sub.plot([example_loc+0.7, example_loc+3], [ex_x, ex_x], c='k')

    top_sub.annotate("$\sigma_{FM}$\n(variance\nexplained by\nfull model)$^{1/2}$", xy=(example_loc+8, ex_full), xycoords='data', fontsize=20, ha='center', va='center')
    top_sub.plot([example_loc+0.7, example_loc+3], [ex_full, ex_full], c='k')

    top_sub.annotate('$\sigma_{QM}$\n(variance\nexplained by\nQTL model)$^{1/2}$', xy=(example_loc-8, ex_qtl), xycoords='data', fontsize=20, ha='center', va='center')
    top_sub.plot([example_loc-0.7, example_loc-3], [ex_qtl, ex_qtl], c='k')

    top_sub.annotate('$\sigma_s$\n(std dev of\nfitness effect)', xy=(example_loc-8, ex_dev), xycoords='data', fontsize=20, ha='center', va='center')
    top_sub.plot([example_loc-0.7, example_loc-3], [ex_dev, ex_dev], c='k')
    
    examples = [
        'GTTTAGCTTCCGTTG',  #in RPL16A
        'TCAAAGCATGAAAAA',  #in SUM1
        'CTTTCTTGTGTATTT',  #in BRR1
        'CCAACACAGGCTTCG',  #in NOP16
        'TGGAGTCTTTGTTGA',  #in NOT3
        'TTATATTTATTTGCT',  #in SIR3
        'CTACTTACAACGGAA',  #in KAP123
        'ATTATCATCGCCATC',  #in PIH1
        'AGTGTTAATCAGACC',  #nearby PAH1
        'GAACTCAGGTTCCAT',  #in MME1
        'AGTGTATGATAATAT',  #nearby KRI1
        'GTACAAGAAATTTTG',  #nearby PDE2
        'TATATTGAACTTTAC', 'TCAAAACGGAGTGTT', 'ACAACCTACCTGCTA' # references
    ]

    c = 0
    for subarr in subps:
        for sub in subarr:
            edge = examples[c]
            c += 1
            if c < 13:
                top_sub.annotate(str(c), ((edges_in_order.index(edge)+0.6)/(len(gene_descrips)+17), -0.22), fontsize=10, xycoords='axes fraction', bbox={"boxstyle" : "circle", "color":"#BBBBBB"}, zorder=5)
                sub.annotate(str(c),(-0.14,0.06), fontsize=12, xycoords='data', bbox={"boxstyle" : "circle", "color":"#BBBBBB"})
            df_row = tp_all.loc[tp_all['Edge']==edge].iloc[0]
            make_single_determinant_plot(sub, df_row, segs, gm, True, True)
            sub.set_xlim([-0.16, 0.12])
            sub.set_ylim([-0.2, 0.09])
            f.add_subplot(sub)
    
    subps[1][0].annotate('Fitness Effect', fontsize=25, xy=(-0.25, -0.16), xycoords='axes fraction', rotation=90, ha='center', va='center')
    subps[3][1].annotate('Background Fitness', fontsize=25, xy=(0.5, -0.2), xycoords='axes fraction', ha='center', va='top')
    sns.despine()
    sns.despine(left=True, bottom=True, ax=top_sub)
            
    f.savefig(outname, background='transparent')



## READING DATA
x_info = pd.read_csv('../accessory_files/Clones_For_Tn96_Experiment.csv')
seg_to_fit = {i[0]: i[1] for i in x_info.as_matrix(['segregant', 'initial fitness, YPD 30C'])}
tp_all = pd.read_csv('../../Analysis/TP_data_by_edge.csv')
tp = tp_all.loc[tp_all['Type']=='Experiment']
bt = pd.read_csv('../../Analysis/BT_data_by_edge.csv')
tp_dfe = pd.read_csv('../../Analysis/TP_DFE_statistics.csv')
bt_dfe = pd.read_csv('../../Analysis/BT_DFE_statistics.csv')
bt['Long.Edge'] = bt['Edge']
bt['Edge'] = bt['Long.Edge'].str[:15]
dats = {'BT': bt, 'TP': tp, 'BT.DFE': bt_dfe, 'TP.DFE': tp_dfe}
exps = {'BT': 'MM', 'TP': 'FM'}
segs_all = {exp: [i.split('.')[0] for i in dats[exp] if '.mean.s' in i] for exp in exps}
segs_use = {exp: [s for s in segs_all[exp] if len(dats[exp].loc[dats[exp][s + '.total.cbcs']>=4])>50] for exp in exps}
gm = get_geno_matrix(segs_all['TP'])

# Make main determinant figure
make_determinants_figure('../../Figures/Genetic_Determinants.png', gm, plot_errors=True)

## Making correlation plots
make_correlation_plot(segs_w_data_in_both_exps, bt, '.rep1.s', '.rep2.s', '.rep1.stderr.s', '.rep2.stderr.s', 'mean bc s rep 1', 'mean bc s rep 2',
                     ('.rep1.cbcs', 2), ('.rep2.cbcs', 2), '../../Figures/all/correlations/' + exps['BT'] + '_replicate_correlations.png')
make_correlation_plot(segs_w_data_in_both_exps, tp, '.rep1.s', '.rep2.s', '.rep1.stderr.s', '.rep2.stderr.s', 'mean bc s rep 1', 'mean bc s rep 2',
                     ('.rep1.cbcs', 2), ('.rep2.cbcs', 2), '../../Figures/all/correlations/' + exps['TP'] + '_replicate_correlation_examples.png')
dm = tp.merge(bt, on='Edge', how='inner', suffixes=('_TP', '_BT'))
make_correlation_plot(segs_w_data_in_both_exps, dm, '.mean.s_BT', '.mean.s_TP', '.stderr.s_BT', '.stderr.s_TP', exps['BT'] + ' s', exps['TP'] + ' s',
                     ('.total.cbcs_BT', 4), ('.total.cbcs_TP', 4), '../../Figures/all/correlations/' + exps['TP'] + '_' +  exps['BT'] + '_correlations.png')
segs_tp = segs_use['TP']
for i in range(int(np.ceil(len(segs_tp)/20))):
    make_correlation_plot(segs_tp[i*20:(i+1)*20], tp, '.rep1.s', '.rep2.s', '.rep1.stderr.s', '.rep2.stderr.s', 'mean bc s rep 1', 'mean bc s rep 2',
                         ('.rep1.cbcs', 2), ('.rep2.cbcs', 2), '../../Figures/all/correlations/' + exps['TP'] + '_replicate_correlations_' + str(i+1) + '.png')

# Making determinant plots for TP
for it in tp.iterrows():
    row = it[1]
    make_single_determinant_plot(None, row, segs_tp, gm, False, True, single_plot_output_name='../../Figures/all/TP_determinants/' + str(row['Gene.Use']).replace(' ', '_') + '_' + row['Edge'] + 'determinants.png')
    make_single_determinant_plot(None, row, segs_tp, gm, True, True, single_plot_output_name='../../Figures/all/TP_determinants/' + str(row['Gene.Use']).replace(' ', '_') + '_' + row['Edge'] + 'determinants_w_error.png')

model_groups = {
    'no_model': tp.loc[~(tp['segregant_model_p']<0.05)].loc[~((tp['qtl_model_p']<0.05) | (tp['x_model_p']<0.05))],
    'seg_model_only': tp.loc[tp['segregant_model_p']<0.05].loc[~((tp['qtl_model_p']<0.05) | (tp['x_model_p']<0.05))],
    'x_model_only': tp.loc[tp['model_comp_p_full_vs_qtl']<0.05].loc[~(tp['model_comp_p_full_vs_x']<0.05)],
    'qtl_model_only': tp.loc[tp['model_comp_p_full_vs_x']<0.05].loc[~(tp['model_comp_p_full_vs_qtl']<0.05)],
    'both_in_model': tp.loc[tp['model_comp_p_full_vs_x']<0.05].loc[tp['model_comp_p_full_vs_qtl']<0.05]
}
for m in model_groups:
    rows = [r[1] for r in model_groups[m].sort_values(by='Gene.Use').iterrows()]
    for j in range(int(np.ceil(len(rows)/20))):
        plot_20_determinants(rows[j*20:(j+1)*20], segs, 'Figures/all/TP_determinants/TP_determinants_' + m + '_' + str(j+1) + '.png')
        plot_20_determinants(rows[j*20:(j+1)*20], segs, 'Figures/all/TP_determinants/TP_determinants_' + m + '_' + str(j+1) + '_error.png', True)

## Making DFE plots
data_by_fit_ranks = {'BT': [[], []], 'TP': [[], []]}
sorted_segs = {exp: sorted(segs_use[exp], key=lambda x: seg_to_fit[x]) for exp in exps}
fit_ranks = {exp: {s: sorted_segs[exp].index(s) for s in segs_use[exp]} for exp in exps}
for exp in exps:
    for seg in segs_use[exp]:
        td = dats[exp].loc[dats[exp][seg + '.total.cbcs']>=4]
        data_by_fit_ranks[exp][0] += [fit_ranks[exp][seg]] * len(td)
        data_by_fit_ranks[exp][1] += list(np.clip(td[seg + '.mean.s'], -0.15, 0.2))
        # plotting individual DFEs
        if exp == 'TP':
            plot_dfe(None, dats[exp], seg, single_plot_output_name='../../Figures/all/DFEs/TP/TP_DFE_' + seg + '.png')
        else:
            plot_dfe(None, dats[exp], seg, n_bins=25, n_inset_bins=25, inset='sig_uncorrected', single_plot_output_name='../../Figures/all/DFEs/BT/BT_DFE_' + seg + '.png')
    #dats[exp]['x_p_any_type'] = dats[exp].apply(lambda row: np.nanmax([row['x_model_p'], row['model_comp_p_full_vs_qtl']]), axis=1)
    #dats[exp]['top_qtl_es'] = dats[exp]['full_model_qtl_effect_sizes'].apply(lambda es: float(str(es).split(';')[0]))
exp = 'TP'
dats[exp]['x_p_any_type'] = dats[exp].apply(lambda row: np.nanmax([row['x_model_p'], row['model_comp_p_full_vs_qtl']]), axis=1)
dats[exp]['top_qtl_es'] = dats[exp]['full_model_qtl_effect_sizes'].apply(lambda es: float(str(es).split(';')[0]))

make_dfe_plot('TP', '../../Figures/' + exps['TP'] + '_dfe_fig.png')
make_dfe_plot('BT', '../../Figures/all/' + exps['BT'] + '_dfe_fig.png')

## Making effect size / genetic determinants plots
make_jointplots(dats['TP'], ['avg_s', 'full_model_x_slope'], '../../Figures/all/X_slopes_fig.png', col_names=['Average Fitness Effect', 'Regression slope'], sig_column='x_p_any_type')

make_jointplots(dats['TP'], ['avg_s', 'full_model_x_effect_size_measure'], '../../Figures/all/X_effect_size_fig.png', col_names=['Average Fitness Effect', 'Background Fitness\nEffect Size'], sig_column='x_p_any_type')

make_jointplots(dats['TP'], ['avg_s', 'top_qtl_es'], '../../Figures/all/QTL_effect_size_fig.png', col_names=['Average Fitness Effect', 'Top QTL\nEffect Size'], sig_column='model_comp_p_full_vs_x')

make_jointplots(dats['TP'], ['full_model_x_effect_size_measure', 'top_qtl_es'], '../../Figures/all/QTL_vs_X_effect_size_fig.png', col_names=['Background Fitness\nEffect Size', 'Top QTL\nEffect Size'], sig_column='model_comp_p_full_vs_x')


## Various other plots

# QTL map
fig, subs = pl.subplots(16, 1, figsize=(16, 16), sharex=True, sharey=True)
for entry in tp['full.model.qtls']:
    if str(entry) != 'nan':
        for qtl in str(entry).split('|'):
            hit, left, right = qtl.split(';')
            chromo = hit[3:5]
            hit_loc = int(hit.split('_')[1])
            if left[3:5] != chromo:
                left_loc = 0
            else:
                left_loc = int(left.split('_')[1])
            if right[3:5] != chromo:
                right_loc = chromo_lens[chromo]
            else:
                right_loc = int(right.split('_')[1])
            random_y = np.random.random()
            subs[int(chromo)-1].plot([left_loc, right_loc], [random_y, random_y], c='k')
            subs[int(chromo)-1].scatter(hit_loc, random_y, c='r')
            
for i in range(16):
    subs[i].plot([0, 0], [0, 0.3], c='b')
    subs[i].plot([chromo_lens['chr' + str(i+1).zfill(2)], chromo_lens['chr' + str(i+1).zfill(2)]], [0, 0.3], c='b')
    subs[i].set_yticks([])
    subs[i].annotate(str(i+1), (0, 0.4), fontsize=20)
sns.despine(left=True)    
fig.savefig('../../Figures/all/QTL_map.png', background='transparent')

# QTL effect sizes
qtl_df = pd.read_csv('../../Analysis/QTL_results.csv')
mhq_df = pd.read_csv('../../Analysis/Multi_hit_QTLs.csv')
f, sub = pl.subplots(1, 1, figsize=(7, 5))
qtl_group_2_simpler = {i[0]: i[0].split('_')[1] + '_' + str(int(i[1])) for i in mhq_df.as_matrix(['QTL', 'median.loc'])}
qtl_df['QTL_group_simple'] = qtl_df['QTL_group'].apply(lambda q: qtl_group_2_simpler.setdefault(q, 'Other'))
sns.swarmplot(x='QTL_group_simple', y='effect_size', data=qtl_df.sort_values(by='QTL_group_simple'), ax=sub)
pl.xticks(rotation='vertical')
sub.tick_params(axis='both', which='major', labelsize=16)
sns.despine()
sub.set_xlabel('Multi-hit QTL', fontsize=16)
sub.set_ylabel('Epistatic effect', fontsize=16)
pl.tight_layout()
f.savefig('../../Figures/all/QTL_group_effect_sizes.png', background='transparent')