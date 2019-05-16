from matplotlib import pyplot as pl
import seaborn as sns
from glob import glob
import numpy as np
from matplotlib.collections import LineCollection
matplotlib.use('Agg')

exps = ['TP', 'BT']
reps = ['rep1', 'rep2']
for exp in exps:
    segs = [i.split('/')[-1] for i in glob('../../S_Estimation/' + exp + '_output/*')]
    fig, subps = pl.subplots(4, 4, sharex=True, sharey=True, figsize=(16, 16))
    for x in range(4):
        for y in range(4):
            found = False
            while (not found):
                seg = np.random.choice(segs)
                rep = np.random.choice(reps)
                fpath = '../../S_Estimation/' + exp + '_output/' + seg + '/' + seg + '_' + rep + '_bc_s_initial.csv'
                try:
                    d = pd.read_csv(fpath)
                    edges_w_outliers = set(d.loc[d['LL.ratio']>15]['Edge'])
                    if len(edges_w_outliers) > 0:
                        found = True
                        tp_names = [i for i in d.columns if '-T' in i]
                        tps = [int(i.split('-T')[-1])*10 for i in tp_names]
                        for tp in tp_names:
                            s = np.sum(d[tp])
                            d[tp + '.log10freq'] = np.log10(np.clip(d[tp] / s, 10**-6, 1))
                        rows = d.as_matrix([tp + '.log10freq' for tp in tp_names])
                        all_lines = np.zeros((len(rows), len(tps), 2), float)
                        for j in range(len(rows)):
                            all_lines[j, :, 1] = rows[j]
                            all_lines[j, :, 0] = tps
                        lines = LineCollection(all_lines, color='k', alpha=0.1, linewidths=1)
                        subps[x][y].add_collection(lines)
                        edge_to_use = np.random.choice(list(edges_w_outliers))
                        td = d.loc[d['Edge']==edge_to_use]
                        for entry in td.loc[td['LL.ratio'] <= 15].as_matrix([tp + '.log10freq' for tp in tp_names]):
                            subps[x][y].plot(tps, entry, c='b')
                        for entry in td.loc[td['LL.ratio'] > 15].as_matrix([tp + '.log10freq' for tp in tp_names]):
                            subps[x][y].plot(tps, entry, c='r')
                        subps[x][y].set_ylim([-6, 0])
                        subps[x][y].set_xlim([0, 40])
                        subps[x][y].set_title(seg + ' ' + rep)
                except FileNotFoundError:
                    print(fpath, 'not found')
    sns.despine()
    pl.tight_layout()
    fig.savefig('../../Figures/supp_figs/' + exp + '_outlier_examples.png', background='transparent', bbox_inches='tight', pad_inches=0.1)
    pl.close("all")
