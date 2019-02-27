import pandas as pd
import numpy as np
from collections import defaultdict
import csv

from goatools.base import download_go_basic_obo
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.associations import read_gaf

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

# Doing QTL and GO term analysis for the TP (Many strains, few mutations) experiment:

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

# Regions with at least 5 detected QTLs (for different edges)
multi_hit_qtls = [('chr05', 360000, 400000), ('chr08', 60000, 140000), ('chr12', 640000, 720000), ('chr13', 40000, 80000), ('chr14', 360000, 420000), 
                  ('chr14', 430000, 500000), ('chr15', 100000, 200000)]

# My rough notes about them
milo_notes = {
    'qtl_chr08_60000_140000': 'Seems likely that there is some interaction here since we get both SIR proteins, but I cannot find any strong candidates.  Is (near) a YPD qtl in Bloom studies.',
    'qtl_chr12_640000_720000': 'Fitness and adaptability qtl in jerison, found in bloom studies too. Cannot figure any strong candidates.',
    'qtl_chr13_40000_80000': 'Near some qtls in Bloom stress environments, no idea what is going on though.  BUL2 (BUL1 paralog) is in the region...',
    'qtl_chr05_360000_400000': 'This is near the Nat marker locus (FLO8), so it may be spurious, but it is possible it is something else.',
    'qtl_chr15_100000_200000': 'Fitness QTL.  Based on Bloom YPD QTLs nearby it could be an IRA2 thing.  The strongest effects are on the ribosomal proteins, where the epistatic effects look like they do not fit with the pattern (mutations are better on the more fit allele background)',
    'qtl_chr14_430000_500000': 'I think this is the MKT1 QTL, based on the confidences ranges on the Bloom QTLs.  This is a fitness QTL, but its epistatic effect on average does not follow the main pattern we see (mutations are better on the more fit allele background)',
    'qtl_chr14_360000_420000': 'KRE33 QTL. Ribosome-associated. This is the major fitness QTL, and its epistatic effect on average does follow the pattern (mutations are better on the less fit allele background)',
}

def has_qtl(qtl_list, search_qtl):
    if str(qtl_list) != 'nan':
        for qtl in str(qtl_list).split('|'):
            hit = qtl.split(';')[0]
            chromo = hit[:5]
            if chromo == search_qtl[0]:
                hit_loc = int(hit.split('_')[1])
                if hit_loc > search_qtl[1] and hit_loc < search_qtl[2]:
                    return True
    return False

tp_all = pd.read_csv('../../Analysis/TP_data_by_edge.csv')
tp = tp_all.loc[tp_all['Type']=='Experiment']
tp_dfe = pd.read_csv('../../Analysis/TP_DFE_statistics.csv')
segs = [i.split('.')[0] for i in tp if '.mean.s' in i]
gm = get_geno_matrix(segs)

rows = []
for entry in tp.as_matrix(['full.model.qtls', 'full_model_coeffs', 'full_model_qtl_effect_sizes', 'Edge', 'Gene.Use']):
    es = str(entry[2]).split(';')
    coeffs = str(entry[1]).split(';')[2:]
    qtl_list = entry[0]
    if str(qtl_list) != 'nan':
        c = 0
        for qtl in str(qtl_list).split('|'):
            hit, left, right = qtl.split(';')
            chromo = hit[:5]
            qtl_group = ''
            for search_qtl in multi_hit_qtls:
                if chromo == search_qtl[0]:
                    hit_loc = int(hit.split('_')[1])
                    if hit_loc > search_qtl[1] and hit_loc < search_qtl[2]:
                        qtl_group = 'qtl_' + '_'.join([str(i) for i in search_qtl])
            rows.append([hit, left, right, float(coeffs[c]), float(es[c])] + list(entry[3:]) + [qtl_group])
            c += 1
qtl_df = pd.DataFrame(rows, columns=['QTL', 'start_conf_interval', 'end_conf_interval', 'coeff', 'effect_size', 'Edge', 'Gene.Use', 'QTL_group'])
qtl_df.to_csv('../../Analysis/QTL_results.csv', index=False)

# just focusing on the multi-hit QTLs now:
qgroups = [i for i in set(qtl_df['QTL_group']) if str(i) != 'nan']
qgroup_median = {g: np.median([int(i.split('_')[1]) for i in qtl_df.loc[qtl_df['QTL_group']==g]['QTL']]) for g in qgroups}

mhq_mat = []
for mhq in multi_hit_qtls:
    mhq_mat.append(['qtl_' + '_'.join([str(i) for i in mhq]), 
                    ';'.join(list( tp.loc[tp['num.measured']>=50].loc[tp['full.model.qtls'].apply(lambda q: has_qtl(q, mhq))]['Gene.Use'].str.split(' ').str[1])),
                    milo_notes['qtl_' + '_'.join([str(i) for i in mhq])], qgroup_median['qtl_' + '_'.join([str(i) for i in mhq])]])
mhq_dat = pd.DataFrame(mhq_mat, columns=['QTL', 'Genes_with_interactions', 'notes', 'median.loc'])

jerison = pd.read_csv('../../QTL_info_old_papers/elife-27167-supp5-v2.csv') #from jerison et al. 2017
bloom2013 = pd.read_excel('../../QTL_info_old_papers/nature11867-s4.xls') #from bloom et al. 2013
bloom2015_pre = pd.read_excel('../../QTL_info_old_papers/ncomms9712-s4.xls') #from bloom et al. 2015
good_cols = [i for i in bloom2015_pre if not str(list(bloom2015_pre[i])[0])=='nan']
bloom2015 = pd.DataFrame(bloom2015_pre.as_matrix(good_cols)[1:], columns=[list(bloom2015_pre[i])[0] for i in good_cols])

roman_dict = {'I':1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10, 'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16}
roman_back = {roman_dict[i]: i for i in roman_dict}

jerison['q.loc'] = jerison['Locus'].apply(lambda L: 'chr' + str(roman_dict[L.split(' ')[0][3:]]).zfill(2) + '_' + L.split(' ')[1])
bloom2013['q.loc'] = bloom2013.apply(lambda row: 'chr' + str(row['Chromosome']).zfill(2) + '_' + str(row['Peak Position (bp)']), axis=1)
bloom2015['q.loc'] = bloom2015.apply(lambda row:  'chr' + str(roman_dict[row['chr'][3:]]).zfill(2) + '_' + str(row['pos']), axis=1)
bloom2013 = bloom2013.loc[bloom2013['Trait']=='YPD']
bloom2015 = bloom2015.loc[bloom2015['trait']=='YPD']

oq_dats = {'Jerison_2017': jerison, 'Bloom_2013': bloom2013, 'Bloom_2015': bloom2015}

jerison_cols = [ 
       'QTL for trait: fitness at OT', 'QTL for trait: fitness at HT',
       'QTL for trait: adaptability at OT',
       'QTL for trait: adaptability at HT',
       'QTL for trait: adaptability at OT beyond fitness',
       'QTL for trait: adaptability at HT beyond fitness',
       'QTL for trait: pleiotropy at OT beyond fitness',
       'QTL for trait: pleiotropy at HT beyond fitness']

c2g = {'Jerison_2017': ['q.loc'] + jerison_cols, 'Bloom_2013': ['q.loc'], 'Bloom_2015': ['q.loc']}

for d in oq_dats:
    qtl_old_info = dict()
    for qg in multi_hit_qtls:
        qtl = 'qtl_' + '_'.join([str(i) for i in qg])
        td = oq_dats[d].loc[oq_dats[d]['q.loc'].apply(lambda q: has_qtl(q, qg))]
        if len(td) > 0:
            qtl_old_info[qtl] = [[] for c in c2g[d]]
            for entry in td.as_matrix(c2g[d]):
                for i in range(len(entry)): 
                    qtl_old_info[qtl][i].append(entry[i])
        else:
            qtl_old_info[qtl] = [[] for c in c2g[d]]
    for i in range(len(c2g[d])):
        mhq_dat[d + '_' + c2g[d][i]] = mhq_dat['QTL'].apply(lambda q: ';'.join([str(c) for c in qtl_old_info[q][i]]))
        
mhq_dat.to_csv('../../Analysis/Multi_hit_QTLs.csv', index=False)

#GO term analysis, modified from https://github.com/tanghaibao/goatools/blob/master/notebooks/goea_nbt3102.ipynb

# Get http://geneontology.org/ontology/go-basic.obo
obo_fname = download_go_basic_obo()
obodag = GODag("go-basic.obo")
geneid2gos_yeast = read_gaf('../accessory_files/gene_association.sgd') #http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz
genename_2_id = dict()
with open('../accessory_files/gene_association.sgd', 'r') as infile:
    for line in infile:
        if line[0] != '!':
            s = line.split('\t')
            genename_2_id[s[2]] = s[1]
            
id_2_genename = {genename_2_id[i]: i for i in genename_2_id}
ids = [i for i in geneid2gos_yeast.keys()]

all_measured_genes = set(tp.loc[tp['num.measured']>=50]['Gene.Use'].apply(lambda s: s.split(' ')[1]))
background_set = [genename_2_id.setdefault(i, 'NA') for i in all_measured_genes]

goeaobj = GOEnrichmentStudy(
        background_set, # List of all genes in analysis
        geneid2gos_yeast, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = False,
        alpha = 0.05, # default significance cut-off
        methods = ['fdr_bh']) # defult multipletest correction method

# I will also look for GO term enrichment in the group with the top 10 (abs value, but they are all negative) background fitness coeffs in the full model
tp['abs_slope_in_full_model'] = np.abs(tp['full_model_x_slope'])
fit_effect_genes = tp.loc[tp['model_comp_p_full_vs_qtl']<0.05].sort_values(by='abs_slope_in_full_model', ascending=False).iloc[:10]['Gene.Use']

#testing my groups
results = dict()
go_groups = list(mhq_dat.as_matrix(['QTL', 'Genes_with_interactions']))
go_groups.append(['Top_x_effects', ';'.join([s.split(' ')[1] for s in fit_effect_genes])])
for entry in go_groups:
    study_set_names = [i for i in entry[1].split(';') if 'None' not in i]
    study_set = [genename_2_id.setdefault(i, 'NA') for i in study_set_names]
    goea_results_all = goeaobj.run_study(study_set, keep_if=lambda x: x.p_uncorrected < 0.05)
    results[entry[0]] = sorted(goea_results_all, key=lambda r: r.p_fdr_bh)[:10]
        
with open('../../Analysis/GO_results.csv', 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Group', 'GO term', 'pval_uncorrected', 'pval_benjamini_hochberg', 'hits'])
    for r in results:
        for i in range(len(results[r])):
              writer.writerow([r, results[r][i].name, results[r][i].p_uncorrected, results[r][i].p_fdr_bh,
                              ';'.join([id_2_genename[i] for i in results[r][i].study_items])])
  
