import pandas as pd
from intermine.webservice import Service
import subprocess
import numpy as np
import csv
import Levenshtein

roman_dict = {
    '01':'I',
    '02':'II',
    '03':'III',
    '04':'IV',
    '05':'V',
    '06':'VI',
    '07':'VII',
    '08':'VIII',
    '09':'IX',
    '10':'X',
    '11':'XI',
    '12':'XII',
    '13':'XIII',
    '14':'XIV',
    '15':'XV',
    '16':'XVI',
    'mt':'mt'
}

CHROMO_DICT = {
    'tpg|BK006937.2|': 'chr03',
    'tpg|BK006938.2|': 'chr04',
    'tpg|BK006948.2|': 'chr15',
    'ref|NC_001224|': 'chrmt',
    'tpg|BK006941.2|': 'chr07',
    'tpg|BK006944.2|': 'chr11',
    'tpg|BK006947.3|': 'chr14',
    'tpg|BK006936.2|': 'chr02',
    'tpg|BK006939.2|': 'chr05',
    'tpg|BK006942.2|': 'chr09',
    'tpg|BK006935.2|': 'chr01',
    'tpg|BK006940.2|': 'chr06',
    'tpg|BK006945.2|': 'chr12',
    'tpg|BK006943.2|': 'chr10',
    'tpg|BK006949.2|': 'chr16',
    'tpg|BK006934.2|': 'chr08',
    'tpg|BK006946.2|': 'chr13',
}

def get_gene2(row):
    # returns "in Gene" if it is in a gene
    # returns "near Gene" if it isn't, and gene is the gene that the insertion is closest to the start of the gene
    if str(row['Gene_ORF']) != 'NA':
        return 'in ' + row['Gene_ORF']
    elif str(row['Gene_ORF.nearby']) != 'NA':
        nearby_list = row['Gene_ORF.nearby'].split('|')
    else:
        return ''
    # if we have a list of nearby genes:
    if len(nearby_list) == 0:
        return 'nearby ' + nearby_list[0]
    else:
        strands = [int(i) for i in row['orf_strand.nearby'].split('|')]
        starts = [int(i) for i in row['start.nearby'].split('|')]
        ends = [int(i) for i in row['end.nearby'].split('|')]
        nearest = None
        nearest_dist = np.inf
        insertion_loc = int(row['insertion_edge'])
        for i in range(len(nearby_list)):
            if strands[i] == -1:
                actual_start = ends[i]
            else:
                actual_start = starts[i]
            if np.abs(actual_start - insertion_loc) < nearest_dist:
                nearest_dist = np.abs(actual_start - insertion_loc)
                nearest = nearby_list[i]
        return 'nearby ' + nearest

def format_bowtie_row(raw_row):
    tmp_dict = {'bowtie_code': raw_row[1], 'chromosome': CHROMO_DICT.setdefault(raw_row[2], 'NA'), 
                'mapq': raw_row[4], 'cigar_match': raw_row[5]}
    # changes insertion location depending on the orientation of alignment
    if raw_row[1] == '16':
        tmp_dict['insertion_strand'] = '-'
        tmp_dict['insertion_edge'] = int(raw_row[3]) + len(raw_row[0]) - 1
    else:
        tmp_dict['insertion_strand'] = '+'
        tmp_dict['insertion_edge'] = int(raw_row[3])
    return tmp_dict

def run_bowtie(seqs, bowtie_base_directory):
    with open('tmp.fasta', 'w') as outfile:
        for seq in seqs:
            outfile.write('>' + seq + '\n' + seq + '\n')
    subprocess.call(['bowtie2', '-x', bowtie_base_directory, '-U', 'tmp.fasta', '-S', 'tmp_bowtie_output.tsv', '-f', '--local'])
    align_info = dict()
    with open('tmp_bowtie_output.tsv', 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        reading = False
        for row in reader:
            if not reading:
                if row[0] == '@PG':
                    reading = True
            else:
                align_info[row[0]] = format_bowtie_row(row)
    return align_info

def gene_orfer(row):
    if str(row['Gene']) == 'None':
        return row['ORF']
    else:
        return row['Gene']

def get_all_gene_annotations():
    service = Service("https://yeastmine.yeastgenome.org:443/yeastmine/service")
    query = service.new_query("Gene")
    col_names = ["briefDescription", "description", "functionSummary",
        "chromosome.primaryIdentifier", "secondaryIdentifier", "symbol",
        "phenotypeSummary", "locations.strand", "locations.end", "locations.start"]
    query.add_view(col_names)
    seen_orfs = set()
    col_dicts = {c: [] for c in col_names}
    for row in query.rows():
        # for some reason rows are repeated in the yeastmine output, so I deduplicate them here
        if row['secondaryIdentifier'] not in seen_orfs:
            for c in col_names:
                 col_dicts[c].append(row[c])
            seen_orfs.add(row['secondaryIdentifier'])
    name_shortener = {
        'chromosome.primaryIdentifier': 'chromosome',
        'secondaryIdentifier': 'ORF',
        'symbol': 'Gene',
        'locations.start': 'start',
        'locations.end': 'end',
        'locations.strand': 'orf_strand'
    }
    td = pd.DataFrame(col_dicts).rename(columns=name_shortener)
    td['Gene_ORF'] = td.apply(lambda row: gene_orfer(row), axis=1)
    return td

def get_genes_hit(chromo, location, ann_df, buf_range):
    chromo_fix = chromo[:3] + roman_dict.setdefault(chromo[3:], 'NA')
    return ann_df.loc[ann_df['chromosome'] == chromo_fix].loc[(ann_df['start'] - buf_range) < location].loc[(ann_df['end'] + buf_range) > location]

def align_and_annotate(edge_set, bowtie_path):
    aligned = run_bowtie(edge_set, bowtie_path)
    ann = get_all_gene_annotations()
    alignment_col_names = ['bowtie_code', 'chromosome', 'insertion_edge', 'insertion_strand', 'mapq', 'cigar_match']
    annotation_col_names = ['briefDescription', 'description', 'functionSummary', 
                            'Gene_ORF', 'phenotypeSummary', 'orf_strand', 'end', 'start']
    dat_dicts = {col: [] for col in ['Edge'] + alignment_col_names + annotation_col_names + [i + '.nearby' for i in annotation_col_names]}
    for edge in edge_set:
        dat_dicts['Edge'].append(edge)
        tmp_align_dict = aligned[edge]
        annotations = get_genes_hit(tmp_align_dict['chromosome'], tmp_align_dict['insertion_edge'], ann, 0)
        ann_nearby = get_genes_hit(tmp_align_dict['chromosome'], tmp_align_dict['insertion_edge'], ann, 500)
        for col in alignment_col_names:
            dat_dicts[col].append(tmp_align_dict[col])
        for col in annotation_col_names:
            if len(annotations) > 0:
                dat_dicts[col].append('|'.join([str(i) for i in annotations[col]]))
            else:
                dat_dicts[col].append('NA')
        for col in annotation_col_names:
            if len(ann_nearby) > 0:
                dat_dicts[col + '.nearby'].append('|'.join([str(i) for i in ann_nearby[col]]))
            else:
                dat_dicts[col + '.nearby'].append('NA')
    return pd.DataFrame(dat_dicts)
    
control_edges = set(['CTAAGTGTGAAGGAGTTGTCTTCTTGCGCT', 'CTGATTTGTGCTGTCTTAGGACCCTCTGAA', 'GCTGCTTATGAGGATATGGATTTAGAGCTA'])
all_edges = set(pd.read_csv('../../BT_Bioinformatic_Work/BT1_output/BT_BC_Assoc/BT_BC_Assoc_filtered_clusters.csv')['Edge'])
all_edges.update(control_edges)
ann_edges = align_and_annotate(all_edges, '../accessory_files/bowtie_build_S288C/S288C')
tncs = list(pd.read_csv('TnCS_Final_Rearray_Data.csv')['Edge.Bases'].str[:30])
ann_edges['Expected.From.TnCS'] = ann_edges['Edge'].isin(tncs)
tp_d = pd.read_csv('../../BT_Bioinformatic_Work/TP_output/TP_BC_Assoc/TP_CS_BC_Assoc_filtered_clusters.csv')[['Edge', 'Total.Counts']].groupby('Edge', as_index=False).sum()
tp = list(tp_d.loc[tp_d['Total.Counts'] > 50000]['Edge'])
ann_edges['In.TP'] = ann_edges['Edge'].str[:15].isin(tp)
ann_edges['TP.Control'] = ann_edges['Edge'].isin(control_edges)
tn96_rearray = pd.read_csv('../accessory_files/Tn96_edges_chosen_final.csv')
tp_chosen = list(tn96_rearray['Edge']) + [i for i in control_edges]
tp_ref = list(tn96_rearray.loc[tn96_rearray['Is.Reference']]['Edge'])
ann_edges['Expected.In.TP'] = ann_edges['Edge'].isin(tp_chosen)
ann_edges['TP.Reference'] = ann_edges['Edge'].isin(tp_ref)
ann_edges['Gene.Use'] = ann_edges.apply(lambda row: get_gene2(row), axis=1)
ann_edges['Edge.ID'] = [i for i in range(len(ann_edges))]
ann_edges.to_csv('../../Mutation_Annotation/Edges_annotated.csv', index=False)
