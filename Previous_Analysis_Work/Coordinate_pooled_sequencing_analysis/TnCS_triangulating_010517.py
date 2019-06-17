#A program to demultiplex reads and count tn insertions for the Tn-seq project
#Milo Johnson
#Modified from SLAMseq_cluster 11/04/2016
#Modified 01/05/2017 (simplified)

import pandas as pd
import csv

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='comma-separated file with edge counts to use as input, with bowtie alignment info')
parser.add_argument('well_id_file', help='intermediate file, gives calls for each well about which insertion is in it')
parser.add_argument('output_file', help='output file')
parser.add_argument('-reads_to_count', type=int, default=2, help='# of reads to count an edge as in that library')
parser.add_argument('-reads_to_ignore', type=int, default=1, help='# of reads that makes us unsure if an edge might be'
                                                                  'in a library - strictest value is 1')

args = parser.parse_args()

###POSITIONAL ARGS
in_csv = args.input_file
well_id_csv = args.well_id_file
output_file = args.output_file

GROUPS = dict()

READS_TO_COUNT = args.reads_to_count
READS_TO_IGNORE = args.reads_to_ignore


def get_group_confidence_and_call(row, group):
    hits = [] #form [hit, reads]
    #looks for reads in any of the columns in the group
    for library in GROUPS[group]:
        reads = row[library]
        if reads > 0:
            hits.append([library, reads])

    if len(hits) == 0:
        return 'none'
    elif len(hits) == 1:
        if hits[0][1] >= READS_TO_COUNT:
            return hits[0][0]
        else:
            return 'none'
    else:
        ambiguous = False
        for hit in hits[1:]:
            if hit[1] >= READS_TO_IGNORE:
                ambiguous = True
        if ambiguous:
            return 'ambiguous'
        else:
            return hits[0][0]

def get_well_identities(Q,PG,R,C,dm):

    hits = []
    for row in dm:
        if 'TnCS_Q' + str(Q) == row[1]:
            if 'TnCS_PG' + str(PG) == row[2]:
                if 'TnCS_' + R == row[3]:
                    if 'TnCS_C' + str(C) == row[4]:
                        hits.append(str(row[0]))
    if len(hits) > 0:
        return ';'.join(hits)
    else:
        return 'none'

def fix_insertion_edges(row, which, offset):
    if row['bowtie_code'] == '16': #reverse strand alignment, need to fix insertion site
        if which == 'strand':
            return '-'
        else:
            return row['insertion_edge']+49
    else:
        if which == 'strand':
            return '+'
        else:
            return row['insertion_edge']


dat = pd.read_csv(in_csv)
dat

TnCS_columns = [i for i in dat.keys() if i[:4] == 'TnCS']
GROUPS['PG'] = [i for i in TnCS_columns if i[:7] == 'TnCS_PG']
GROUPS['C'] = [i for i in TnCS_columns if i[:6] == 'TnCS_C' and len(i) > 6]
GROUPS['Q'] = [i for i in TnCS_columns if i[:6] == 'TnCS_Q']
GROUPS['R'] = [i for i in TnCS_columns if len(i) == 6]
SUMS = {i:sum(dat[i]) for i in TnCS_columns}

print('Triangulating...')
for group in GROUPS.keys():
    #adds column for group hits
    dat[group+'_hit'] = dat.apply(lambda row: get_group_confidence_and_call(row, group), axis=1)


print('Calling wells...')
gs = ['Q', 'PG', 'R', 'C']
dm = dat.as_matrix(['Edge.Bases'] + [g+'_hit' for g in gs])

well_identities = []
plate_counter = 0
well_counter = 0
for PG in range(1,9):
    for Q in range(1,5):
        plate_counter += 1
        for R in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
            for C in range(1, 13):
                well_identities.append(['P' + str(plate_counter) + 'R' + str(R) + 'C' + str(C),
                                        'P' + str(plate_counter).zfill(2) + str(R) +  str(C),
                                        get_well_identities(Q,PG,R,C,dm)])

wells_with_hits = 0
wells_with_single_hit = 0
with open(well_id_csv, 'w') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['Well', 'Well.Simple', 'Hits', 'Good'])
    for row in well_identities:
        row.append(row[2])
        if row[2] != 'none':
            wells_with_hits += 1
            if ';' not in row[2]:
                wells_with_single_hit += 1
                row[3] = 'good'
        writer.writerow(row)

print('Called', wells_with_single_hit, 'out of', len(well_identities), 'wells.', wells_with_hits-wells_with_single_hit,
      'wells had multiple insertions mapped to them. \nCombining data for wells with a single called hit...')
well_dat = pd.read_csv(well_id_csv)

final_calls = well_dat.merge(dat, left_on='Hits', right_on='Edge.Bases')

redundant_columns = ['Hits', 'Good', 'Reads.For.Clustering'] + [g+'_hit' for g in gs]

final_calls[[i for i in list(final_calls.keys()) if i not in redundant_columns]].to_csv(output_file, index=False)





