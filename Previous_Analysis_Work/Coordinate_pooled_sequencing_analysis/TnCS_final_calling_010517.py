#Made 01/05/2017

import pandas as pd
import re
import numpy as np

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='comma-separated file with called wells and edge info')
parser.add_argument('genome_in', help='location of S288C chromosomes (+ chr11.fsa for example)')
parser.add_argument('filter_file', help='intermediate file, gives filtering decision for each called well')
parser.add_argument('output_file', help='output file')
parser.add_argument('-sali_range', type=int, default=800, help='bp range for exclusion based on being close to a SALI site')

args = parser.parse_args()

###POSITIONAL ARGS
in_csv = args.input_file
genome_in = args.genome_in
filter_file = args.filter_file
output_file = args.output_file
sali_range = args.sali_range

def find_all_substrings(big_str, search_str):
    #from http://stackoverflow.com/questions/4664850/find-all-occurrences-of-a-substring-in-python
    return [m.start() for m in re.finditer(search_str, big_str)]

def is_rdna(row):
    r = 0
    if row['chromosome'] == 'chr12':
        if row['insertion_edge'] > 445000 and row['insertion_edge'] < 495000:
            r = 1
    return r

def near_sali(row, sali_d, cut_range):
    chromo = row['chromosome']
    location = row['insertion_edge']
    if chromo in sali_d.keys():
        for cut_site in sali_d[chromo]:
            if abs(location-cut_site) < cut_range:
                return 1
    return 0

#FIND ALL SALI sites in genome
SALI_site = 'GTCGAC'
chrs = ['chr01', 'chr02', 'chr03', 'chr04', 'chr05', 'chr06', 'chr07', 'chr08', 'chr09', 'chr10', 'chr11', 'chr12',
        'chr13', 'chr14', 'chr15', 'chr16', 'chrmt']

SALI_dict = dict()
for chr in chrs:
    with open(genome_in + chr + '.fsa', 'r') as infile:
        rc = 0
        sequence = ''
        for line in infile:
            rc += 1
            if rc > 1:
                if line[len(line)-1] == '\n':
                    sequence += line[:len(line)-1]
                else: #for last line
                    sequence += line
    SALI_dict[chr] = find_all_substrings(sequence, SALI_site)

well_calls = pd.read_csv(in_csv)

well_calls['rDNA_loose'] = well_calls.apply(lambda row: is_rdna(row), axis=1)
well_calls['sali_likely'] = well_calls.apply(lambda row: near_sali(row, SALI_dict, sali_range), axis=1)
mapq_filt = well_calls.loc[well_calls['mapq'] == 44]
rDNA_filt = mapq_filt.loc[mapq_filt['rDNA_loose'] == 0]
sali_filt = rDNA_filt.loc[rDNA_filt['sali_likely'] == 0]
print('Looked at', len(well_calls), 'wells.', len(mapq_filt), 'aligned uniquely in the genome.\n',
      len(rDNA_filt), 'after excluding insertions near the rDNA array.', len(sali_filt),
      'were also not cut by SalI.')

well_calls.to_csv(filter_file, index=False)

possible_wells = []
for plate_counter in range(1, 13):
    for R in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:
        for C in range(1, 13):
         possible_wells.append('P' + str(plate_counter).zfill(2) + str(R) +  str(C))

print('So there will be', len(possible_wells)-len(sali_filt), 'blanks in the 12 plates after re-array.')

blanks = np.random.choice(possible_wells, len(possible_wells)-len(sali_filt), replace=False)

sali_filt.sort('Well')

print('blanks:', ','.join(blanks))

sali_filt['Rearray.Well'] = [i for i in possible_wells if not i in blanks]

sali_filt[['Rearray.Well', 'Well.Simple', 'Edge.Bases', 'Total.Count', 'Example.Full.Read', 'Index', 'strand',
           'chromosome', 'insertion_edge', 'mapq', 'cigar_match']].to_csv(output_file, index=False)