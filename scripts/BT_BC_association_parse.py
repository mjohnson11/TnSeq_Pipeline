""" 
A program to parse reads from random barcode transposon association sequencing

Can use single-end or paired-end sequencing - if paired-end is given, there is a step to try to make sure 
the tagmentation was not so close to the Tn7 edge that the bases at the edge are mixed with tagmentation bases.
Otherwise, simply looks for and counts barcodes and sequences at the edge of the transposon ("edges").
"""

import gzip
import csv
import regex
from milo_tools import FourLineFastq
from collections import defaultdict, Counter
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input_base', help='input fastq data directory')
parser.add_argument('output_base', help='output data directory base')
parser.add_argument('library_map_file', 
                    help='file with library names, file names, and primer offsets')
parser.add_argument('-genomic_bp_cutoff', type=int, default=15, help='number of bases to read after the edge of the transposon')
parser.add_argument('-paired_end', action='store_true')
args = parser.parse_args()

# Library map file input should have three comma-separated columns: library name, primary file names, and primer offsets
# The primary file names can be ; separated if there are multiple, and they should correspond to the reads starting inside Tn7
# If it is paired-end sequencing, the seconday files must be the same with R1 and R2 switched (primary does not have to be R1)

# POSITIONAL ARGS
input_base = args.input_base
output_base = args.output_base
library_map_file = args.library_map_file
GENOMIC_BP_CUTOFF = int(args.genomic_bp_cutoff)
paired_end = args.paired_end

# POSITIONS FOR SEARCHING FOR STUFF
BC_RANGE_START = 9
BC_RANGE_END = 56
# a list of increasingly lenient regular expressions to extract the barcode
BC_REGEX_LIST = [regex.compile('(GTAGAA)(\D{28})(GCGTCA)'),
                 regex.compile('(GTAGAA)(\D{26,30})(GCGTCA)'),
                 regex.compile('(GTAGAA){e<=1}(\D{28})(GCGTCA){e<=1}'),
                 regex.compile('(GTAGAA){e<=1}(\D{26,30})(GCGTCA){e<=1}')]

EDGE_SEQ = 'TCCGCCCACA'
EDGE_SEQ_REGEX = regex.compile('(TCCGCCCACA){e<=1}')
EDGE_RANGE_START = 108-len(EDGE_SEQ)-5

# Because I have had issues with fragments that are too small (tagmentation was too close to the Tn7 edge)
# I have various checks to try to avoid recording "edges" that have tagmentation bases in them
# These are the first 60 bp of Tn7 that will be read by the "secondary" read.  
# If I have a secondary read, I will use it to tell how many bp of genomic dna I have.
EDGE_SEQ_RT = 'TGTGGGCGGACAAAAAAGTCTCAAACTGGACAAAAAAGATCCGTAACGTAGGTCTCTGAC'

def get_len_genomic_bp(read):
    """function to find the number of bp between the tn7 edge and tagmentation"""
    len_gbp = Counter()
    # searching for 6 10bp sequences
    for i in range(6):
        short_seq = EDGE_SEQ_RT[i*10:(i+1)*10]
        if short_seq in read:
            len_gbp[read.index(short_seq) - i*10] += 1
    if len(len_gbp) > 0:
        # if there were hits, and at least 2 hits agreed, return the most common hit
        best_hit = len_gbp.most_common(1)[0]
        if best_hit[1] > 1:
            return best_hit[0]
    # otherwise we'll assume the fragment is large
    return 150

# main data structure, basically a dictionary where keys are bc_edge and values are dicts
# with keys for different libraries and counts as values.
# the defaultdict just makes it so I don't have to initialize anything.
cd = defaultdict(Counter)

stats_rec = []
stat_titles = ['Reads', 'Tagmentation.Too.Close', 'BC.Found', 'Edge.Found', 'Usable']

def do_it(dem, primary_files, offset):
    # because we use primers with variable lengths, we include the 'offset' variable to change the positions where
    # we look for the barcode and the edge
    st = {stat_title: 0 for stat_title in stat_titles}
    st['Library'] = dem
    for primary_read_file in primary_files:
        print(dem, primary_read_file)
        # For sequencing reasons, some of my libraries have R1 as the primary read and some have R2 instead
        if paired_end:
            if '.R1.fastq' in primary_read_file:
                secondary_read_file = primary_read_file.replace('.R1.fastq', '.R2.fastq')
            else:
                secondary_read_file = primary_read_file.replace('.R2.fastq', '.R1.fastq')
            r2in = gzip.open(secondary_read_file, 'rt')
            r2_generator = FourLineFastq(r2in)
        r1in = gzip.open(primary_read_file, 'rt')
        for r1_title, r1_seq, r1_qual in FourLineFastq(r1in):
            st['Reads'] += 1
            tag_too_close = False
            if paired_end:
                r2_title, r2_seq, r2_qual = next(r2_generator)
                # checking if tagmentation was too close to the edge
                if get_len_genomic_bp(r2_seq) < GENOMIC_BP_CUTOFF:
                    st['Tagmentation.Too.Close'] += 1
                    tag_too_close = True
            if not tag_too_close:
                bc, genome = None, None
                # looking for barcode
                for reg_search in BC_REGEX_LIST:
                    reg_hit = reg_search.search(r1_seq[BC_RANGE_START + offset:BC_RANGE_END + offset])
                    if reg_hit:
                        bc = reg_hit.group(2)
                        break
                if bc:
                    st['BC.Found'] += 1
                    edge_onward = r1_seq[EDGE_RANGE_START + offset:]
                    # looking for transposon edge
                    if EDGE_SEQ in edge_onward[:20]:
                        genome = edge_onward[edge_onward.index(EDGE_SEQ)+len(EDGE_SEQ):len(edge_onward)]
                    else:
                        edge_reg = EDGE_SEQ_REGEX.search(edge_onward[:20])
                        if edge_reg:
                            genome = edge_onward[edge_reg.end():len(edge_onward)-1]
                    if genome:
                        st['Edge.Found'] += 1
                        edge = genome[:GENOMIC_BP_CUTOFF]
                        if len(edge) == GENOMIC_BP_CUTOFF:
                            combined_bc_edge = bc + '_' + edge
                            cd[combined_bc_edge][dem] += 1
                            st['Usable'] += 1
        r1in.close()
        if paired_end:
            r2in.close()

    print(dem, ':', st['Reads'], 'reads.')
    print(st['Tagmentation.Too.Close'], 'reads excluded for tagmentation being too close to the tn edge.')
    print(st['BC.Found'], 'found the bc by regex.', st['Edge.Found'], 'found the tn edge.')
    print(st['Usable'], 'reads were used.')
    stats_rec.append(st)

        
### BEGIN MAIN ###
file_info = pd.read_csv(library_map_file)
libs_in_order = list(file_info['Library'])
lib_info = {r[0]: r[1:] for r in file_info.as_matrix(['Library', 'File(s)', 'Offset'])}

print('To count and split:', ', '.join(libs_in_order), '\n')
for lib1 in libs_in_order:
    do_it(lib1, [input_base + i for i in lib_info[lib1][0].split(';')], int(lib_info[lib1][1]))

# OUTPUT STATS
with open(output_base + 'read_stats.csv', 'w') as outstats:
    writer = csv.writer(outstats)
    cols = ['Library'] + stat_titles
    writer.writerow(cols)
    for stat in stats_rec:
        writer.writerow([stat[c] for c in cols])

# TALLY TOTALS AND OUTPUT DATA
for cbe in cd:
    td = cd[cbe]
    td['Total.Counts'] = sum([td[k] for k in td])

with open(output_base + 'counts.csv', 'w') as outcounts:
    writer = csv.writer(outcounts)
    tcols = ['BC', 'Edge', 'Total.Counts'] + libs_in_order
    writer.writerow(tcols)
    for cbe in sorted([k for k in cd], key=lambda x: cd[x]['Total.Counts'], reverse=True):
        td = cd[cbe]
        tmp_row = cbe.split('_')
        for col in tcols[2:]:
            # not using defaultdict simplicity here to save memory by not storing tons of zeros
            if col in td:
                tmp_row.append(td[col])
            else:
                tmp_row.append(0)
        writer.writerow(tmp_row)

print('Done!')
