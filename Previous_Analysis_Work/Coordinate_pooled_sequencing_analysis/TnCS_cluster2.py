#A program to demultiplex reads and count tn insertions for the Tn-seq project
#Milo Johnson
#Modified from SLAMseq_cluster 11/04/2016
#Modified 01/05/17

import csv
from SLAMseq_DeletionClustering import DeletionClusterator
from TnCS_ClusterParser2 import ClusterParser
import pandas as pd
import subprocess

#the argument tells us which primer index pair to use to demultiplex
#it indicates the row of the primer_map_file to use, and can be used for jabba rays


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('input_file', help='comma-separated file with bc counts to use as input')
parser.add_argument('output_file', help='intermediate .csv file with bc clustering information (each row goes centroid, clustered_bc, clustered_bc, etc.)')
parser.add_argument('output_file_base', help='output file base directory for final clustered counts and clustering statistics')
parser.add_argument("-bowtie_base", type=str, default='none', help="bowtie genome file base directory, if not included, bowtie will not run,"
                                                                   "i.e. bowtie_build_S288C/S288C ")
parser.add_argument("-edit_dist_thresh", type=int, default=3, help='threshold for edit distance clustering within UMI families, default=3')
parser.add_argument("-centroid_count_thresh", type=int, default=8, help='threshold of counts for a centroid to be considered real, default=8')
parser.add_argument("-cluster_count_thresh", type=int, default=8, help='threshold of counts for a cluster to appear in the final dataset, default=8')
parser.add_argument("-cols_for_total", type=str, default='Total.Count', help='commma-separated list of column headers of'
                                                                             'the count data columns to be used for clustering,'
                                                                             'so the sum of the reads in these columns will be used'
                                                                             'in reference to the count thresholds, and the list of '
                                                                             'bcs will be sorted based on these columns, '
                                                                             'default= Total.Count')
args = parser.parse_args()

###POSITIONAL ARGS
#input file: should have columns with headers that say "Edge.Bases", and "Total.Count" (followed by additional info)
in_csv = args.input_file
#cluster output file
out_csv = args.output_file
#output base for cluster info
out_base = args.output_file_base

###OPTIONAL ARGS

BOWTIE_BASE = args.bowtie_base

#max edits for errors to cluster
EDITDIST_THRESH = args.edit_dist_thresh

#potential centroids with less than or equal to this number of total reads are not included
CENTROID_COUNT_THRESH = args.centroid_count_thresh

#clusters with less than or equal to this many reads will not be included
CLUSTER_THRESH = args.cluster_count_thresh

COLS_FOR_TOTAL = args.cols_for_total.split(',')

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


def output_seq_correction_info(ec, full_data_list):
    # going through the original file and correcting, outputting a file with row indices like:
    # row index, row index, row index, etc.

    clusters_by_index = dict()

    for r in range(len(full_data_list)):
        row = full_data_list[r]
        dna_region = row[0]
        if dna_region in ec.corrector.keys():
            real_region = ec.corrector[dna_region]
            if real_region in clusters_by_index.keys():
                clusters_by_index[real_region].append(r)
            else:
                clusters_by_index[real_region] = [r]


    print('Writing error correction cluster output to:', out_csv)
    with open(out_csv, 'w') as outfile:
        writer = csv.writer(outfile)
        for bc_combo in sorted(clusters_by_index.keys(), key = lambda x: clusters_by_index[x][0]):
            writer.writerow(clusters_by_index[bc_combo])

    print('Parsing these clusters now with ClusterParser...')

    main_out = out_base + '_counts_clustered.csv'
    cluster_details = out_base + '_cluster_info.csv'
    unclustered_out = out_base + '_unclustered.csv'
    my_clus_pars = ClusterParser(out_csv, in_csv, COLS_FOR_TOTAL, CLUSTER_THRESH)
    unclustered_bcs, unclustered_reads = my_clus_pars.output_unclustered(unclustered_out)
    clustered_bcs, clustered_reads = my_clus_pars.output_clusters(main_out)
    my_clus_pars.output_cluster_stats(cluster_details, 0.01)

    print("Looked at", my_clus_pars.total_initial_barcodes, "dna regions.")
    print("Made", len(my_clus_pars.centroids_dict.keys()), "clusters from", clustered_bcs, "bcs, totaling", clustered_reads, 'reads.')
    print("Could not match", unclustered_bcs, "bcs, totaling", unclustered_reads, 'reads.')



class BowtieAction:

    def __init__(self, clustered_data_file, bowtie_out):
        self.bowtie_out = bowtie_out
        self.clustered_data_in_order = []
        rc = 0
        with open(clustered_data_file, 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                rc += 1
                if rc == 1:
                    self.data_title_row_tmp = row
                    self.position_of_index = row.index('Index')
                else:
                    self.clustered_data_in_order.append(row)

        with open('tmp.fasta', 'w') as outfile:
            for row in self.clustered_data_in_order:
                outfile.write('>' + row[self.position_of_index] + '\n' + row[0] + '\n')

        subprocess.call(['bowtie2', '-x', BOWTIE_BASE, '-U', 'tmp.fasta', '-S', self.bowtie_out, '-f', '--local'])

    def parse_output(self):
        self.index_to_bowtie_out = dict()
        with open(self.bowtie_out, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            reading = False
            for row in reader:
                if not reading:
                    if row[0] == '@PG':
                        reading = True
                else:
                    #fixing chromosome
                    if not row[2] == '*':
                        row[2] = CHROMO_DICT[row[2]]
                    #fixing insertion location based on bowtie code which indicates the strand direction
                    if row[1] == '16':
                        row[1] = '-'
                        row[3] = str(int(row[3])+len(row[9])-1)
                    else:
                        row[1] = '+'
                    self.index_to_bowtie_out[row[0]] = row


    def write_full_output(self, f_out):
        with open(f_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(self.data_title_row_tmp + ['strand', 'chromosome', 'insertion_edge', 'mapq', 'cigar_match'])
            for row in self.clustered_data_in_order:
                tmp_index = row[self.position_of_index]
                writer.writerow(row + self.index_to_bowtie_out[tmp_index][1:6])


#
#MAIN CALLS
#


data = pd.read_csv(in_csv)
data['Reads.For.Clustering'] = data.apply(lambda row: sum([row[col] for col in COLS_FOR_TOTAL]), axis=1)
use_list = data.as_matrix(['Edge.Bases', 'Reads.For.Clustering'])
print('Correcting edge base sequences')
seq_corrector = DeletionClusterator(use_list)
seq_corrector.get_deletion_clusters(True)
seq_corrector.cluster_within_delnets(EDITDIST_THRESH, CENTROID_COUNT_THRESH, True)
output_seq_correction_info(seq_corrector, use_list)

if not BOWTIE_BASE == 'none':
    bt = BowtieAction(out_base + '_counts_clustered.csv', out_base + '_bowtie_output.tsv')
    bt.parse_output()
    bt.write_full_output(out_base + '_clustered_with_alignments.csv')