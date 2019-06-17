import csv
import pandas as pd


#Parses EditDistanceCluster output to make final read count files
# *** The input file should be ***
# arg1: A csv with barcodes, sorted (probably from high abundance to low abundance)
#  (this MUST be the same file used for the clustering step)
# arg2: The clustering output csv
# arg3: a file base for output (see output types below)
#
# Revising 11/03/2016 - simplifying etc. - no clashing clustering events are possible anymore, so I removed that part



def make_num(s):
    #a simple thing to decide to caste something as an int or float
    if '.' in s:
        return float(s)
    else:
        return int(s)


class ClusterParser:

    def __init__(self, f_in, bc_in, cols_for_total, total_count_thresh):
        self.total_count_thresh = total_count_thresh
        print('\nReading cluster file:', f_in)
        self.raw_rows = []
        with open(f_in, 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                self.raw_rows.append([int(i) for i in row])

        self.raw_rows = sorted(self.raw_rows, key=lambda x: x[0])
        #print('Sorting check: should read 2345 ... ', self.raw_rows[2345][0], '\n')

        print('\nReading barcode count data file:', bc_in)
        self.bcdata = pd.read_csv(bc_in)
        self.bcdata['Reads.For.Clustering'] = self.bcdata.apply(lambda row: sum([row[col] for col in cols_for_total]), axis=1)
        self.cols_to_add = self.bcdata.select_dtypes(include=['int']).keys()
        self.bcdata['Index'] = range(len(self.bcdata)) #make index column
        print('Read', len(self.bcdata), 'rows.')
        self.total_initial_barcodes = len(self.bcdata)

        print('\nMaking clustering dictionaries')
        self.centroids_dict = dict()
        self.error_cor_dict = dict()
        for unsorted_row in self.raw_rows:
            #sort row by total counts in the clustering columns - the one with the most counts is called the 'centroid' here
            row = sorted(unsorted_row, key=lambda x: self.bcdata['Reads.For.Clustering'][x], reverse=True)
            if row[0] in self.error_cor_dict.keys():
                print('Centroid was already clustered that is odd...', row)
            else:
                self.centroids_dict[row[0]] = set(row[1:])
                for index in row[1:]:
                    #note that every centroid corrects to itself
                    self.error_cor_dict[index] = row[0]

        print('\nNot using clusters with <=', self.total_count_thresh, 'reads in the clustering columns.',
              'Starting with', len(self.centroids_dict), 'clusters')
        tmp_cen_dict = dict()
        self.officially_clustered = set(self.centroids_dict.keys()) #A set of all bcs that are clustered
        for centroid in sorted(self.centroids_dict.keys()):
            tmp_count = self.bcdata['Reads.For.Clustering'][centroid]
            for index in self.centroids_dict[centroid]:
                if not index == centroid:
                    #add read counts to centroid row
                    tmp_count += self.bcdata['Reads.For.Clustering'][index]
            if tmp_count > self.total_count_thresh:
                tmp_cen_dict[centroid] = self.centroids_dict[centroid]
                self.officially_clustered.update(self.centroids_dict[centroid])

        #print(self.officially_clustered)
        self.centroids_dict = tmp_cen_dict #update the centroids_dict to only include the ones that passed the test
        print('New number of clusters:', len(self.centroids_dict))

    def write_cluster_fasta(self, fasta_out, bc_index_in_question):
        #Using this you can write a fasta file output to look at a cluster by eye
        #I use .txt just so my mac spotlight will work
        fasta_file = fasta_out + str(bc_index_in_question) + '.txt'
        print('\nOutputting cluster sequence files to:', fasta_out + 'INDEX.txt\n')
        with open(fasta_file, 'w') as outfile:
            centroid_reads = self.bcdata['Reads.For.Clustering'][bc_index_in_question]
            outfile.write('>' + str(bc_index_in_question) + ':' + str(centroid_reads) + ':Centroid\n')
            outfile.write(self.bcdata['Edge.Bases'][bc_index_in_question] + '\n')
            for index in sorted(self.centroids_dict[bc_index_in_question]):
                if not index == bc_index_in_question:
                    #header like >Index:Total Reads
                    outfile.write('>' + str(index) + ':' + str(self.bcdata['Reads.For.Clustering'][index]) + ':' +
                                  str(self.bcdata['Reads.For.Clustering'][index]/centroid_reads) + '\n')
                    outfile.write(self.bcdata['Edge.Bases'][index] + '\n')

    def output_centroids_without_errors(self, f_out):
        print('\nOutputting centroid barcode rows (no errors) to:', f_out)
        self.bcdata.loc[self.bcdata['Index'].isin(self.centroids_dict.keys())].to_csv(f_out, index=False)
        print('Wrote', len(self.bcdata.loc[self.bcdata['Index'].isin(self.centroids_dict.keys())]), 'rows.')

    def output_clusters(self, f_out):
        print('\nOutputting clustered lineages to:', f_out)

        #GOING INTO MATRIX FORM TO MAKE ADDING FASTER
        self.bc_rows = self.bcdata.as_matrix()
        titlerow = [i for i in self.bcdata.keys()]
        columns_to_add = [titlerow.index(col) for col in self.cols_to_add]
        total_count_col = titlerow.index('Total.Count')
        total_reads_clust = 0
        total_bc_clust = 0
        with open(f_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(titlerow)
            for centroid in sorted(self.centroids_dict.keys()):
                tmp_row = self.bc_rows[centroid]
                total_bc_clust += len(self.centroids_dict[centroid])
                for index in self.centroids_dict[centroid]:
                    if not index == centroid:
                        #add read counts to centroid row
                        for i in columns_to_add:
                            tmp_row[i] += self.bc_rows[index][i]
                        #redivide for the umi percent column
                writer.writerow(tmp_row)
                total_reads_clust += tmp_row[total_count_col]

        print('Wrote', len(self.centroids_dict), 'rows.')
        return total_bc_clust, total_reads_clust


    def output_cluster_stats(self, f_out, error_percent_thresh):
        print('\nOutputting stats on clusters to:', f_out)
        total_error_percent_too_high = 0
        with open(f_out, 'w') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Centroid Index', 'Edge.Bases', 'Number of clustered barcodes', 'Number above ' + str(error_percent_thresh*100) + ' % of centroid reads'])
            for centroid in sorted(self.centroids_dict.keys()):
                error_percent_too_high = []
                for index in self.centroids_dict[centroid]:
                    if not index == centroid:
                        if (self.bcdata['Reads.For.Clustering'][index] / self.bcdata['Reads.For.Clustering'][centroid]) > error_percent_thresh:
                            error_percent_too_high.append(index)
                writer.writerow([centroid, self.bcdata['Edge.Bases'][centroid], len(self.centroids_dict[centroid])-1, len(error_percent_too_high)] + error_percent_too_high)
                total_error_percent_too_high += len(error_percent_too_high)

        print('Found', total_error_percent_too_high, 'clustering events in which the clustered barcode had more than',
              error_percent_thresh*100, '% of the total reads of the centroid barcode.')

    def output_unclustered(self, f_out):
        print('\nOutputting unclustered barcodes to:', f_out)
        self.bcdata.loc[~self.bcdata['Index'].isin(self.officially_clustered)].to_csv(f_out, index=False)
        num_unclust_rows = len(self.bcdata.loc[~self.bcdata['Index'].isin(self.officially_clustered)])
        print('Wrote', num_unclust_rows, 'rows.')
        return num_unclust_rows, sum(self.bcdata.loc[~self.bcdata['Index'].isin(self.officially_clustered)]['Total.Count'])
