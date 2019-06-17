
# DeletionClusterator is a class for clustering barcodes by linking single deletions that overlap into networks

import Levenshtein
import numpy as np


def get_deletion_neighborhood(stringer):
    return set([stringer] + [stringer[:x] + stringer[x+1:] for x in range(len(stringer))])

class DeletionClusterator:

    def __init__(self, list_of_bcs):
        #input: list of [bc, count]
        self.entries = dict()
        #adds counts for identical bcs (these exist because there might be the same env bc with different div bcs or vice versa)
        for entry in list_of_bcs:
            bc = entry[0]
            if bc in self.entries.keys():
                self.entries[bc][0] += entry[1]
            else:
                self.entries[bc] = [entry[1]]


    def get_deletion_clusters(self, verbose):
        # Makes clusters of barcodes linked by common single deletions
        # this accounts for single edits - deletions, insertions, and substitions
        # (not inversions, but that doesn't seem relevant here)
        # these clusters are kept in the self.hood_lists (a dictionary with indices that access the lists)

        self.hood_lists = dict()
        # this is similar to hood_lists, but includes all deletions and the values are sets
        # it is needed for merges, when all these deletions need to be re-pointed at the merge-to hood
        self.hood_all_dels = dict()
        # dictionary that points deletions to hoods (indices)
        self.del_to_hood = dict()
        # list of merged indices
        merged = []
        # counter for hood indices
        hood_count = 0
        # sort bcs by count
        sorted_bcs = sorted(self.entries.keys(), key=lambda x: self.entries[x][0], reverse=True)
        for bc in sorted_bcs:
            entry = [bc, self.entries[bc][0]]
            del_net = get_deletion_neighborhood(bc)
            hoods_hit = set()
            for d in del_net:
                if d in self.del_to_hood.keys():
                    hoods_hit.add(self.del_to_hood[d])
            if len(hoods_hit) == 0:
                #new hood
                hood_count += 1
                self.hood_lists[hood_count] = [entry + ['peak']]
                self.hood_all_dels[hood_count] = del_net
                use_hood = hood_count
            elif len(hoods_hit) > 1:
                #multiple hits, merge hoods
                hood_to_merge_to = min(hoods_hit)
                for hood in hoods_hit:  #done in a loop in case it is >2 hoods
                    if not hood == hood_to_merge_to:
                        for d in self.hood_all_dels[hood]:
                            #point each deletion to the merge-to hood
                            self.del_to_hood[d] = hood_to_merge_to
                        #combine hood lists and all_dels sets, then delete the ones that got merged
                        self.hood_lists[hood_to_merge_to] += self.hood_lists[hood]
                        self.hood_all_dels[hood_to_merge_to].update(self.hood_all_dels[hood])
                        del self.hood_lists[hood]
                        del self.hood_all_dels[hood]
                        merged.append(hood)
                use_hood = hood_to_merge_to
                #add current entry to hood list and its del_net to all_dels set
                self.hood_lists[use_hood].append(entry + ['not_peak'])
                self.hood_all_dels[use_hood].update(del_net)
            else:
                use_hood = hoods_hit.pop() #takes the only one in the set
                #add current entry to hood list and its del_net to all_dels set
                self.hood_lists[use_hood].append(entry + ['not_peak'])
                self.hood_all_dels[use_hood].update(del_net)

            #assign dels to the matched hood
            for d in del_net:
                self.del_to_hood[d] = use_hood
        if verbose:
            print('Made', len(self.hood_lists), 'deletion neighborhoods from', len(self.entries), 'barcodes.')



    def cluster_within_delnets(self, max_edits, centroid_thresh, verbose):
        # For each deletion neighboorhood, this function clusters barcodes based explicitly on edit distance
        # centroids are only allowed when they have at least centroid_thresh reads (just the centroid, no errors)
        # and when they have no single-deletion neighbors with more counts than them
        # then bcs in the hood are clustered to these centroid if they are within max_edits away
        # if not, and if they pass the read threshold, they become a new centroid in that hood
        # if any bc matches multiple centroids in its hood, that is recorded

        self.corrector = dict()
        singletons = 0 #counts only singletons that pass the threshold and become real clusters
        centroid_cluster_sizes = dict()
        excluded_clusters = 0  #excluded because the most common bc in the hood didn't have enough reads
        excluded_others = 0  #excluded because it was out of edit distance range of centroids and didn't have enough reads to be a centroid
        matched_mult_centroids = 0 #counter for when a bc matches multiple centroids in its del hood
        peaks_clustered = 0
        for hood in sorted(self.hood_lists.keys()):

            #sorts the entries in a deletion hood from high counts to low
            entries = sorted(self.hood_lists[hood], key=lambda x: x[1], reverse=True)

            if len(entries) == 1:
                bc = entries[0][0]
                reads = entries[0][1]
                if reads <= centroid_thresh:
                    excluded_clusters += 1
                else:
                    centroid_cluster_sizes[bc] = reads
                    self.corrector[bc] = bc
                    singletons += 1
            else:
                #CRITERIA FOR A CENTROID:
                # 1) has no single-deletion neighbors with more counts (is a "peak")
                # 2) has more than CENTROID_COUNT_THRESH reads
                # 3) is more than max_edits away from centroids with more reads
                cluster_is_cool = True
                peaks = [i for i in entries if (i[2] == 'peak' and i[1] > centroid_thresh)]
                centroids = []
                for i in range(len(peaks)):
                    peak_is_good = True
                    within_range = []
                    for j in range(len(peaks)):
                        if not i == j:
                            if Levenshtein.distance(peaks[i][0], peaks[j][0]) <= max_edits:
                                within_range.append(j)
                    if len(within_range) > 0:
                        if min(within_range) < i: #this peak is within clustering range of a higher peak
                            peak_is_good = False
                            peaks_clustered += 1
                            #print('Clustered a peak:')
                            #print(peaks[i])
                            #print(peaks[min(within_range)])
                    if peak_is_good:
                        centroids.append(peaks[i])

                for entry in centroids:
                    self.corrector[entry[0]] = entry[0]
                    centroid_cluster_sizes[entry[0]] = entry[1]
                for entry in entries:
                    if not entry in centroids:
                        tmp_bc, tmp_reads, peak_jnk = entry
                        matches_centroids = []
                        for centroid in centroids:
                            dist = Levenshtein.distance(centroid[0], tmp_bc)
                            if dist <= max_edits:
                                matches_centroids.append(centroid[0])
                        if len(matches_centroids) == 0:
                            #wasn't within edit distance of the centroids, no good
                            excluded_others += 1
                        else:
                            if len(matches_centroids) > 1:
                                matched_mult_centroids += 1
                                cluster_is_cool = False
                            #in the case of more than one matches, use the one with the most reads
                            #this also handles the normal, one match case
                            self.corrector[entry[0]] = matches_centroids[0]
                            centroid_cluster_sizes[matches_centroids[0]] += entry[1]

                if not cluster_is_cool:
                    print('Centroids with overlapping cluster-ers, here are the multiple centroids:')
                    for c in centroids:
                        print(c)

        if verbose:
            print('Out of', len(self.hood_lists), 'deletion neighborhoods we made', len(centroid_cluster_sizes), 'clusters.')
            print('There were', singletons, 'singleton neighborhoods, we excluded', excluded_clusters, 'clusters based on a centroid read threshold')
            print(excluded_others, 'bcs were out of the edit-distance range and were not centroids and', matched_mult_centroids, 'matched multiple centroids')