import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from simple_s_estimation_v2 import s_estimation
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('demult_id_key', help='key to identify which file to parse (which row of the demult file) - for jabba-rays (job arrays)')
args = parser.parse_args()

demult_id_key = int(args.demult_id_key)

ll_cutoff = 15 # log-likelihood ratio cutoff for excluding outliers
bc_target_num = 7 # maximum # of bcs to reduce to (see function in simple_s_estimation)

experiment = 'BT'
input_base = '../BT_Bioinformatic_Work/BT1_output/BT1_BFA/'
output_base = 'BT_output/'
replicate_info_file = 'accessory_files/BT_segregant_replicate_info.csv'
excluded_segs = list(pd.read_csv('accessory_files/BT_segregants_excluded.csv')['segregant'])

# Reading extra info
rep_d = pd.read_csv(replicate_info_file)
rep_info = {i[0]: i[1:] for i in rep_d.as_matrix(['segregant', 'Replicate_1', 'Replicate_2'])}
segs = list(rep_d['segregant'])

s_estimation([s for s in segs[demult_id_key*5:(demult_id_key+1)*5] if s not in excluded_segs], rep_info, output_base, input_base, experiment, ll_cutoff, bc_target_num, [], consider_all_edges_neutral=True)
