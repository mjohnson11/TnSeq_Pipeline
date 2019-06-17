"""
Simulating data from the TP experiment.  I'm using parameters from a first-pass estimation to make it as realistic as possible:
Mean of the DFE: mean: -0.026, std dev: 0.012 (going to be a random variable)
Std Dev of the DFE: mean: 0.041 (just going to use this, not going have it be an RV)
# of bcs: big range, not normal at all.  I will just do an exponential distribution with an expectation of 1200.  For inital frequency I will just do a multinomial draw for now.
# reads / tp: mean: ~120,000, std dev: ~67000 (going to be a random var)
"""

import pandas as pd
import numpy as np

def wf_step(f, s, N):
    # Wright-Fisher multinomial sampling
    f_alt = f*(np.exp(s)) # applying s values
    f_use = f_alt / sum(f_alt)
    new_counts = np.random.multinomial(N, f_use)
    return new_counts / sum(new_counts)
    
def simulate_evolution(initial_frequencies, s_vec, N_list):
    # Simulating evolution according to WF dynamics
    # The only weird thing here is that I could the bottleneck as a step, 
    # so the indices don't exactly correspond to generations
    f_rec = []
    current_f = initial_frequencies
    neut_s = np.array([0]*len(s_vec))
    for i in range(len(N_list)):
        if i == 0 or (i > 0 and (N_list[i] < N_list[i-1])):  #bottleneck step, no s effects here
            current_f = wf_step(current_f, neut_s, N_list[i])
        else:
            current_f = wf_step(current_f, s_vec, N_list[i])
        f_rec.append(current_f)
    return f_rec

def do_sim(bottle, real_edge_s_vals, outname):
    outlier_frequency = np.random.uniform(0, 0.1) # % of outliers ranges from 0% to 10%
    num_bcs = int(np.random.exponential(1200))
    initial_freqs = np.random.multinomial(bottle, np.array([1/num_bcs]*num_bcs))
    bc_edge_assignments = np.random.choice(np.arange(100), size=num_bcs)
    bc_true_edge_s = np.array([real_edge_s_vals[bc_edge_assignments[i]] for i in range(num_bcs)])
    bc_s = []
    for i in range(num_bcs):
        if np.random.random() < outlier_frequency:
            bc_s.append(np.random.uniform(-0.1, 0.1))
        else:
            bc_s.append(bc_true_edge_s[i])
            
    # note that I'm counting the bottleneck step as a generation, so each cycle is 11 generations here
    N_L = [bottle*(2**i) for i in range(11)]*5
    sampling_times = [i for i in range(len(N_L)) if N_L[i]==bottle*(2**10)]
    f_rec = simulate_evolution(initial_freqs, np.array(bc_s), N_L)
    # Now drawing reads from those frequencies, # of reads is normally distributed RV, but I cut it off so it can't go below 1000
    # In my analysis tps with <5000 reads are excluded, so this data won't be used anyways
    reads = [np.random.multinomial(int(np.clip(np.random.normal(120000, 67000), 1000, 1000000)), f_rec[i]) for i in sampling_times]
    full_array = [[int(i) for i in bc_edge_assignments], [i for i in range(num_bcs)]] + reads
    td = pd.DataFrame(np.array(full_array).T, columns=['Edge', 'BC'] + ['sim-T'+str(i) for i in range(5)], dtype=int)
    td['Edge.s'] = bc_true_edge_s
    td['Used.s'] = bc_s
    td.to_csv(outname, index=False)
            
def do_seg(b_neck, simnum):
    dfe_mean = np.random.normal(-0.026, 0.012)
    dfe_std_dev = 0.041
    edge_s_vals = np.concatenate([np.zeros(5), np.random.normal(dfe_mean, dfe_std_dev, 95)])
    do_sim(b_neck, edge_s_vals, 'sim/TP-sim' + str(simnum) + '_r1_counts.csv')
    do_sim(b_neck, edge_s_vals, 'sim/TP-sim' + str(simnum) + '_r2_counts.csv')
    
# Doing simulations for 100 fake segregants
for sn in range(1, 101):
    do_seg(200000, sn)

