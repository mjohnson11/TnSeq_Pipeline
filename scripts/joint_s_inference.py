"""
A class for calculating s for BFAs.  Uses a multinomial model to jointly infer s values.

The 0-index / first lineage in the reads is taken as the reference and is assigned s=0.
"""
import numpy as np
from scipy import optimize

class BfaParamEstimator:
    def __init__(self, R, lineage_info, lineage_column_names, tps, record_process=False):
        # the first lineage will be used as the reference
        self.R = np.copy(R) # Read counts
        self.L = len(R[0]) # Number of lineages
        self.N = np.sum(self.R, axis=1) # read sums by timepoint
        self.lineage_info = lineage_info # identifying lineage information
        self.lineage_column_names = lineage_column_names # column names of identifying info
        # timepoints are shifted so that the lowest one is 0
        self.tps = np.array([t-min(tps) for t in tps])
        self.i_s, self.i_fnot = self.get_initial_params()
        self.current_slopes = None
        self.s_standard_errors = np.zeros(self.L-1)
        self.converged = -1
        self.iterations = 0
        self.record_process = record_process
        if record_process:
            self.ss_rec = []  # records step sizes at each step
            self.ll_rec = []  # records log-likelihoods at each step
            self.s_diffs = []  # records mean and max s changes at each step
            self.fnot_diffs = []  # records mean and max fnot changes at each step
        # get initial parameters
        self.start_params = np.concatenate([self.i_s, self.i_fnot])
        self.current_params = np.copy(self.start_params)
        self.current_ll = self.ll_simple(self.current_params)

    def get_initial_params(self):
        # fnot guess = log(reads at t0 / reference (first lineage) reads at t0)
        start_fnot = np.clip(np.log(self.R[0]/self.R[0][0])[1:], np.log(1/self.R[0][0]), 100)
        s_guesses = []
        for i in range(1, self.L):
            tmp_s_vals = []
            for t in range(1, len(self.tps)-1):
                if self.R[t][i] > 3 and self.R[t][0] > 3:
                    tmp_s_vals.append((np.log(self.R[t][i]/self.R[t][0]) - start_fnot[i-1]) / self.tps[t])
            if len(tmp_s_vals) > 0:
                s_guesses.append(np.mean(tmp_s_vals))
            else:
                s_guesses.append(0)
        start_s = np.clip(s_guesses, -0.5, 0.5)
        return start_s, start_fnot

    def deterministic_freqs(self, params):
        s = np.concatenate([np.zeros(1), np.array(params[:self.L-1])])
        # fnot is clipped at to avoid crazy frequencies
        fnot = np.clip(np.concatenate([np.zeros(1), np.array(params[self.L-1:])]), -100, 100)
        f_unscaled = np.exp(fnot + s*np.tile(self.tps, (self.L, 1)).T)
        f_unscaled_sums = np.sum(f_unscaled, axis=1)
        return (f_unscaled.T/f_unscaled_sums).T

    def ll_simple(self, params):
        # returns the negative of the log-likelihood
        assert len(params) == 2*self.L - 2
        # calculating frequencies using a deterministic equation
        f = self.deterministic_freqs(params)
        # multinomial log probability of reads given frequencies
        ll = np.nansum(self.R*np.log(f), axis=(0,1))
        return -ll

    def ll_derivative(self, params, all_s_neutral=False):
        assert len(params) == 2*self.L - 2
        f = self.deterministic_freqs(params)
        # calculating partial derivatives
        Nf = (f.T*self.N).T
        fnot_deriv_terms =  (self.R - Nf).T
        fnot_derivs = np.sum(fnot_deriv_terms, axis=1)
        if all_s_neutral:
            s_derivs = np.zeros(len(fnot_derivs))
        else:
            s_derivs = np.sum(fnot_deriv_terms*self.tps, axis=1)
        return np.concatenate([s_derivs[1:], fnot_derivs[1:]])
    
    def get_s_standard_errors(self, params):
        # NOTE I don't ever use this function - this is the analytical result, but (maybe because the multinomial assumption is wrong) this outputs
        # unreasonably small errors
        assert len(params) == 2*self.L - 2
        f = self.deterministic_freqs(params)
        # calculating partial second derivatives
        Nf_one_minus_f = ((1-f).T*f.T*self.N)
        self.s_standard_errors = 1 / np.sqrt(np.sum(Nf_one_minus_f*(self.tps*self.tps), axis=1))

    def ll_step_deriv(self, step_size, *args):
        # the derivative of the log likelihood with respect to step size (given a search / slope vector)
        p, p_slopes = args
        s_slopes = np.concatenate([np.zeros(1), np.array(p_slopes[:self.L-1])])
        fnot_slopes = np.concatenate([np.zeros(1), np.array(p_slopes[self.L-1:])])
        params = p + step_size*p_slopes
        s = np.concatenate([np.zeros(1), np.array(params[:self.L-1])])
        fnot = np.clip(np.concatenate([np.zeros(1), np.array(params[self.L-1:])]), -100, 100)
        f_unscaled = np.exp(fnot + s*np.tile(self.tps, (self.L, 1)).T)
        f_unscaled_sums = np.sum(f_unscaled, axis=1)
        f = (f_unscaled.T/f_unscaled_sums).T
        slope_coefficient = fnot_slopes + s_slopes*np.tile(self.tps, (self.L, 1)).T
        Nf = (f.T*self.N).T
        ll_deriv = np.nansum((slope_coefficient * (self.R - Nf)), axis=(0, 1))
        return ll_deriv

    def ll_step(self, step_size, *args):
        # a function to minimize to find the optimal step size for gradient ascent (not using root finding)
        p, p_slopes = args
        tmp_params = p + step_size*p_slopes
        return self.ll_simple(tmp_params)

    def run_ll_max(self, param_guesses='NA', max_iters=1000, s_converge_thresh=1e-12, use_root=True, max_ss=1e-6, recording=True, neutral_assumption=False):
        if param_guesses != 'NA':
            self.current_params = param_guesses
        if neutral_assumption:
            self.current_params = np.array([0, self.current_params[1]])
        for i in range(max_iters):
            self.iterations += 1
            # finding slope
            self.current_slopes = self.ll_derivative(self.current_params, all_s_neutral=neutral_assumption)
            # finding optimal step size
            if use_root:
                tmss = max_ss  # ensuring that the max_ss is beyond the root
                while (self.ll_step_deriv(tmss, self.current_params, self.current_slopes) > 0):
                    tmss = tmss*10
                ss = optimize.brentq(self.ll_step_deriv, 0, tmss, args=(self.current_params, self.current_slopes))
            else:
                ss = optimize.minimize_scalar(self.ll_step, args=(self.current_params, self.current_slopes))['x']
            param_diff = self.current_slopes*ss
            self.current_params = self.current_params + param_diff
            self.current_ll = self.ll_simple(self.current_params)
            if self.record_process:
                self.ss_rec.append(ss)
                self.ll_rec.append(self.ll_simple(self.current_params))
                self.s_diffs.append([np.mean(np.abs(param_diff[:self.L-1])), np.max(param_diff[:self.L-1])])
                self.fnot_diffs.append([np.mean(np.abs(param_diff[self.L-1:])), np.max(param_diff[self.L-1:])])
            if np.max(np.abs(param_diff[:self.L-1])) < s_converge_thresh:
                self.converged = self.iterations
                break