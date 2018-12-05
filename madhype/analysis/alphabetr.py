

'''
Solver for ALPHABETR method
'''


#------------------------------------------------------------------------------# 

# standard libraries
import random
import sys
import itertools as it
import multiprocessing as mp
from multiprocessing import cpu_count

# nonstandard libraries
import numpy as np
import scipy.optimize, scipy.misc, scipy.cluster
from tqdm import tqdm

# homegrown libraries
from ..defaults import alphabetr_options as default_options

#------------------------------------------------------------------------------# 


def extract_chains(seq_data):
    alphas_per_well, betas_per_well = seq_data['well_data']['A'],seq_data['well_data']['B']
    return (
        list(set.union(*map(set, alphas_per_well))),
        list(set.union(*map(set, betas_per_well)))
    )


def solve(seq_data,**kwargs):

    """ Solver using ALPHABETR method """

    print 'ALPHABETR solver'

    # choose options
    options = default_options.copy()
    options.update(kwargs)

    # transfer dictionary to local namespace
    iters = options['iters']
    pair_threshold = options['pair_threshold']
    silent = options['silent']

    # set seed, if specified
    if 'seed' in options:
        random.seed(options['seed'])
   
    # Extract all distinct alpha- and beta-chains observed
    # TODO: might be better to extract the chains directly from the cells in the system
    all_alphas, all_betas = extract_chains(seq_data)
    
    # Create dictionary to look up alpha/beta index in constant time
    label_to_idx = {'A':{a: i for i,a in enumerate(all_alphas)},'B':{b: i for i,b in enumerate(all_betas)}}

    # Transform all well data to reference alpha- and beta-chains by index
    #well_data = [[[alpha_to_idx[a] for a in data[0]], [beta_to_idx[b] for b in data[1]]] for data in seq_data['well_data']]
    well_data = zip(*[[[label_to_idx[chain][l] for l in well_labels] 
            for well_labels in seq_data['well_data'][chain]] for chain in ('A','B')])

    if options['num_cores'] is None:
        overall_pairing_counts = _solve_singleprocessing(well_data, all_alphas, all_betas, **options)
    else:
        overall_pairing_counts = _solve_multiprocessing(well_data, all_alphas, all_betas, **options)

    overall_good_pairs = [pair for pair in overall_pairing_counts if overall_pairing_counts[pair]>=pair_threshold*iters]

    pairs = [(all_alphas[a], all_betas[b]) for a,b in overall_good_pairs]

    if options['dual_clones']:
        # Turns pairs of associated alpha- and beta- chains into cells that may have dual alpha chains
        cells, cell_freqs, cell_freqs_CI = pairs_to_cells(seq_data, pairs) 
        thresholds = [overall_pairing_counts[p]/float(iters) if len(p[0])==1 and len(p[1])==1 else None for p in cells]
    else:
        # Assume no dual clones
        cells = [((a,),(b,)) for a,b in pairs]
        thresholds = [overall_pairing_counts[p]/float(iters) for p in overall_good_pairs]
        if options['estimate_frequencies']:
            cell_freqs, cell_freqs_CI = estimate_cell_frequencies(seq_data, cells)
        else:
            cell_freqs, cell_freqs_CI = [0.0]*len(cells), [0.0]*len(cells)

    # NOTE: not using confidence intervals atm
    results = [(c,t,{'i':0.,'j':0.,'ij':f}) for c,t,f in zip(cells,thresholds,cell_freqs)]

#    print 'HERE:'
    for r in sorted(results,key=lambda x: -x[1]):
            pass#print r
    

    '''
    # OUTDATED
    results = {
        'cells': cells,
        'cell_frequencies': cell_freqs,
        'cell_frequencies_CI': cell_freqs_CI,
        'cell_thresholds': thresholds
    }
    '''

    return results

## Auxiliary solve functions
def _solve_singleprocessing(well_data, all_alphas, all_betas, **kwargs):
    solve_args = (well_data, len(all_alphas), len(all_betas), kwargs)
    iters = kwargs['iters']

    overall_pairing_counts = {}
    for iter in tqdm(range(iters), desc="Chain pairs", bar_format="{desc}: |{bar}| {percentage:3.0f}% [{elapsed}<{remaining}, {rate_fmt}{postfix}]"):
        pairs = _solve_iter(solve_args + (random.random(),))
        for p in pairs:
            overall_pairing_counts[p] = overall_pairing_counts.get(p, 0) + 1
    return overall_pairing_counts
def _solve_multiprocessing(well_data, all_alphas, all_betas, **kwargs):
    iters = kwargs['iters']
    num_cores = kwargs['num_cores']

    if num_cores == 0:
        num_cores = cpu_count()

    pool = mp.Pool(num_cores)
    pool_args = [(well_data, len(all_alphas), len(all_betas), kwargs, random.random()) for _ in range(iters)]
    results_iter = pool.imap_unordered(_solve_iter, pool_args)
    pool.close()
    
    overall_pairing_counts = {}
    for pairs in tqdm(results_iter, total=iters, desc="Chain pairs", bar_format="{desc}: |{bar}| {percentage:3.0f}% [{elapsed}<{remaining}, {rate_fmt}{postfix}]"):
        for p in pairs:
            overall_pairing_counts[p] = overall_pairing_counts.get(p, 0) + 1
    return overall_pairing_counts
def _solve_iter(args):
    well_data, num_alphas, num_betas, options, subprocess_seed = args

    # Choose random subset of wells for this iter
    # Each random subset is constant fraction (0.75) of all wells
    random.seed(subprocess_seed)
    wells_idx = random.sample(range(len(well_data)), int(0.75*len(well_data)))
    
    # Calculate association scores
    S = [[0 for j in range(num_betas)] for i in range(num_alphas)]
    for well_idx in wells_idx:
        well_alpha_idx, well_beta_idx = well_data[well_idx]
        for a_idx in well_alpha_idx:
            for b_idx in well_beta_idx:
                increment = 1./len(well_alpha_idx) + 1./len(well_beta_idx)
                S[a_idx][b_idx] += increment


    # Compute well pairings for every well
    # Then accumulate the number of times each pair has been assigned in a well pairing
    pairing_counts = {}
    for idx, well_idx in enumerate(wells_idx):
        well_pairings = _compute_well_pairings(*well_data[well_idx], scores=S)
        for a,b in well_pairings:
            pairing_counts[(a,b)] = pairing_counts.get((a,b), 0) + 1

    # Compute filter cutoff (average of all nonzero pair counts)
    cutoff = np.mean(pairing_counts.values())

    # Extract all pairs with counts exceeding the cutoff
    good_pairs = [pair for pair in pairing_counts if pairing_counts[pair]>cutoff]

    return good_pairs

## Computes a solution to the alpha-beta pairing problem, using the methods in Lee et al. (2017)
def _compute_well_pairings(alpha_idx, beta_idx, scores):

    if len(alpha_idx)==0 or len(beta_idx)==0:    
        return [] # scipy hungarian implementation doesn't handle this edge case

    # Reformulate problem as a general assignment problem
    # Then apply Hungarian algorithm
    # Indices in the results of running Hungarian are transformed back into alpha/beta chain ids
    ratings = [[-scores[i][j] for j in beta_idx] for i in alpha_idx]
    pairings = [(alpha_idx[i], beta_idx[j]) for i,j in zip(*scipy.optimize.linear_sum_assignment(ratings))]

    return pairings
 


def estimate_cell_frequencies(seq_data, cells):

    def log_likelihood_func(f, N, W, K, Q_memo = {}, error_rate=0.15, is_dual=False):
        # Note: See Eqs (3) and (4) in Lee et al. for explanation of variables
    
        # Compute Q vector if not previously computed
        Q_key = (tuple(N), tuple(W), error_rate, f, is_dual)
        if Q_key not in Q_memo:
            Q = []
            if not is_dual:
                prefactor = lambda m: 2*error_rate**m - error_rate**(2*m)
            else:
                prefactor = lambda m: 3*error_rate**m - 3*error_rate**(2*m) + error_rate**(3*m)
            for n,w in zip(N,W):
                q = (1-f)**n + sum([
                        prefactor(m) * scipy.misc.comb(n,m) * (f**m) * (1-f)**(n-m) 
                for m in range(1, n+1)])
                Q.append(q)
            Q_memo[Q_key] = Q
    
        # Retrieve Q from memoized dict of Q's
        # Note that Q only depends on N, W, error_rate, f, and is_dual
        Q = Q_memo[Q_key]
    
        # Compute log likelihood as sum of Binomial probabilities
        # Note that the "combinations" in the binomial PDF is ignored as it does not affect
        # the location of the maximum
        return sum([np.log(scipy.misc.comb(w,k)) + k*np.log((1-q)) + (w-k)*np.log(q) for w,k,q in zip(W,K,Q)])

    cells_per_well, N, W = extract_cells_per_well(seq_data)

    K = extract_cell_counts(seq_data, cells, cells_per_well, N, W)

    cell_freqs = []
    cell_freq_CIs = []

    with tqdm(total = len(cells), desc="Frequencies", bar_format="{desc}: |{bar}| {percentage:3.0f}% [{elapsed}<{remaining}, {rate_fmt}{postfix}]", unit='cell') as pbar:
        for (alist, blist), k in zip(cells, K):
            L_func = lambda f: log_likelihood_func(f, N, W, k, is_dual=len(alist)>1)
    
            # Find maximal likelihood
            f_opt = scipy.optimize.minimize_scalar(lambda f: -L_func(f), method='Bounded', bounds=(0,1)).x
            L_max = L_func(f_opt)
    
            # Find confidence interval, as specified in the paper
            f_min = scipy.optimize.minimize_scalar(lambda f: (L_max-1.96-L_func(f))**2, method='Bounded', bounds=(0,f_opt)).x
            f_max = scipy.optimize.minimize_scalar(lambda f: (L_max-1.96-L_func(f))**2, method='Bounded', bounds=(f_opt,1)).x
    
            cell_freqs.append(f_opt)
            cell_freq_CIs.append((f_min, f_max))
    
            pbar.update()

    return cell_freqs, cell_freq_CIs

def pairs_to_cells(seq_data, pairs):
    def find_duals_likelihood(candidate_duals, freqs_dict, well_size_cutoff = 50, error_rate=0.15):
        cells_per_well, N, W = extract_cells_per_well(seq_data)

        duals = []
        for alist, blist in candidate_duals:
            cells_temp = [((alist[0],), blist), ((alist[1],), blist), (alist, ()), (alist, blist)]
            #print "", cells_temp
            K = extract_cell_counts(seq_data, cells_temp, cells_per_well, N, W)
            #print "", K

            # Extract individual cell counts (see Lee et al., SI, Section 5 for explanation of variables)
            K_1 = K[0]
            K_2 = K[1]
            K_3 = K[2]
            K_d = K[3]
            K_o = [w-k1-k2-k3-kd for w,k1,k2,k3,kd in zip(W, K_1, K_2, K_3, K_d)]

            # Extract relevant cell frequencies
            f_q = freqs_dict[cells_temp[0]]
            f_r = freqs_dict[cells_temp[1]]
            f_d = freqs_dict[cells_temp[3]]

            # Null hypothesis (no dual clone)
            log_fact = lambda x: scipy.special.gammaln(x+1)

            # disgusting calculations
            null_P_a1b = [
                sum([
                    scipy.misc.comb(n, k) * f_q**k * (1-f_q-f_r)**(n-k) * (1-error_rate**k)**2
                    for k in range(1, n+1)
                ]) +
                sum([
                    np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2) + n_1*np.log(f_q) + n_2*np.log(f_r) + (n-n_1-n_2)*np.log(1-f_q-f_r) + 2*np.log(1-error_rate**n_1) + n_2*np.log(error_rate))
                    for n_1 in range(1, n) for n_2 in range(1, n-n_1+1)
                ]) +
                sum([
                    np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2) + n_1*np.log(f_q) + n_2*np.log(f_r) + (n-n_1-n_2)*np.log(1-f_q-f_r) + np.log(1-error_rate**n_1) + n_1*np.log(error_rate) + np.log(1-error_rate**n_2) + n_2*np.log(error_rate))
                    for n_1 in range(1, n) for n_2 in range(1, n-n_1+1)
                ])
                for n in N if n<well_size_cutoff
            ]
            null_P_a2b = [
                sum([
                    scipy.misc.comb(n, k) * f_r**k * (1-f_q-f_r)**(n-k) * (1-error_rate**k)**2
                    for k in range(1, n+1)
                ]) +
                sum([
                    np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2) + n_1*np.log(f_r) + n_2*np.log(f_q) + (n-n_1-n_2)*np.log(1-f_q-f_r) + 2*np.log(1-error_rate**n_2) + n_1*np.log(error_rate))
                    for n_2 in range(1, n) for n_1 in range(1, n-n_2+1)
                ]) +
                sum([
                    np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2) + n_1*np.log(f_r) + n_2*np.log(f_q) + (n-n_1-n_2)*np.log(1-f_q-f_r) + np.log(1-error_rate**n_1) + n_1*np.log(error_rate) + np.log(1-error_rate**n_2) + n_2*np.log(error_rate))
                    for n_2 in range(1, n) for n_1 in range(1, n-n_2+1)
                ])
                for n in N if n<well_size_cutoff
            ]
            null_P_a1a2 = [
                sum([
                    np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2) + n_1*np.log(f_q) + n_2*np.log(f_r) + (n-n_1-n_2)*np.log(1-f_q-f_r) + n_1*np.log(error_rate) + np.log(1-error_rate**n_1) + n_2*np.log(error_rate) + np.log(1-error_rate**n_2))
                    for n_1 in range(1,n) for n_2 in range(1, n-n_1+1)
                ])
                for n in N if n<well_size_cutoff
            ]
            null_P_a1a2b = [
                sum([
                    np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2)) * f_q**n_1 * f_r**n_2 * (1-f_q-f_r)**(n-n_1-n_2) * (
                        error_rate**n_1*(1-error_rate**n_1)*(1-error_rate**n_2)**2 +
                        (1-error_rate**n_1)**2*(1-error_rate**n_2)**2 +
                        (1-error_rate**n_1)**2*error_rate**n_2*(1-error_rate**n_2)
                    )
                    for n_1 in range(1,n) for n_2 in range(1, n-n_1+1)
                ])
                for n in N if n<well_size_cutoff
            ]
            null_P_other = [1-p1-p2-p3-pd for p1,p2,p3,pd in zip(null_P_a1b, null_P_a2b, null_P_a1a2, null_P_a1a2b)]
                
            null_log_likelihood = sum([
                log_fact(w) - log_fact(k1) - log_fact(k2) - log_fact(k3) - log_fact(kd) - log_fact(ko) + k1*np.log(p1) + k2*np.log(p2) + k3*np.log(p3) + kd*np.log(pd) + ko*np.log(po)
                for w,k1,k2,k3,kd,ko,p1,p2,p3,pd,po in zip(W, K_1, K_2, K_3, K_d, K_o, null_P_a1b, null_P_a2b, null_P_a1a2, null_P_a1a2b, null_P_other)
            ])

            # Alternative hypothesis (dual clone a1a2b)
            alt_P_2 = [
                sum([
                    scipy.misc.comb(n, k) * f_d**k * (1-f_d)**(n-k) * error_rate**k * (1-error_rate**k)**2
                    for k in range(1, n+1)
                ])
                for n in N if n<well_size_cutoff
            ]
            alt_P_3 = [
                sum([
                    scipy.misc.comb(n,k) * f_d**k * (1-f_d)**(n-k) * (1-error_rate**k)**3
                    for k in range(1, n+1)
                ])
                for n in N if n<well_size_cutoff
            ]
            alt_P_other = [1 - 3*p1 - p2 for p1,p2 in zip(alt_P_2, alt_P_3)]

            alt_log_likelihood = sum([
                log_fact(w) - log_fact(k1) - log_fact(k2) - log_fact(k3) - log_fact(kd) - log_fact(ko) + (k1+k2+k3)*np.log(p1) + kd*np.log(p2) + ko*np.log(po)
                for w,k1,k2,k3,kd,ko,p1,p2,po in zip(W, K_1, K_2, K_3, K_d, K_o, alt_P_2, alt_P_3, alt_P_other)
            ])

            if alt_log_likelihood - null_log_likelihood >= 10:
                duals.append(cells_temp[3])
            ## It appears Lee et al. include an additional restriction (see https://github.com/edwardslee/alphabetr/blob/f8fe7cdc1630e89ba173d0652b59b03724b1cade/R/dual_top.R#L160), that the null_log_likelihood be between 40 and 100.
            ## Not implemented here, as it's not part of the published methodology.

        return duals

    def find_duals_clustering(candidate_duals, freqs_dict):
        if len(candidate_duals) < 2:
            return []

        # Preliminaries
        cells_per_well, N, W = extract_cells_per_well(seq_data)
        K = extract_cell_counts(seq_data, candidate_duals, cells_per_well, N, W)

        # Compute ratio statistics
        R = []
        for (alist, blist), K_row in zip(candidate_duals, K):
            f1 = freqs_dict[((alist[0],), blist)]
            f2 = freqs_dict[((alist[1],), blist)]
            expected = sum([
                w*(1 - (1-f1)**n - (1-f2)**n + (1-f1-f2)**n)
                for n,w in zip(N,W)
            ])
            #print (alist, blist), f1, f2, sum(K_row), expected
            R.append(float(sum(K_row))/expected)

        # Perform clustering based on R
        centroids = sorted(scipy.cluster.vq.kmeans(R, 2)[0])
        if len(centroids)==1:    return []

        C_nondual, C_dual = centroids
        
        # Filter cells based on how close R is to C_dual vs. C_nondual
        duals = []
        #print centroids
        for candidate, r in zip(candidate_duals, R):
            if np.abs(r-C_dual) < np.abs(r-C_nondual):
                duals.append(candidate)
        
        return duals
        

    candidate_non_duals = [((a,),(b,)) for a,b in pairs] # list of all cells, will be modified as duals are found
    cells = candidate_non_duals[:]
    #print pairs

    # Determine candidate dual cells
    # Note that only dual alpha-chains are considered by the Lee et al. technique
    candidate_duals = []
    beta_pairings = {}
    for a,b in pairs:
        beta_pairings[b] = beta_pairings.get(b, []) + [a]
    for b,alist in beta_pairings.iteritems():
        if len(alist) >= 2:
            alist = sorted(alist)
            candidate_duals.extend(it.product(it.combinations(alist, 2), [(b,)]))

    #print beta_pairings

    # Estimate frequencies of all potential cells
    freqs_list, freqs_CI_list = estimate_cell_frequencies(seq_data, candidate_non_duals + candidate_duals)
    freqs_dict = {c: f for c,f in zip(candidate_non_duals+candidate_duals, freqs_list)}
    freqs_CI_dict = {c: f for c,f in zip(candidate_non_duals+candidate_duals, freqs_CI_list)}
    
    # Find duals using likelihood method, which is computationally infeasible with wells with >50 cells
    likelihood_duals = find_duals_likelihood(candidate_duals, freqs_dict)
    #print "Likelihood duals", likelihood_duals

    # Find duals using clustering method, which works better lower-frequency cells
    clustering_duals = find_duals_clustering(candidate_duals, freqs_dict)
    #print "Clustering duals", clustering_duals

    # Remove non-dual counterparts for each dual cell found and add in corresponding dual cell
    duals = list(set(likelihood_duals + clustering_duals))
    for alist,blist in duals:
        if ((alist[0],), blist) in cells:
            cells.remove(((alist[0],), blist))
        if ((alist[1],), blist) in cells:
            cells.remove(((alist[1],), blist))
        cells.append((alist, blist))
    
    cell_freqs = [freqs_dict[c] for c in cells]
    cell_freqs_CI = [freqs_CI_dict[c] for c in cells]

    return cells, cell_freqs, cell_freqs_CI
    
## Auxiliary functions to help out pairs_to_cells() and estimate_cell_frequencies()
def extract_cells_per_well(seq_data):
    cpw = seq_data['options']['cpw']
    num_wells = seq_data['options']['num_wells']
    cells_per_well = [c for c,w in zip(cpw,num_wells) for _ in xrange(w)]
    '''
    # OUTDATED
    if cpw_distro == 'constant':
        cells_per_well = [cpw_params['cells_per_well']]*len(seq_data.well_data)
    elif cpw_distro == 'poisson':
        cells_per_well = [cpw_params['lam']]*len(seq_data.well_data) # This is an approx. Not sure if it'll affect results
    elif cpw_distro == 'explicit':
        cells_per_well = cpw_params['cells_per_well']
    else:
        print "Unknown cell/well distribution: {0} with parameters {1}".format(cpw_distro, cpw_params)
        return None, None, None
    '''

    # Gather distinct #s of cells per well (N) and count number of wells w/ each cell count
    N_dict = {}
    for cpw in cells_per_well:
        N_dict[cpw] = N_dict.get(cpw, 0) + 1
    N,W = zip(*sorted(N_dict.iteritems()))
    return cells_per_well, N, W

def extract_cell_counts(seq_data, cells, cells_per_well, N, W):
    K = [[0]*len(N) for i in range(len(cells))]

    for well_size,well_data_a,well_data_b in zip(cells_per_well,seq_data['well_data']['A'],seq_data['well_data']['B']):
        well_alphas = set(well_data_a)
        well_betas = set(well_data_b)
        N_idx = N.index(well_size)
        for i,(alist,blist) in enumerate(cells):
            if all([a in well_alphas for a in alist]) and all([b in well_betas for b in blist]):
                K[i][N_idx] += 1

    return K

