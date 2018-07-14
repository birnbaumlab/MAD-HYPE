

'''
Solver for ALPHABETR method
using the R package published by Lee et al. (2017)
'''


#------------------------------------------------------------------------------# 

# standard libraries
import random
import inspect
import os,sys
import itertools as it
import multiprocessing as mp
from multiprocessing import cpu_count
import tempfile
import subprocess
import shutil

# nonstandard libraries
import numpy as np
import scipy.optimize, scipy.misc, scipy.cluster
import matplotlib.pyplot as plt
from tqdm import tqdm

# homegrown libraries
from ..defaults import alphabetr_options as default_options

#------------------------------------------------------------------------------# 


def extract_chains(seq_data):
    alphas_per_well, betas_per_well = seq_data['well_data']['A'],seq_data['well_data']['B']
    return list(set.union(*alphas_per_well)),list(set.union(*betas_per_well))


def solve(seq_data,**kwargs):

    """ Solver using ALPHABETR method """

    print 'ALPHABETR solver (R)...'

    # choose options
    options = default_options.copy()
    options.update(kwargs)

    # transfer dictionary to local namespace
    iters = options['iters']
    pair_threshold = options['pair_threshold']
    silent = options['silent']

    pairs, pair_scores = run_alphabetr_R(seq_data, iters=iters, pair_threshold = pair_threshold)

    # note: dual clones and frequency estimates not implemented
    cells = [((a,),(b,)) for a,b in pairs]

    # NOTE: not using confidence intervals atm
    results = [(c,s,{'i':0.,'j':0.,'ij':0.}) for c,s in zip(cells,pair_scores)]

    return results

## Auxiliary solve functions
def run_alphabetr_R(seq_data, iters, pair_threshold):
    # Extract all distinct alpha- and beta-chains observed
    all_alphas, all_betas = extract_chains(seq_data)

    # Transform all well data to reference alpha- and beta-chains by index
    label_to_idx = {'A':{a: i for i,a in enumerate(all_alphas)},'B':{b: i for i,b in enumerate(all_betas)}}
    well_data = zip(*[[[label_to_idx[chain][l] for l in well_labels] 
            for well_labels in seq_data['well_data'][chain]] for chain in ('A','B')])

    # Create temporary folder and relevant datafile paths
    tempdir = tempfile.mkdtemp(dir = os.getcwd())
    alpha_path = os.path.join(tempdir, 'alpha_welldata.csv')
    beta_path = os.path.join(tempdir, 'beta_welldata.csv')
    R_script_path = os.path.join(tempdir, 'alphabetr.R')
    output_path = os.path.join(tempdir, 'output.csv')

    # Export well_data into format for R implementation to read
    num_wells = sum(seq_data['options']['num_wells'])
    num_alphas = len(all_alphas)
    num_betas = len(all_betas)
    export_well_data_alphabetr_R(well_data, num_alphas, num_betas, alpha_path, beta_path)

    # Call R implementation via command line
    make_R_script(R_script_path)
    subprocess.call(['rscript', R_script_path, str(iters), str(pair_threshold), alpha_path, beta_path, output_path])

    # Parse output file (all alpha/beta specified are indices into all_alphas/all_betas)
    output = import_output_alphabetr_R(output_path)

    # Clean up temporary files
    shutil.rmtree(tempdir)

    # Return two lists, a list of pairs and a list of scores for each pair
    pairs = [(all_alphas[a], all_betas[b]) for (a,b),_ in output]
    pair_scores = [r[1] for r in output]

    return pairs, pair_scores

def make_R_script(path):
    src = """
library(alphabetr)

args = commandArgs(trailingOnly=TRUE)
reps = as.integer(args[1])
pair_threshold = as.numeric(args[2])
alpha_path = args[3]
beta_path = args[4]
output_path = args[5]

# Read well data files
adata=as.matrix(read.csv(alpha_path))
bdata=as.matrix(read.csv(beta_path))

# Run alphabetr
bagpipe_res = bagpipe(adata, bdata, replicates=reps)
bagpipe_filtered = bagpipe_res[bagpipe_res[, "prop_replicates"] >= pair_threshold,]

# Remove extra columns
bagpipe_filtered = bagpipe_filtered[,c(-2,-4)]

write.table(bagpipe_filtered, file=output_path, sep=",", row.names=FALSE, col.names=FALSE)
    """

    out = open(path, 'w')
    out.write(src)
    out.close()

def export_well_data_alphabetr_R(well_data, num_alphas, num_betas, alpha_path, beta_path):
    num_wells = len(well_data)

    out_matrix_alpha = np.zeros((num_wells, num_alphas), dtype=np.int)
    out_matrix_beta = np.zeros((num_wells, num_betas), dtype=np.int)
    for well_idx, (alist, blist) in enumerate(well_data):
        for a in alist:
            out_matrix_alpha[well_idx][a] = 1
        for b in blist:
            out_matrix_beta[well_idx][b] = 1

    out_alpha = open(alpha_path, 'w')
    out_alpha.write(
        '\n'.join(
            ','.join(str(v) for v in adata)
            for adata in out_matrix_alpha
        )
    )
    out_alpha.close()

    out_beta = open(beta_path, 'w')
    out_beta.write(
        '\n'.join(
            ','.join(str(v) for v in bdata)
            for bdata in out_matrix_beta
        )
    )
    out_beta.close()

def import_output_alphabetr_R(output_path):
    indata = open(output_path)
    lines = indata.readlines()
    indata.close()

    output = []
    for l in lines:
      bstr, astr, scorestr = l.split(",")
      output.append(((int(astr)-1, int(bstr)-1), float(scorestr)))

    return output

