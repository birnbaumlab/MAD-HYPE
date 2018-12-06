import random
import itertools as it

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Import MAD-HYPE package
import madhype

def collect_data():
    # Set up run parameters
    solvers = ['madhype', 'alphabetr']
    solver_options = [{'fdr': 0.05},{'estimate_frequencies': False}]
    
    # Set up parameters that apply to all solvers/simulations
    general_options = {
            'num_cells': 3000,
    
            'cell_freq_constant': 2,
            'cell_freq_max': 0.05,
            'alpha_sharing_probs': None,
            'beta_sharing_probs': None,
    
            'num_wells': (96,),
            'cpw': (300,),
    
            'visual': False,
    
            'fdr': 0.05,
            }
    
    repeats = 20
    
    madhype_results = []
    alphabetr_results = []

    for rep in xrange(repeats):
      seed = random.randint(0,1e6)
      print "Iteration {}: SEED {}".format(rep, seed)

      solver_options[1]['seed'] = seed
      _, results = madhype.simulate_run(solvers, solver_options, seed = seed, **ops)

      madhype_results.append(results[0])
      alphabetr_results.append(results[1])

def plot_data(madhype_results, alphabetr_results):
    solvers = ['MAD-HYPE', 'ALPHABETR']
    results = [madhype_results, alphabetr_results]

    colors = [(0.0,0.5,0.18), (1.0,0.85,0.4)]

    min_freq = min(min(res[0]['freqs']) for res in results)
    max_freq = max(max(res[0]['freqs']) for res in results)
   
    plt.figure()
    plt.xscale('log')
    for c,m in zip(colors,results):
      pmean = np.mean(np.array([r['pattern'] for r in m]), axis=0)
      freqs = m[0]['freqs']
   
      w=1.5
      x = freqs[1:]
      y = [np.mean([pmean[i] for i in range(len(pmean)) if f_min/w<freqs[i]<f_min*w]) for f_min in x]
      plt.plot(x, y, color=c, linewidth=4)
    
    plt.xlabel('Clone frequency')
    plt.ylabel('Fraction identified')
    plt.legend(solvers)
    
    plt.savefig('figure_S5.pdf')
    
    mpos = [r['positives'] for r in results[0]]
    mfr = [r['frac_repertoire'] for r in results[0]]
    apos = [r['positives'] for r in results[1]]
    afr = [r['frac_repertoire'] for r in results[1]]
    
    print "MAD-HYPE:"
    print "  # clones:", np.mean(mpos), "+/-", np.std(mpos)
    print "  % rep.:", np.mean(mfr), "+/-", np.std(mfr)
    print "ALPHABETR:"
    print "  # clones:", np.mean(apos), "+/-", np.std(apos)
    print "  % rep.:", np.mean(afr), "+/-", np.std(afr)
   
madhype_results, alphabetr_results = collect_data()
plot_data(madhype_results, alphabetr_results)
