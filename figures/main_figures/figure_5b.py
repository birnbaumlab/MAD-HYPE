import random
import itertools as it

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#plt.ion()

# Import MAD-HYPE package
import madhype

def collect_data():
    # Set up run parameters
    solvers = ['madhype', 'alphabetr']
    solver_options = [{'silent': True}, {'silent': True, 'estimate_frequencies': False}]
    
    # Set up parameters that apply to all solvers/simulations
    general_options = {
            'num_cells': 1000,
    
            'cell_freq_constant': 2,
            'cell_freq_max': 0.01,
            'alpha_sharing_probs': None,
            'beta_sharing_probs': None,
    
            'visual': False,
    
            'silent': True,
            }
    
    total_cells = 9600
    total_wells = 96
    
    w_range = np.arange(12,96,12)
    c_range = np.ceil(np.logspace(0,np.log10(total_cells/total_wells),11)).astype(np.int)
    repeats = 10
    
    madhype_results = {(c,w): [] for c,w in it.product(c_range, w_range)}
    alphabetr_results = {(c,w): [] for c,w in it.product(c_range, w_range)}
    
    print "Running simulations with the following settings:"
    print "  Well partitions:", w_range
    print "  CPW:", c_range
    print "  # reps/condition", repeats
    
    for i,c in enumerate(c_range):
      for j,w in enumerate(w_range):
        num_wells = (w, total_wells - w)
        cpw = (c, int((total_cells - c*w) / (total_wells - w)))
        print "Number of wells: {}, cells per well: {}".format(num_wells, cpw)
        for rep in xrange(repeats):
          seed = random.randint(0,1e6)
          print "Iteration {}: SEED {}".format(rep, seed)
    
          _, results = madhype.simulate_run(solvers, solver_options, num_wells = num_wells, cpw=cpw, seed=seed, **general_options)
    
          madhype_results[(c,w)].append(results[0]['frac_repertoire'])
          alphabetr_results[(c,w)].append(results[1]['frac_repertoire'])

    return madhype_results, alphabetr_results
    
def plot_data(madhype_results, alphabetr_results):
    vmin = 0.0
    vmax = 0.7

    solvers = ['MAD-HYPE', 'ALPHABETR']
    
    for solver, res in zip(solvers, [madhype_results, alphabetr_results]):
      cpw,num_wells = zip(*res.keys())
      cpw = sorted(set(cpw))
      num_wells = sorted(set(num_wells))

      matvals = np.array([[np.mean(res[(c,w)]) for w in num_wells] for c in cpw])

      plt.figure()
      plt.imshow(matvals, interpolation='nearest', vmin=vmin, vmax=vmax)
      plt.colorbar()
      plt.xticks(range(len(num_wells)), num_wells)
      plt.xlabel('# wells in first partition')
      plt.yticks(range(len(cpw)), cpw)
      plt.ylabel('# cells/well in first partition')
      plt.title('{}, 9600 cells'.format(solver))
      plt.savefig('Figure5B_{}.pdf'.format(solver))

madhype_results, alphabetr_results = collect_data()
plot_data(madhype_results, alphabetr_results)
