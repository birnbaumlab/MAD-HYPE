import random
import itertools as it
import multiprocessing as mp

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "serif"


# Import MAD-HYPE package
import madhype

def run_iter(args):
    solvers, solver_options, general_options, seed = args

    for op in solver_options:  op['seed'] = seed

    times = []
    for s,s_ops in zip(solvers, solver_options):
      tstart = time.time()
      print "{} START TIME: {}".format(s, tstart)
      madhype.simulate_run([s], [s_ops], seed = seed, **general_options)
      tstop = time.time()
      times.append(tstop - tstart)

    return times

def collect_data():
    # Set up parameters that apply to all solvers/simulations
    general_options = {
            'num_cells': 1000,
    
            'cell_freq_constant': 2,
            'cell_freq_max': 0.05,
            'alpha_sharing_probs': None,
            'beta_sharing_probs': None,
    
            'fdr': 0.05,
    
            'visual': False,
    
            'silent': True,
            }
    
    w_range = [96]
    c_range=[5,10,20,30,50,100, 200, 300, 500, 1000, 2000]
    repeats = 15
    
    madhype_times = {(c,w): [] for c,w in it.product(c_range, w_range)}
    alphabetr_times = {(c,w): [] for c,w in it.product(c_range, w_range)}

    print "Running simulations with the following settings:"
    print "  Well partitions:", w_range
    print "  CPW:", c_range
    print "  # reps/condition", repeats
    
    for c in c_range:
      for w in w_range:
        num_wells = (w,)
        cpw = (c,)
        print "Number of wells: {}; cells per well: {}".format(num_wells, cpw)
    
        ops = general_options.copy()
        ops['num_wells'] = num_wells
        ops['cpw'] = cpw
    
        seeds = [random.randint(0,1e6) for _ in xrange(repeats)]
    
        # Get compute times
        args_list = [(['madhype', 'alphabetr_R'],[{'silent': True,'num_cores': None}, {'iters': 100}],ops,seed) for seed in seeds]
        pool = mp.Pool(mp.cpu_count())
        times = pool.map(run_iter, args_list)
        madhype_results[(c,w)], alphabetr_results[(c,w)] = zip(*times)
    
        times_arr = np.array([madhype_results[(c,w)], alphabetr_results[(c,w)]])
        print "num_wells = {}, cpw={}:".format(w, c)
        print "  MH times:", times_arr[0,:], np.mean(times_arr[0,:])
        print "  AB times:", times_arr[1,:], np.mean(times_arr[1,:])
    
    return madhype_results, alphabetr_results
    
def plot_data(madhype_results, alphabetr_results):
    solvers = ['MAD-HYPE', 'ALPHABETR']
    results = [madhype_results, alphabetr_results]
    
    colors = [(.24, .52, .78), (1.0, .7, .4)]
    
    plt.figure()
    for solver,res,c in zip(solvers, results, colors):
      keys = sorted(res.keys())
      x = [cpw for cpw,_ in keys]
    
      y_raw = np.array([res[k] for k in keys])
      y = np.mean(y_raw, axis = 1)
      yerr = np.std(y_raw, axis=1)
    
      plt.errorbar(x,y, yerr=yerr, capsize=3, color=c, linewidth=2)
    
      print solver, ":"
      print "  x:", x
      print "  y:", list(y)
      print "  yerr:", list(yerr)
    
    plt.legend(solvers)
    plt.xlabel('Cells per Well')
    plt.ylabel('Computation Time (sec)')
    plt.yscale('log')

    plt.savefig('figure_S8.pdf')

madhype_results, alphabetr_results = collect_data()
plot_data(madhype_results, alphabetr_results)
