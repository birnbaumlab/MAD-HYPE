
"""
Tests out the new and improved variable solver
"""

# standard libraries
import time
import copy
from collections import Counter
from operator import mul
from sys import argv

# nonstandard libraries
from scipy.misc import comb
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# intrapackage libraries
from madhype import simulate_run
#from defaults import general_options as default_options

        

#------------------------------------------------------------------------------# 

""" Main Callable Method """

#------------------------------------------------------------------------------# 

if __name__ == "__main__":

    """ Use command line arguments """
    
    if len(argv) > 1:
        mode = int(argv[1])
        print 'User set mode to {}...'.format(mode)
        time.sleep(1.0)
    else:
        mode = 1
        print 'Default mode set to {}...'.format(mode)
        time.sleep(1.0)

    solvers = ['madhype']
    solver_options = [{}]

    ''' Generates conserved cell count experiment '''
    options = {
            'num_cells':1000,
            'cell_freq_max':0.01,
            'cell_freq_constant':2.,
            # visual cues
            'silent':True,
            'visual':False
            }

    T = 96000
    W = 96
    #w_range = [24,36,48]
    w_range = [0,6,12,18,24,30,36,42,48] # real well range
    c_range = np.logspace(0,3,16,dtype=int) # 16 segments
    repeats = 3 # 3 repeats

    id_map = np.zeros((len(c_range),len(w_range)))

    for i,c in enumerate(c_range):
        for j,w in enumerate(w_range):
            for seed in xrange(repeats):
                #print 'Starting cpw = {}, num_wells = {}'.format((c,int((T - c*w)/(W-w))),(w,W-w))
                specific_options = copy.copy(options)

                specific_options['seed'] = seed
                specific_options['cpw'] = (c,int((T - c*w)/(W-w)))
                specific_options['num_wells'] = (w,W-w)

                _,results = simulate_run(solvers,solver_options,**specific_options)

                id_map[i,j] += results[0]['frac_repertoire']

            print 'Finished w1 = {}!'.format(w)
        print 'Finished c1 = {}!'.format(c)

    # normalize fraction of repertoires
    id_map /= repeats

    c_labels = [v for v in [1,10,100] if v <= max(c_range)]
    c_inds = [min(range(len(c_range)), key=lambda i: abs(c_range[i]-v)) for v in c_labels]

    fig, ax = plt.subplots() # create blank figure
    cax = ax.imshow(id_map,interpolation='nearest',vmin=0,vmax=1.0) # fill with heatmap
    ax.set_aspect(aspect='auto') # fix aspect ratios

    plt.title('Fixed cell count = {}'.format(T))
    plt.xlabel('# wells in partition')
    plt.ylabel('# cells/well in partition')

    ax.set_xticks([i for i in xrange(len(w_range))])
    ax.set_xticklabels(w_range)

    ax.set_yticks(c_inds)
    ax.set_yticklabels(c_labels)

    # colorbar in action
    cbar = plt.colorbar(cax,ax=ax,ticks=[0,1])
    cbar.ax.set_yticklabels(['0%','100%'])  # vertically oriented colorbar

    plt.savefig('figure_4a.png', format='png', dpi=300)

#    # show, allow user to exit
#    plt.show(block=False)
#    raw_input('Press enter to close...')
#    plt.close()


