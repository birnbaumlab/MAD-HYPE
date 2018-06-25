
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
    solvers = ['madhype']
    solver_options = [{}]

    # system parameters
    repeats = 50 #50
    cpw_range = np.logspace(0,4,41,dtype=int) # 41 breaks
    num_wells_range = [12,24,36,48]
    #num_wells_range = [60,72,84,96]

    # simulation parameters
    options = {
            'num_cells':1000,
            'num_wells':(24,),
            'seed':1,
            }

    id_map = np.zeros((len(cpw_range),options['num_cells']))

    for ii,num_wells in enumerate(num_wells_range): 
        for i,cpw in enumerate(cpw_range):
            for seed in xrange(repeats):

                specific_options = copy.copy(options)

                specific_options['cpw'] = (cpw,) 
                specific_options['num_wells'] = (num_wells,) 
                specific_options['seed'] = seed 

                _,results = simulate_run(solvers,solver_options,**specific_options)

                id_map[i,:] += results[0]['pattern']

            print 'Finished cpw = {}!'.format(cpw)

        # display graph
        c_labels = [v for v in [1,10,100,1000,10000] if v <= max(cpw_range)]
        c_inds = [min(range(len(cpw_range)), key=lambda i: abs(cpw_range[i]-v)) 
                    for v in c_labels]

        fig, ax = plt.subplots()

        cax = ax.imshow(id_map,interpolation='nearest')
        ax.set_aspect(aspect='auto')

        ax.set_yticks(c_inds)
        ax.set_yticklabels(c_labels)

        plt.title('Identification of clones with {} wells'.format(num_wells))
        plt.xlabel('Clonal Index')
        plt.ylabel('Cells per well (#)')
        cbar = fig.colorbar(cax, ticks=[0,repeats])
        cbar.ax.set_yticklabels(['0%','100%'])  # vertically oriented colorbar
        
        fig.savefig('{}_wells.eps'.format(num_wells),format='eps',dpi=1000)

    # show, allow user to exit
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()

