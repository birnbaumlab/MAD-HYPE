
# Testing happens here

# standard libraries


# nonstandard libraries
from matplotlib import pyplot as plt
import numpy as np

# homegrown libraries

from methods.madhype import solve as solve_madhype
from datasim.seq_generator import *

from variable import * # functions: test_settings
from analysis import * # functions: process_results


#cpw_distro = [(10,20),(20,20),(30,20),(40,20),(50,20)]
#cpw_distro = [(1,20),(4,20),(16,20),(64,20),(256,20)]
#cpw_distro = [(2777,96)]
cpw_distro = [(12,24),(25,24),(50,24),(50,24)]

cpw_by_well = [wt[0] for wt in cpw_distro for i in xrange(wt[1])]

options = {}

sg = SequencingGenerator(**options)

sg.set_cells_per_well(distro_type = 'explicit',
                      cells_per_well = cpw_by_well)

# generate cell data

'''
Options: 
    > num_wells
    > cells_per_well_distibution
    > num_cells
    > cells
    > cell_frequency_distribution
    > chain_misplacement_prob
    > chain_deletion_prob
    > alpha_sharing_probs
    > beta_sharing_probs
'''

options = {
            'num_cells':100,
            'chain_deletion_prob':0.1,
            'alpha_sharing_probs':0.0,
            'beta_sharing_probs':0.0
          }


sg.set_options(**options)

data = sg.generate_data()

# generate results from solver, from data
results = solve_madhype(data,pair_threshold=2.0)

# assertion check

#
<<<<<<< HEAD
#results_processed = process_results(results,data)

=======
repertoire_coverage(results,data)
#results_processed = visualize_results(results,data)
>>>>>>> a18e16b57e8d76afed664be80425424d0e5f9075


#test_settings(cpw_distro=cpw_distro)


