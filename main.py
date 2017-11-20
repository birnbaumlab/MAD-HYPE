
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


cpw_distro = [(10,20),(100,20),(1000,20),(10000,20),(10000,20)]
#cpw_distro = [(2777,96)]

cpw_by_well = [wt[0] for wt in cpw_distro for i in xrange(wt[1])]

options = {}

sg = SequencingGenerator(**options)

sg.set_cells_per_well(distro_type = 'explicit',
                      cells_per_well = cpw_by_well)

# generate cell data
sg.generate_cells(100)
data = sg.generate_data()

# generate results from solver, from data
results = solve_madhype(data,pair_threshold=2.0)

# assertion check

#
#results_processed = process_results(results,data)



#test_settings(cpw_distro=cpw_distro)


