
# Testing happens here

# standard libraries


# nonstandard libraries
from matplotlib import pyplot as plt
import numpy as np

# homegrown libraries
from datasim.seq_generator import *
from variable import *
from methods.madhype import solve as solve_madhype


cpw_distro = [(10,24),(100,24),(1000,24),(10000,24)]
#cpw_distro = [(10000,96)]

cpw_by_well = [wt[0] for wt in cpw_distro for i in xrange(wt[1])]

options = {}

sg = SequencingGenerator(**options)

sg.set_cells_per_well(distro_type = 'explicit',
                      cells_per_well = cpw_by_well)

sg.generate_cells(100)
data = sg.generate_data()

results = solve_madhype(data,pair_threshold = 0.99,verbose=0,real_data=False)


#test_settings(cpw_distro=cpw_distro)


