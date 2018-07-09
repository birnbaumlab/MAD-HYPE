
"""
Supplemnetary Figure 1
"""

# nonstandard libraries
import matplotlib.pyplot as plt
from scipy.optimize import minimize,root
import numpy as np

# homegrown libraries
import madhype
from madhype.simulation import *  # we should explicitly import here
from madhype import simulate_run

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 16 

def main(*args,**kwargs):

    """ Loads data, does some pretty basic simulations """


    labels = ['Peripheral Blood (9-25 y)','Peripheral Blood (61-66 y)']

    specific_settings = {
            # ranges 0.1% to 3%
            'Peripheral Blood (9-25 y)':{
                'cell_freq_max':0.0025,
                'cell_freq_constant':1 + 1./1.2,
                'num_cells':10000,
                'threshold':2.0,
                'block':False,
                'num_wells':(48,48),
                'cpw':(500,10000)
                },
            'Peripheral Blood (61-66 y)':{
                'cell_freq_max':0.086,
                'cell_freq_constant':1 + 1./1.2,
                'num_cells':10000,
                'threshold':2.0,
                'block':True,
                'num_wells':(48,48),
                'cpw':(500,10000)
                },
            }

    # figure specific properties
    fig,axes = plt.subplots(nrows=1, ncols=len(labels), figsize=(15,4))
    #plt.subplots_adjust(left=0.15,right=0.9,top=0.85,bottom=0.3,hspace=0.5,wspace=0.5)

    # some default settings
    solvers = ['madhype']
    solver_options = [{}]

    for ind,label in enumerate(labels):

        settings = {
            'plot_repertoire': {
                'ax':axes[ind],
                'fig':fig,
                },
            'visual':                True,
            'silent':               False,
            'legend':           bool(ind),
            }


        specific_settings[label].update(settings)
        data,results = madhype.simulate_run(solvers, solver_options, **specific_settings[label])

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()

#------------------------------------------------------------------------------# 

if __name__ == '__main__':
    main()















