
"""
Makes the set of preliminary figures
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# homegrown libraries
from madhype.analysis.methods import match_probability
import madhype

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 16 


def main():

    settings = {
            
               }

    # update settings
    #for arg in args: settings.update(arg)
    #settings.update(kwargs)

    """ Match probability as a function of w_a,w_b """

    fig1,ax1 = plt.subplot(3,1,figsize=(12,10))
    fig2,ax2 = plt.subplot(3,1,figsize=(12,5))

    fig = plt.gcf()
    plt.subplots_adjust(wspace=0.45, hspace=0.3)

    # some default settings
    solvers = ['madhype']
    solver_options = [{}]

    base_settings = {
        'plot_repertoire': {
            'ax':ax1,
            'fig':fig1,
            },
        'plot_frequency_estimation': {
            'ax':ax2,
            'fig':fig2,
            },
        'visual':                True,
        'silent':               False,
        }


    settings_list = [
            {
                'num_wells':(42,54),
                'cpw':(25,1758),
            },
            {
                'num_wells':(48,48),
                'cpw':(250,1750),
            },
            {
                'num_wells':(96,),
                'cpw':(1000,),
            },
            ]


    for i,sub_settings in enumerate(settings_list):

         

        _,results = madhype.simulate_run(solvers, solver_options, **settings)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()
    

if __name__ == '__main__':
    main()
    

