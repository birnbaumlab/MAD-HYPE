
"""
Makes the set of preliminary figures
"""

# standard libraries

# nonstandard libraries
import numpy as np
import copy
import matplotlib
matplotlib.use('Agg')
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

    fig1,ax1 = plt.subplots(3,1,figsize=(8,12))
    fig2,ax2 = plt.subplots(1,3,figsize=(15,4))

    fig = plt.gcf()
    plt.subplots_adjust(wspace=0.45,bottom=0.2)
    #plt.subplots_adjust(left=0.1,right=0.9,wspace=0.25, hspace=0.2)

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
            'figsize':  (12,9),
            'xlim':((2.)*(10.**-4),(6./5)*(10.**-2)),
            'ylim':((2.)*(10.**-4),(6./5)*(10.**-2)),
            },
        'num_cells':1000,
        'cell_freq_max':         0.01, # 0.01
        'cell_freq_constant':       2,
        'visual':                True,
        'silent':               False,
        }


    settings_list = [
            {
                'num_wells':(42,54),
                'cpw':(25,1758),
            },
            {
                'num_wells':(60,36),
                'cpw':(159,2402),
            },
            {
                'num_wells':(96,),
                'cpw':(1000,),
            },
            ]


    for i,sub_settings in enumerate(settings_list):

        specific_options = copy.copy(base_settings)

        specific_options.update(sub_settings)
        specific_options['plot_repertoire']['ax'] = ax1[i]
        specific_options['plot_frequency_estimation']['ax'] = ax2[i]

        _,results = madhype.simulate_run(solvers, solver_options, **specific_options)

    fig1.savefig('Figure 5B.png', format='png', dpi=300)
    fig2.savefig('Figure 5C.png', format='png', dpi=300)

#    plt.show(block=False)
#    raw_input('Press enter to close...')
#    plt.close()
    

if __name__ == '__main__':
    main()
    
    




