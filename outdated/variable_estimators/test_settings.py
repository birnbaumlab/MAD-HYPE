
"""
Attempt to test variable cells per well
"""

# standard libraries

# nonstandard libraries
import numpy as np

# homegrown libraries
from methods import *

def test_settings(*args,**kwargs):

    # default settings parameters
    settings = {
                'distribution':'power-law',
                'alpha':1.0,
                'clone_count':100000,
                'cpw_distro':[(1000,96)],
                'pattern_repeats':10
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)
    
    # change cpw_distro into dictionary if not there already
    if all([isinstance(s,tuple) for s in settings['cpw_distro']]):
        settings['cpw_distro'] = [{'cpw':a,'well_count':b} for a,b in settings['cpw_distro']]
    
    # assertion check
    assert all([isinstance(s,dict) for s in settings['cpw_distro']]),"cpw_distro incorrect type"

    # frequencies tested on plot
    plot_freqs = np.logspace(-1,-6,21)


    """ Runs tests on variable cells/well """
    
    # build patterns for specific frequencies of clones (for plotting)
    patterns = generate_simulated_patterns(plot_freqs,settings['cpw_distro'],settings['pattern_repeats'])

    # get repertoire frequencies
    freqs = generate_power_law_clones(settings['alpha'],settings['clone_count'])

    # get data for all the available patterns 
    plot_results = get_batch_pattern_probabilities(patterns,settings['cpw_distro'],freqs)

    # visualized create vectors
    visualize_results(plot_freqs,plot_results)
    

    #get_probability_pattern_match([10,10,10],settings['cpw_distro'],freqs)


if __name__ == '__main__':
    test_settings()
