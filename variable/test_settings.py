
"""
Attempt to test variable cells per well
"""

# standard libraries

# nonstandard libraries
import numpy as np

# homegrown libraries
from methods import *

def test_settings(*args,**kwargs):

    # create list of tuples for (cpw,well_counts)
    cpw_distro = [(10,24),(100,24),(1000,24),(10000,24)]
    #cpw_distro = [(10,32),(100,32),(1000,32)]
    cpw_distro = [(10000,96)]
    
    # default settings parameters
    settings = {
                'distribution':'power-law',
                'alpha':1.0,
                'clone_count':100000,
                'cpw_distro':[{'cpw':a,'well_count':b} for a,b in cpw_distro],
                'pattern_repeats':20
               }
    
    plot_freqs = np.logspace(-1,-6,41)

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

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
