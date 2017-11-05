
"""

File: variable/methods.py

Collection of common methods used in MAD-HYPE

Included functions:
    generate_power_law_clones(a,n)
    get_probability_pattern_match(pattern,well_distro,freqs)
    generate_simulated_patterns(freqs,well_distro,repeats=1)
    get

"""

#----------------------------------------------------------------#

# standard libraries

# nonstandard libraries
import numpy as np
from matplotlib import pyplot as plt
from scipy.misc import comb
from scipy.stats import binom
from scipy.misc import comb

# homegrown libraries


# library modifications
plt.rcParams["font.family"] = "serif"


#----------------------------------------------------------------#
# Main Functions
#----------------------------------------------------------------#

def generate_power_law_clones(a,n):
    """ Generates n clones under a frequency distribution w/ alpha = a """ 
    freqs = 10.**(-a*np.log10(np.arange(1,n+1)))
    freqs = freqs/sum(freqs)
    return freqs

#----------------------------------------------------------------#

def get_batch_pattern_probabilities(patterns,well_distro,freqs):
    """ Checks entire batch of patterns for their intersection likelihood """
    # initialize avg,std results
    results_dict = {'avg':[],'std':[]}

    # assertion check
    assert isinstance(patterns,np.ndarray),"submitted patterns are not the right type"
    
    # iterate across each pattern set
    for i in xrange(patterns.shape[1]):
        # get clonal intersection probability
        clonal_probs = [get_probability_pattern_match(patterns[j,i,:],well_distro,freqs)
                for j in xrange(patterns.shape[0])]

        # store results_dict in expanding dictionary
        results_dict['avg'].append(np.mean(clonal_probs))
        results_dict['std'].append(np.std(clonal_probs))

    # return results_dict dict
    return results_dict 

#----------------------------------------------------------------#

def generate_simulated_patterns(freqs,well_distro,repeats=1):
    """ Simulates well distro patterns for frequencies, given well_distro and repeats"""
    """
    Outputs:
        patterns - np.array, of size (repeats,len(freqs),len(distros))
    """

    # assertion check
    assert isinstance(repeats,int),"repeats needs to be <type int>"

    # initialize pattern array
    patterns = np.zeros((repeats,len(freqs),len(well_distro)))

    # iterate across well distribution partitions
    for r in xrange(repeats):
        for i,distro in enumerate(well_distro):
            patterns[r,:,i] = np.random.binomial(distro['well_count'],_calculate_well_freqs(freqs,distro['cpw']))

    # return pattern array
    return patterns.astype(int)

#----------------------------------------------------------------#

def get_probability_pattern_match(pattern,well_distro,freqs):

    """

    Calculates a probability of a repertoire matching a particular distribution pattern

    Input:
        pattern - tuple or list of number of wells occupied, in order of well_distro
        well_distro - tuple or list of dictionaries, each with kwargs (cpw,well_count)
        freqs - tuple or list of clonal frequencies

    Output:
        prob - probability that there is atleast one clone that occupies same pattern

    """

    # assertion checks
    assert len(pattern) == len(well_distro), "input pattern not same length as well_distribution"
    assert all([isinstance(distro,dict) for distro in well_distro]), "not all entries in well_distro are dicts"
    assert all([isinstance(p,int) for p in pattern]), "not all entries in pattern are ints"

    # initialize variables
    prob_clonal_pattern = np.ones((len(freqs),))

    # iterate through well count partitions
    for p,distro in zip(pattern,well_distro):
        well_freqs = _calculate_well_freqs(freqs,distro['cpw'])
        subset_prob = binom.pmf(p,distro['well_count'],well_freqs)*(1/comb(distro['well_count'],p))
        prob_clonal_pattern = np.multiply(prob_clonal_pattern,subset_prob)

    # returns the product of probabilities
    return np.prod(1. - prob_clonal_pattern) 

#----------------------------------------------------------------#
# Visualization functions
#----------------------------------------------------------------#

def visualize_results(freqs,results_dict,*args,**kwargs):

    options = {
              'title':'My data'
              }
    
    # update graphing options with user-inputs
    for arg in args: options.update(arg)
    options.update(kwargs)

    # assertion check
    assert isinstance(results_dict,dict),"results dict is not dict"
    assert all([k in results_dict.keys() for k in ['avg','std']]),"avg and std not keys in dict"

    # just plot the results using errorbars
    plt.errorbar(freqs,results_dict['avg'],yerr=results_dict['std'])
    
    # making things pretty
    plt.xscale('log')
    plt.title(options['title'])
    plt.xlabel('Clonal frequency')
    plt.ylabel('% Successful Identification')
    plt.ylim([-0.1,1.1])


    # show plot, with user-input exit
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()


#----------------------------------------------------------------#
# Internal functions
#----------------------------------------------------------------#


def _calculate_well_freqs(freqs,cpw):
    """ Calculates likelihood of observing clones given clonal freqs and cpw """
    return 1. - ((1. - freqs)**cpw)

#----------------------------------------------------------------#




#----------------------------------------------------------------#




#----------------------------------------------------------------#




#----------------------------------------------------------------#


