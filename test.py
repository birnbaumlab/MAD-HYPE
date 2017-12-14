
import numpy as np
from scipy.optimize import root

def _power_law_distribution(num_cells,max_freq,alpha):

    # Lower bound
    if max_freq <= 1./num_cells:
        print 'Max. freq too low! Returning uniform distribution...'
        return [1./num_cells for _ in xrange(num_cells)]
    
    # Upper bound
    if max_freq >= 1.:
        print 'Max. freq too high! Returning delta distribution...'
        return [1] + [0 for _ in xrange(num_cells-1)]
 
    # Find a shift
    shift = root(_get_max_freq_diff, 1.0, args=(num_cells,max_freq,alpha)).x
    
    # Find best
    return _get_freq_distribution(shift,num_cells,max_freq,alpha)

def _get_max_freq_diff(shift,num_cells,max_freq,alpha):
    """ Function for finding diff. b/w max_freq and current distribution """
    freqs = _get_freq_distribution(shift,num_cells,max_freq,alpha)
    return max_freq - freqs[0]

def _get_freq_distribution(shift,num_cells,max_freq,alpha):
    """ Generate a normalized power-law distribution """
    freqs = np.arange(shift,num_cells+shift) ** -alpha
    return freqs/sum(freqs)


