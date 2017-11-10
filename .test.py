from scipy.stats import binom
import numpy as np
import time


def generate_power_law_clones(a,n,shift):
    """ Generates n clones under a frequency distribution w/ alpha = a """ 
    freqs = 10.**(-a*np.log10(np.arange(1,n+1)))
    freqs = freqs/sum(freqs)
    return freqs


print generate_power_law_clones(1,10)
print generate_power_law_clones2(1,10)


