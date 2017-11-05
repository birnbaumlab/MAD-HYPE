
from scipy.stats import binom
import numpy as np

def _calculate_well_freqs(freqs,cpw):
    return 1. - ((1. - freqs)**cpw)

f = np.array([0.2,0.4,0.6])



