
### This program attempts to estimate the parameter regimes for full-repertoire and TIL sequencing ###

# standard libraries


# nonstandard libraries
import numpy as np
from matplotlib import pyplot as plt
from datasim import *

# homegrown libraries



#distributions = {'full-repertoire':1.0,'TIL':2.0}
def view_parameter_estimates():


    #freqs = np.random.pareto(2.0, 1000) 
    #freqs = freqs/sum(freqs)
    freqs = 10.**(-3*np.log10(np.arange(1,100001)))
    freqs = freqs/sum(freqs)
    bins = np.logspace(np.log10(min(freqs)),np.log10(max(freqs)),6)

    sg= seq_generator.SequencingGenerator()
    
    print "bins: ", bins
    plt.xscale('log')
    plt.hist(freqs,bins=bins) 
    plt.yscale("log", nonposx='clip')
    plt.show()
    

#for sample,scale in distributions.items():
#    pass
    #generate power_law_clones(scale)


def generate_power_law_clones(slope):
    pass    

