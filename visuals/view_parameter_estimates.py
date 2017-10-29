
### This program attempts to estimate the parameter regimes for full-repertoire and TIL sequencing ###

# standard libraries


# nonstandard libraries
import numpy as np
from matplotlib import pyplot as plt
from datasim import *
from scipy.misc import comb

# homegrown libraries



def view_parameter_estimates(*args,**kwargs):

    distributions = {'full-repertoire':1.0,'TIL':2.0}
    settings = {
                'clone_count':100000,
                'cells_per_well':1000,
                'well_total':48
               }
    
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    for sample,alpha in distributions.items():

        freqs = generate_power_law_clones(alpha,settings['clone_count'])
        bins = generate_bins(freqs)
        count_dict = generate_well_count_dict(freqs,settings)

        print 'Population:',sample
        print 'Observed clones:',sum(count_dict.values())-count_dict[0]

        for k,v in count_dict.items():
            if v == 0: continue
            print '>',k,' - ',(1./(1. + (v-1)*(1./comb(settings['well_total'],k))))
        

        """
        plt.xscale('log')
        plt.hist(freqs,bins=bins) 
        plt.title(sample)

        plt.show(block = False)
        raw_input('Press enter to close...')
        plt.close()
        """

### Factory Methods

def generate_well_count_dict(freqs,settings):
        well_freqs = 1. - ((1. - freqs)**settings['cells_per_well'])
        x = np.random.binomial(settings['well_total'],well_freqs,len(well_freqs))
        y = np.bincount(x)
        ii = np.nonzero(y)[0]
        return dict(zip(ii,y[ii]))

def generate_power_law_clones(alpha,count):
    freqs = 10.**(-alpha*np.log10(np.arange(1,count+1)))
    freqs = freqs/sum(freqs)
    return freqs


def generate_bins(freqs):
    return np.logspace(np.log10(min(freqs)),np.log10(max(freqs)),6)

