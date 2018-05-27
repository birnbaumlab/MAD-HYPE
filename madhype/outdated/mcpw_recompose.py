
"""
This program is designed to try and figure out if there's
ways to predict frequency well in multi-cpw scenarios
"""

# standard libraries

# nonstandard libraries
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fmin 

# homegrown libraries
from variable.methods import *

# library modifications
plt.rcParams["font.family"] = "serif"

#----------------------------------------------------------------#

def test(*args,**kwargs):

    settings = {
                'distribution':'power-law',
                'alpha':1.0,
                'clone_count':100000,
                'cpw_distro':[(1000,96)],
                #'cpw_distro':[(10,48),(100,48)],
                'pattern_repeats':1,
                'well_total':96
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    # change cpw_distro into dictionary if not there already
    if all([isinstance(s,tuple) for s in settings['cpw_distro']]):
        settings['cpw_distro'] = [{'cpw':a,'well_count':b} for a,b in settings['cpw_distro']]

    # generate power law frequencies
    freqs = generate_power_law_clones(1,100,0)

    # generate patterns from clones
    patterns = generate_simulated_patterns(freqs,settings['cpw_distro'],settings['pattern_repeats'])
    
    est_freqs = []

    for i,f,p in zip(xrange(len(freqs)),freqs,patterns[0]):
        est = estimate(p,settings)
        print 'Clone {}:'.format(i)
        print '> Clonal frequency - {}'.format(f)
        print '> Pattern:'.format(p)
        print '> Estimated frequency:',est
        est_freqs.append(est)

    plt.plot((min(freqs),1),(min(freqs),1),'-b')
    plt.xscale('log')
    plt.yscale('log')
    plt.scatter(freqs,est_freqs)
    plt.xlabel('Actual clonal frequencies')
    plt.ylabel('Predicted clonal frequencies')
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()

def estimate(p,settings):
    ''' Gets a frequency estimate for a given pattern '''
    ff = FindFreq(p,settings)
    return fmin(ff.func,0.0)
    
class FindFreq:
    def __init__(self,pattern,settings):
        ''' Initialization with settings,pattern '''
        self.pattern = pattern
        self.alpha = settings['alpha']
        self.cpws = [i['cpw'] for i in settings['cpw_distro']]
        self.w_tots = [i['well_count'] for i in settings['cpw_distro']]
    
    def func(self,f):
        ''' Outputs a probability for a given frequency estimate '''
        val = f**-self.alpha 
        for cpw,w_tot,w_i in zip(self.cpws,self.w_tots,self.pattern):
            val *= ((1.-(1.-f)**cpw)**w_i)*(((1.-f)**cpw)**(w_tot-w_i))
        return -val


if __name__ == "__main__":
    test()
