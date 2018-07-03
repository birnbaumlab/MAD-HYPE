
"""
Supplemnetary Figure 1
"""

# nonstandard libraries
import matplotlib.pyplot as plt
from scipy.optimize import minimize,root
import numpy as np

# homegrown libraries
import madhype
from madhype.simulation import *  # we should explicitly import here
from madhype import simulate_run

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 16 

def main(*args,**kwargs):
    """ Loads data, does some pretty basic simulations """

    fnames = {
            'Adjacent Tissue':'PTC.txt',
            'Tumor Tissue':'TTC.txt'
            }

    specific_settings = {
            'Adjacent Tissue':{
                'num_cells':1000,
                'num_wells':(48,48),
                'cpw':(50,1000)
                },
            'Tumor Tissue':{
                'num_cells':1000,
                'num_wells':(48,48),
                'cpw':(50,1000)
                }
            }
    
    # figure specific properties
    fig,axes = plt.subplots(nrows=2, ncols=len(fnames), figsize=(12,12), sharey=False)
    plt.tight_layout()
    #plt.subplots_adjust(left=0.15,right=0.9,top=0.85,bottom=0.3,hspace=0.5,wspace=0.5)
    fs = 18

    # iterate across sample/data sets
    for ind,(sample,fname) in enumerate(fnames.items()):
        
        ax = axes[ind][0]
        plt.sca(ax)

        with open('./figures/'+fname) as f:
            data = f.readlines()

        data = [float(x.strip()) for x in data] 
        data = [x/sum(data) for x in data]

        options = {'maxiters':10000}
        alpha = minimize(_error,1.0,data,method = 'Nelder-Mead',tol=1e-12,options=options)

        alpha = 1. + 1./alpha.x

        print 'Sample: {} / Alpha: {}'.format(sample,alpha)

        [i.set_linewidth(3) for i in axes[ind][0].spines.itervalues()]

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Clonal Index',fontsize=fs,fontweight='bold')
        ax.set_ylabel('Clonal Frequency',fontsize=fs,fontweight='bold')
        ax.tick_params(width = 3,length=8,labelsize=18)
        #ax.tick_params(axis='y',pad=22)

        ax.scatter(xrange(1,len(data)+1),
                data,c='blue',s=55,marker='s')
        ax.scatter(xrange(1,len(data)+1),
                _power_law_distribution(len(data),max(data),alpha),c='red',s=30,marker='D')

        plt.yticks(rotation='vertical')

        # default settings
        settings = {
                'cell_freq_constant': alpha,
                'cell_freq_max':    max(data),
                'plot_frequency_estimation': {
                    'ax':axes[ind][1],
                    'fig':fig,
                    },
                'visual':                True,
                'silent':               False,
                }

        # update settings with sample specific hyperparameters
        settings.update(specific_settings[sample])

        solvers = ['madhype']
        solver_options = [{}]

        data,results = madhype.simulate_run(solvers, solver_options, **settings)
        
    plt.subplots_adjust(left=0.125,right=0.9,top=0.9,bottom=0.1,hspace=0.4,wspace=0.4)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()

def _error(alpha,data): 
    simulation = _power_law_distribution(len(data),max(data),alpha)
    return sum([abs(a-b) for a,b in zip(data,simulation)])

        
        
#------------------------------------------------------------------------------# 

def _power_law_distribution(num_cells,max_freq,alpha):
    """ Returns power law distribution using given parameters """ 
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

#------------------------------------------------------------------------------# 

def _get_max_freq_diff(shift,num_cells,max_freq,alpha):
    """ Function for finding diff. b/w max_freq and current distribution """
    freqs = _get_freq_distribution(shift,num_cells,max_freq,alpha)
    return max_freq - freqs[0]

#------------------------------------------------------------------------------# 

def _get_freq_distribution(shift,num_cells,max_freq,alpha):
    """ Generate a normalized power-law distribution """
    freqs = np.arange(shift,num_cells+shift) ** -alpha
    return freqs/sum(freqs)

#------------------------------------------------------------------------------# 

if __name__ == '__main__':
    main()

