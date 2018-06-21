
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
plt.rcParams['xtick.labelsize'] = 16

def main(*args,**kwargs):
    """ Loads data, does some pretty basic simulations """

    fnames = {
            'Adjacent Tissue':'PTC.txt',
            'Tumor Tissue':'TTC.txt'
            }
    

    # figure specific properties
    fig,axes = plt.subplots(nrows=1, ncols=len(fnames), figsize=(14,6), sharey=False)
    plt.subplots_adjust(left=0.15,right=0.9,top=0.85,bottom=0.3,hspace=0.5,wspace=0.5)
    fs = 18

    specific_settings = {
            # ranges 0.1% to 3%
            'Peripheral Blood (9-25 y)':{
                'cell_freq_max':0.0025,
                'alpha':1.2,
                'num_cells':10000,
                'threshold':2.0,
                'block':False,
                'num_wells':(48,48),
                'cpw':(500,10000)
                },
            'Peripheral Blood (61-66 y)':{
                'cell_freq_max':0.086,
                'alpha':1.2,
                'num_cells':10000,
                'threshold':2.0,
                'block':True,
                'num_wells':(48,48),
                'cpw':(500,10000)
                },
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

    # iterate across sample/data sets
    for ind,(sample,fname) in enumerate(fnames.items()):
        
        with open('./figures/'+fname) as f:
            data = f.readlines()

        data = [float(x.strip()) for x in data] 
        data = [x/sum(data) for x in data]

        options = {'maxiters':10000}
        alpha = minimize(_error,1.0,data,method = 'Nelder-Mead',tol=1e-12,options=options)
        print 'Sample: {} / Alpha: {}'.format(sample,alpha.x)

        axes[ind].scatter(xrange(1,len(data)+1),
                data,c='blue',s=55,marker='s')
        axes[ind].scatter(xrange(1,len(data)+1),
                _power_law_distribution(len(data),max(data),alpha.x),c='red',s=30,marker='D')

        [i.set_linewidth(3) for i in axes[ind].spines.itervalues()]

        axes[ind].set_xscale('log')
        axes[ind].set_yscale('log')
        axes[ind].set_xlabel('Clonal Index',fontsize=fs)
        axes[ind].set_ylabel('Clonal Frequency',fontsize=fs)
        axes[ind].set_title(sample,fontsize=fs+4)

        settings = {
                'cell_freq_constant':alpha.x,
                'cell_freq_max':max(data),
                'visual':True,
                'display':True,
                'silent':False
                }

        solvers = ['madhype']
        solver_options = [{}]

        data,results = madhype.simulate_run(solvers, solver_options, **settings)
        

    settings = {
        'analysis':('madhype',),
        'visual':True,
        'silent':False
        }

    print 'Starting finals...'
    #simulate_system(settings,**specific_settings['Peripheral Blood (9-25 y)'])
    #simulate_system(settings,specific_settings['Peripheral Blood (61-66 y)'])


    '''
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig1S.png', format='png', dpi=300)
    plt.close()
    #'''

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

