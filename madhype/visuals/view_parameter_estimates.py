
### This program attempts to estimate the parameter regimes for full-repertoire and TIL sequencing ###

# standard libraries


# nonstandard libraries
import numpy as np
from matplotlib import pyplot as plt
from datasim import *
from scipy.misc import comb
from scipy.stats import binom

# homegrown libraries


# library modifications
plt.rcParams["font.family"] = "serif"


def view_parameter_estimates(*args,**kwargs):

    distributions = {'full-repertoire':1.16}
    settings = {
                'clone_count':10000,
                'cell_per_well_distribution':'constant',
                'cells_per_well':1000,
                'well_total':96,
                'error_repeats':20,
                'success_repeats':1000
               }
    
    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)
    
    for sample,alpha in distributions.items():
        
        print 'Sample:',sample
        cells_per_well = [1,10,100,1000]
        labels = ['{} cells/well'.format(c) for c in cells_per_well]
        styles = ['-','--','-.',':']

        freqs = generate_power_law_clones(alpha,settings['clone_count'])

        for l,s,cpw in zip(labels,styles,cells_per_well):
            settings['cells_per_well'] = cpw
            error_dict = generate_error_rate_dict(freqs,settings)
            successes = generate_success_rate_list(freqs,error_dict,settings)
            plt.plot(freqs,successes,s,label=l,linewidth = 2.0)
            print 'Finished {} cells/well!'.format(cpw)

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Clonal Frequency')
        plt.ylabel('Successful Identification Rate')
        plt.xlim((min(freqs),max(freqs)))
        plt.legend(loc='lower right')

        plt.show(block = False)
        raw_input('Press enter to close...')
        plt.close()


### Factory Methods

'''
def binomial_probability(w_t,f_well,partitions = [96]):
    """ Returns a <c/w subset> by <w_t + 1> binomial matrix """
    result = [1. for i in
    for i,partition in enumerate(partitions):
        distro = [binom.pmf(xrange(partition + 1),partition,f) for f in f_well]
        if i != 0:
'''
            

def generate_success_rate_list(freqs,error_dict,settings):
    
    success_rates = np.zeros((len(freqs),)) 
    w_t = settings['well_total']
    errors = [error_dict['avg'][i] for i in xrange(w_t + 1)]
    f_well = calculate_well_freqs(freqs,settings['cells_per_well'])

    for i,f in enumerate(freqs):
        f_well = 1. - (1.-f)**settings['cells_per_well']
        probs = binom.pmf(xrange(w_t + 1),w_t,f_well)
        success_rates[i] = sum([e*p for e,p in zip(errors,probs)])

    return success_rates 

def error_rate(count,wells,well_total):
    if count == 0: return 1.0
    return (1./(1. + (count-1)*(1./comb(well_total,wells))))

def generate_error_rate_dict(freqs,settings):
    count_dict_list = []
    w_t = settings['well_total']

    # create a bunch of count dicts
    for i in xrange(settings['error_repeats']):
        count_dict_list.append(generate_well_count_dict(freqs,settings))
    
    total_count = dict([(i,sum([c[i] for c in count_dict_list])) 
                        for i in  xrange(w_t + 1)]) 

    error_dict_avg = dict([(i,np.sum([c[i]*error_rate(c[i],i,w_t) 
                            for c in count_dict_list])/total_count[i]) 
                            for i in  xrange(w_t + 1)]) 
    error_dict_std = dict([(i,np.std([error_rate(c[i],i,w_t) 
                            for c in count_dict_list])) 
                            for i in  xrange(w_t + 1)]) 

    error_dict_avg = dict([(k,v) if not np.isnan(v) else (k,1.) for k,v in error_dict_avg.items()])
    error_dict_std = dict([(k,v) if not np.isnan(v) else (k,1.) for k,v in error_dict_std.items()])

    return {'avg':error_dict_avg,'std':error_dict_std}

def calculate_well_freqs(freqs,cpw):
    return [1. - ((1. - freqs)**s for s in cpw)]

def generate_well_count(freqs,settings):
    well_freqs = calculate_well_freqs(freqs,settings['cells_per_well'])
    return np.random.binomial(settings['well_total'],well_freqs,len(well_freqs))

def generate_well_count_dict(freqs,settings):
    x = generate_well_count(freqs,settings)
    y = np.bincount(x)
    result = dict([(i,0) for i in xrange(0,settings['well_total'] + 1)])
    for i,val in enumerate(y): result[i] = val 
    return result

def generate_power_law_clones(alpha,count):
    freqs = 10.**(-alpha*np.log10(np.arange(1,count+1)))
    freqs = freqs/sum(freqs)
    return freqs

def generate_bins(freqs):
    return np.logspace(np.log10(min(freqs)),np.log10(max(freqs)),6)


