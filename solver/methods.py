
"""

Author: Patrick V. Holec
Created: 4 December 2017

"""

"""

Common solver methods for MAD-HYPE

Included:
    > Callable <
    estimate_match_frequencies(well_data)
    estimate_match_probability(well_data,freqs)
    estimate_nonmatch_frequencies(well_data)
    estimate_nonmatch_probability(well_data,freqs)
    
    > Internal <
    _find_freq(well_data)
    _prob_func(f,well_dict)
    _add(a,b)

"""

#------------------------------------------------------------------------------# 

""" Library Importation """

# standard libraries
import operator
from math import log10

# nonstandard libraries
import numpy as np
from scipy.optimize import fmin
from scipy.misc import comb

# homegrown libraries

#------------------------------------------------------------------------------# 
""" Main Callable Methods """
#------------------------------------------------------------------------------# 

@profile
def match_probability(well_data,prior = 1.0):

    """ Calculates match probability given well_data and a prior ratio """

    freqs_match = estimate_match_frequencies(well_data)
    freqs_nonmatch = estimate_nonmatch_frequencies(well_data)
   
    p_match = estimate_match_probability(well_data,freqs_match)
    p_nonmatch = estimate_nonmatch_probability(well_data,freqs_nonmatch)

    if p_match == 0.:
        print well_data
        print p_nonmatch
        print p_match
        raw_input()
        return None
        
    return prior*p_nonmatch/p_match


def estimate_match_frequencies(data):

    """ Estimates fa,fb,fab given well data """

    # assertion check
    keys = ('w_i','w_j','w_ij','w_o','w_tot','cpw','alpha')
    assert isinstance(data,dict),'input incorrect type <{}>'.format(type(data))
    assert all(k in data for k in keys),'keys missing in dict'
    assert all(isinstance(data[k],(list,tuple)) for k in keys[:-1]),'entries in dict not all <tuple>'

    # identify noise rates
    freqs = {'i':_find_freq({'w':data['w_i'],
                             'w_tot':_add(data['w_i'],data['w_o']),
                             'cpw':data['cpw'],
                             'alpha':data['alpha']}),
             'j':_find_freq({'w':data['w_j'],
                             'w_tot':_add(data['w_j'],data['w_o']),
                             'cpw':data['cpw'],
                             'alpha':data['alpha']})}

    # find w_ij resulting from noise, estimate clonal f_ij after
    w_ij_clonal = _w_ij_adjustment(data['w_ij'],data['w_tot'],data['cpw'],freqs)
    freqs['ij'] = _find_freq({'w':_subtract(data['w_ij'],w_ij_clonal),
                              'w_tot':_subtract(data['w_tot'],w_ij_clonal),
                              'cpw':data['cpw'],
                              'alpha':data['alpha']})


    # return result
    return freqs

#------------------------------------------------------------------------------# 

def estimate_nonmatch_frequencies(data):

    """ Estimates fa,fb given well data """

    # assertion check
    keys = ('w_i','w_j','w_ij','w_o','w_tot','cpw','alpha')
    assert isinstance(data,dict),'input incorrect type <{}>'.format(type(data))
    assert all(k in data for k in keys),'keys missing in dict'
    assert all(isinstance(data[k],(list,tuple)) for k in keys[:-1]),'entries in dict not all <tuple>'

    # identify noise rates
    freqs = {'i':_find_freq({'w':_add(data['w_i'],data['w_ij']),
                             'w_tot':data['w_tot'],
                             'cpw':data['cpw'],
                             'alpha':data['alpha']}),
             'j':_find_freq({'w':_add(data['w_j'],data['w_ij']),
                             'w_tot':data['w_tot'],
                             'cpw':data['cpw'],
                             'alpha':data['alpha']})}

    # return result
    return freqs

#------------------------------------------------------------------------------# 

def estimate_match_probability(data,freqs):

    """ Estimate match probability given data,freqs """

    p_total = 0. 
    keys = ['w_i','w_j','w_ij','w_o','w_tot']
    for w_i,w_j,w_ij,w_o,w_tot in zip(*[data[k] for k in keys]):
        for w_ij_clonal in xrange(0,w_ij):
            # calculate match probability
            p_total += _binomial_pdf(w_tot,w_ij_clonal,freqs['ij'])* \
                       _prob_distribution(w_i,w_j,w_ij-w_ij_clonal,w_o,freqs)

    return p_total


#------------------------------------------------------------------------------# 

def estimate_nonmatch_probability(data,freqs):

    """ Estimate nonmatch probability given data,freqs """
    
    p_total = 0. 
    keys = ['w_i','w_j','w_ij','w_o']
    for w_i,w_j,w_ij,w_o in zip(*[data[k] for k in keys]):
        # calculate nonmatch probability
        p_total += _prob_distribution(w_i,w_j,w_ij,w_o,freqs)

    return p_total

#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _w_ij_adjustment(w_ij,w_tot,cpw,freqs):
    """ Remove the contribution of noise to w_ij counts """
    w_ij_est = [w_t*(1 - (1 - freqs['i'])**c)*(1 - (1 - freqs['j'])**c) 
                for w_t,c in zip(w_tot,cpw)]
    w_ij_est = _subtract(w_ij,w_ij_est) # remove explained wells
    return tuple(max(0,w) for w in w_ij_est) 



#------------------------------------------------------------------------------# 

def _prob_distribution(w_i,w_j,w_ij,w_o,freqs):
    """ Estimates probability based on given data/freqs """
    p = _binomial_pdf(w_i+w_j+w_ij+w_o,w_i+w_ij,freqs['i']) # place clone i on i,ij wells
    p *= _binomial_pdf(w_i+w_ij,w_ij,freqs['j']) # place clone j on ij wells
    p *= _binomial_pdf(w_j+w_o,w_j,freqs['j'])   # place clone j on j,o wells
    return p

#------------------------------------------------------------------------------# 

def _find_freq(well_data):
    """ Estimates single frequency given well data """
    return (fmin(_prob_func,0.01,(well_data,),disp=False)) # maximize PDF

#------------------------------------------------------------------------------# 

def _prob_func(f,well_dict):
    ''' Outputs a probability for a given frequency estimate '''
    val = f**-well_dict['alpha']
    for cpw,w_tot,w in zip(well_dict['cpw'],well_dict['w_tot'],well_dict['w']):
        val *= ((1.-(1.-f)**cpw)**w)*(((1.-f)**cpw)**(w_tot-w))
    return -val

#------------------------------------------------------------------------------# 
""" Factory Methods """
#------------------------------------------------------------------------------# 

def _add(a,b):
    """ add iterables a,b """
    return tuple(map(operator.add, a, b))

def _subtract(a,b):
    """ add iterables a,b """
    return tuple(map(operator.sub, a, b))

def _product(a,scalar):
    """ multiply iterable a by scalar """
    return tuple([scalar*x for x in a])

def _binomial_pdf(n,k,f):
    """ binomial pdf calculator, n choose k at f """
    return comb(n,k,exact=False)*(f**k)*((1.-f)**(n-k))

#------------------------------------------------------------------------------# 
""" Testing Methods """
#------------------------------------------------------------------------------# 

# namespace catch
if __name__ == '__main__':
    # unit tests

    # TESTING: estimate_match_frequencies
    well_data = {'w_i':(20,20),
                 'w_j':(30,50),
                 'w_ij':(30,25),
                 'w_o':(20,5),
                 'w_tot':(100,100),
                 'cpw':(10,30),
                 'alpha':1}

    well_data = {'w_i':(16,),
                 'w_j':(16,),
                 'w_ij':(4,),
                 'w_o':(64,),
                 'w_tot':(100,),
                 'cpw':(1,),
                 'alpha':1}

    freqs_match = estimate_match_frequencies(well_data)
    freqs_nonmatch = estimate_nonmatch_frequencies(well_data)
    
    print freqs_match
    print freqs_nonmatch
    print estimate_match_probability(well_data,freqs_match)
    print estimate_nonmatch_probability(well_data,freqs_nonmatch)

    print 'Match probability:',match_probability(well_data)
