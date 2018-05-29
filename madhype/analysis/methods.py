
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
from operator import mul
from math import log10

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.optimize import brentq
from scipy.misc import comb,factorial

# homegrown libraries

#------------------------------------------------------------------------------# 
""" Main Callable Methods """
#------------------------------------------------------------------------------# 

def match_probability(well_data,prior = 1.0,memory={}):

    """ Calculates match probability given well_data and a prior ratio """

    key = (well_data['w_i']) + \
           (well_data['w_j']) + \
           (well_data['w_ij']) + \
           (well_data['w_o']) + \
           (well_data['w_tot']) + \
           (well_data['cpw']) + \
           (well_data['alpha'],)
    try:
        # backdoor check memory for probability
        return memory[key]

    except KeyError:
        """ MLE ESTIMATES """
        # generate MLE for clonal match
        freqs_match = estimate_match_frequencies(well_data)

        # generate MLE for clonal nonmatch
        freqs_nonmatch = estimate_nonmatch_frequencies(well_data)

        #
        if freqs_match['ij'] == 0.:
            return None,(freqs_match,freqs_nonmatch)

        """ PROBABILITY CALCULATIONS """
        # calculate probability for clonal match
        p_match = estimate_probability(well_data,freqs_match)


        # calculate probability for clonal nonmatch
        p_nonmatch = estimate_probability(well_data,freqs_nonmatch)


        #"""#
        #TESTING
        if False:#well_data['w_ij'] == (19,) and well_data['w_o'] == (36,):
            print 'W_i:',well_data['w_i']
            print 'W_j:',well_data['w_j']
            print 'W_ij:',well_data['w_ij']
            print 'W_o:',well_data['w_o']
            print 'Match freqs:'
            for k,v in freqs_match.items(): print '{}:{}'.format(k,v)
            print 'Nonmatch freqs:'
            for k,v in freqs_nonmatch.items(): print '{}:{}'.format(k,v)
            print 'Nonmatch:',p_nonmatch
            print 'Match:',p_match
            raw_input()
        #"""#

        memory[key] = prior*p_match/p_nonmatch,(freqs_match,freqs_nonmatch)
        return memory[key]


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
    freqs['ij'] = _find_freq({'w':w_ij_clonal,
                              'w_tot':_add(_subtract(data['w_tot'],data['w_ij']),w_ij_clonal),
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
                             'alpha':data['alpha']}),
             'ij':0.0}

    # return result
    return freqs

#------------------------------------------------------------------------------# 

def estimate_probability(data,freqs):

    """ Estimate match probability given data,freqs """

    p_total = 1.
    keys = ['w_i','w_j','w_ij','w_o','w_tot','cpw']
    for w_i,w_j,w_ij,w_o,w_tot,cpw in zip(*[data[k] for k in keys]):
        f_i,f_j,f_ij = _well_freq(freqs['i'],cpw),_well_freq(freqs['j'],cpw),_well_freq(freqs['ij'],cpw)
        p_total *= _multinomial(w_i,w_j,w_ij,w_o)* \
                    (f_i**w_i)*(1.-f_i)**(w_o + w_j)* \
                    (f_j**w_j)*(1.-f_j)**(w_o + w_i)* \
                    ((1 - (1-f_i*f_j)*(1-f_ij))**w_ij)*(1.-f_ij)**(w_i + w_j + w_o)
    return p_total


#------------------------------------------------------------------------------# 

def estimate_nonmatch_probability(data,freqs):

    """ Estimate nonmatch probability given data,freqs """

    p_total = 1.
    keys = ['w_i','w_j','w_ij','w_o','cpw']
    for w_i,w_j,w_ij,w_o,cpw in zip(*[data[k] for k in keys]):
        # calculate nonmatch probability
        p = _prob_distribution(w_i,w_j,w_ij,w_o,freqs,cpw)
        p_total *= p

    return p_total

#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _adjust_freq(f,cpw):
    """ Adjust repertoire frequency by cpw """
    return (1. - (1. - f)**cpw)


#------------------------------------------------------------------------------# 

def _w_ij_adjustment(w_ij,w_tot,cpw,freqs): 
    """ Remove the contribution of noise to w_ij counts """
    w_ij_est = [w_t*(1 - (1 - freqs['i'])**c)*(1 - (1 - freqs['j'])**c)
                for w_t,c in zip(w_tot,cpw)]
    w_ij_est = tuple(float(w) for w in _subtract(w_ij,w_ij_est)) # remove explained wells

    return tuple(max(0,w) for w in w_ij_est)



#------------------------------------------------------------------------------# 

def _prob_distribution(w_i,w_j,w_ij,w_o,freqs,cpw):
    """ Estimates probability based on given data/freqs """

    """
    #TESTING
    print 'bin1:',_binomial_pdf(w_i+w_j+w_ij+w_o,w_i+w_ij,freqs['i']) # place clone i on i,ij wells
    print 'v1:',w_i+w_j+w_ij+w_o,w_i+w_ij,freqs['i']
    print 'bin2:',_binomial_pdf(w_i+w_ij,w_ij,freqs['j']) # place clone j on ij wells
    print 'v2:',w_i+w_ij,w_ij,freqs['j']
    print 'bin3:',_binomial_pdf(w_j+w_o,w_j,freqs['j'])   # place clone j on j,o wells
    print 'v3:',w_j+w_o,w_j,freqs['j']
    """

    p = _binomial_pdf(w_i+w_j+w_ij+w_o,w_i+w_ij,_adjust_freq(freqs['i'],cpw)) # place clone i on wells 
    p *= _binomial_pdf(w_i+w_ij,w_ij,_adjust_freq(freqs['j'],cpw)) # place clone j on ij wells
    p *= _binomial_pdf(w_j+w_o,w_j,_adjust_freq(freqs['j'],cpw))   # place clone j on j,o wells
    return p

#------------------------------------------------------------------------------# 

def _find_freq(well_data,memory={}):
    """ Estimates single frequency given well data """
    key = (well_data['w']) + \
           (well_data['w_tot']) + \
           (well_data['cpw']) + \
           (well_data['alpha'],)
    
    try:
        return memory[key]
    except KeyError:
        memory[key] = \
                brentq(_derivative_prob_func,0,1,(well_data,),disp=False)
        return memory[key]

"""
def _find_freq_ez(well_data):
    x,val = 0.,1
    cw = _product(well_data['cpw'],well_data['w'])
    sum_cW = sum(_product(well_data['cpw'],well_data['w_tot']))
    while abs(val) > 1e-6:
        val = well_data['alpha']*(x/(1.-x)) + sum_cW - \
            sum((c*(1./(1.-x**ci)) for c,ci in cw,well_data['cpw']))

    return 1-x
"""

#------------------------------------------------------------------------------# 

def _derivative_prob_func(f,d):
    if f == 0.:
        return -(d['alpha']*(1.-f)) + \
            sum([f*c*(W-w) - w
                for c,W,w in zip(d['cpw'],d['w_tot'],d['w'])])
    elif f == 1.:
        return sum([f*c*(W-w)
                for c,W,w in zip(d['cpw'],d['w_tot'],d['w'])])
    else:
        return -(d['alpha']*(1.-f)) + \
                f*sum([c*((W-w) + w*((1.-f)**c)/(((1.-f)**max(c,0.001)) - 1)) # FIXME: need rigorous solution to max function
                for c,W,w in zip(d['cpw'],d['w_tot'],d['w'])])


def _prob_func(f,well_dict):
    ''' Outputs a probability for a given frequency estimate '''
    val = f**-well_dict['alpha']
    for cpw,w_tot,w in zip(well_dict['cpw'],well_dict['w_tot'],well_dict['w']):
        val *= ((1.-(1.-f)**cpw)**w)*(((1.-f)**cpw)**(w_tot-w))
    return -val

#------------------------------------------------------------------------------# 
""" Factory Methods """
#------------------------------------------------------------------------------# 

def _well_freq(f,c):
    return (1.-(1.-f)**c)

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

def _multinomial(*args):
    """ binomial pdf calculator, n choose k at f """
    return factorial(sum(args),exact=False)/reduce(mul,[factorial(arg,exact=False) for arg in args],1)

#------------------------------------------------------------------------------# 
""" Testing Methods """
#------------------------------------------------------------------------------# 

# namespace catch
if __name__ == '__main__':
    # unit tests

    mode = 'normal'

    if mode == 'find_freq':

        y = []
        x_span = np.logspace(-20,0,101)

        inputs = {'w':(1,),
                  'w_tot':(99,),
                  'cpw':(1,),
                  'alpha':1}

        for i in x_span:
            y.append(_derivative_prob_func(i,inputs))

        plt.plot(x_span,y)
        plt.xscale('log')
        plt.show(block=False)
        raw_input('Press enter to close...')
        plt.close()

    if mode == 'normal':
        # TESTING: estimate_match_frequencies
        well_data = {'w_i':(5,),
                     'w_j':(7,),
                     'w_ij':(23,),
                     'w_o':(61,),
                     'w_tot':(96,),
                     'cpw':(10,),
                     'alpha':1}
        #"""

        freqs_nonmatch = estimate_nonmatch_frequencies(well_data)
        freqs_match = estimate_match_frequencies(well_data)
        
        print 'Nonmatch frequency:',freqs_nonmatch
        print 'Match frequency:',freqs_match
        print 'Nonmatch probability:',estimate_probability(well_data,freqs_nonmatch)
        print 'Nonmatch probability:',estimate_nonmatch_probability(well_data,freqs_nonmatch)
        print 'Match probability:',estimate_probability(well_data,freqs_match)


