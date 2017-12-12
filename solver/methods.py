
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
import matplotlib.pyplot as plt
from scipy.optimize import fmin
from scipy.optimize import brentq
from scipy.misc import comb

# homegrown libraries

#------------------------------------------------------------------------------# 
""" Main Callable Methods """
#------------------------------------------------------------------------------# 

def match_probability(well_data,prior = 1.0,memory={}):

    """ Calculates match probability given well_data and a prior ratio """

    #print well_data

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
    p_match = estimate_match_probability(well_data,freqs_match)

    # calculate probability for clonal nonmatch
    p_nonmatch = estimate_nonmatch_probability(well_data,freqs_nonmatch)

    #"""#
    #TESTING
    if prior*p_match/p_nonmatch > 1e35:
        #if p_match > 1.0 or p_nonmatch > 1.0:
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

    return prior*p_match/p_nonmatch,(freqs_match,freqs_nonmatch)


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
                             'alpha':data['alpha']})}

    # return result
    return freqs

#------------------------------------------------------------------------------# 

def estimate_match_probability(data,freqs):

    """ Estimate match probability given data,freqs """

    p_total = 1. 
    keys = ['w_i','w_j','w_ij','w_o','w_tot']
    for w_i,w_j,w_ij,w_o,w_tot in zip(*[data[k] for k in keys]):
        p = 0.
        #print 'Start new:'
        for w_ij_clonal in xrange(0,w_ij+1):
            # calculate match probability
            #print 'Val:',_binomial_pdf(w_tot,w_ij_clonal,freqs['ij'])* \
            #           _prob_distribution(w_i,w_j,w_ij-w_ij_clonal,w_o,freqs)
            p += _binomial_pdf(w_tot,w_ij_clonal,freqs['ij'])* \
                       _prob_distribution(w_i,w_j,w_ij-w_ij_clonal,w_o,freqs)
            #print 'Match p:',p
        p_total *= p

    return p_total


#------------------------------------------------------------------------------# 

def estimate_nonmatch_probability(data,freqs):

    """ Estimate nonmatch probability given data,freqs """
    
    p_total = 1. 
    keys = ['w_i','w_j','w_ij','w_o']
    for w_i,w_j,w_ij,w_o in zip(*[data[k] for k in keys]):
        # calculate nonmatch probability
        #p = _prob_distribution(w_i,w_j,w_ij,w_o,freqs)
        p = _prob_distribution(w_i,w_j,w_ij,w_o,freqs)
        p_total *= p
        #print 'NonMatch p:',p

    return p_total

#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _w_ij_adjustment(w_ij,w_tot,cpw,freqs):
    """ Remove the contribution of noise to w_ij counts """
    w_ij_est = [w_t*(1 - (1 - freqs['i'])**c)*(1 - (1 - freqs['j'])**c)
                for w_t,c in zip(w_tot,cpw)]
    w_ij_est = tuple(float(w) for w in _subtract(w_ij,w_ij_est)) # remove explained wells

    return tuple(max(0,w) for w in w_ij_est)



#------------------------------------------------------------------------------# 

def _prob_distribution(w_i,w_j,w_ij,w_o,freqs):
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

    p = _binomial_pdf(w_i+w_j+w_ij+w_o,w_i+w_ij,freqs['i']) # place clone i on i,ij wells
    p *= _binomial_pdf(w_i+w_ij,w_ij,freqs['j']) # place clone j on ij wells
    p *= _binomial_pdf(w_j+w_o,w_j,freqs['j'])   # place clone j on j,o wells
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
            f*sum([c*((W-w) + w*((1.-f)**c)/(((1.-f)**c) - 1))
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
        well_data = {'w_i':(3,0),
                     'w_j':(13,11),
                     'w_ij':(3,36),
                     'w_o':(29,1),
                     'w_tot':(48,48),
                     'cpw':(10,100),
                     'alpha':0}

        """
        well_data = {'w_i':(16,),
                     'w_j':(16,),
                     'w_ij':(4,),
                     'w_o':(64,),
                     'w_tot':(100,),
                     'cpw':(1,),
                     'alpha':1}
        """

        freqs_match = estimate_match_frequencies(well_data)
        freqs_nonmatch = estimate_nonmatch_frequencies(well_data)
        

        print freqs_match
        print freqs_nonmatch
        val = estimate_match_probability(well_data,freqs_match)
        print val
        print estimate_nonmatch_probability(well_data,freqs_nonmatch)

        print 'Match probability:',match_probability(well_data)

