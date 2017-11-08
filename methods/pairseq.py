
'''

Project: 
pairSEQ (something-something Howie)

Class(s): 
(none)

Function: 
Alternative Bayesian system for evaluating paired chain information

'''

'''
Library Importation
'''


# standard libraries
import math
import time
from math import log

# nonstandard library import
import numpy as np
from scipy.misc import comb
from scipy.special import factorial
import matplotlib.pyplot as plt



'''
Factory Methods
'''

def log_factorial_approx(n):
    return log((2*3.1415*n)**0.5) + n*log(n/2.718)

# Mean of PDF
def mean_pdf(x,y):
    return sum([i*j for i,j in zip(x,y)])

# Variance of PDF
def var_pdf(x,y):
    avg = mean_pdf(x,y)
    return sum([((i-avg)**2)*j for i,j in zip(x,y)])

# N choose K
def nCk(n, r):
    return comb(n,r)

# Beta Function
def beta(x,y):
    #return math.gamma(x)*math.gamma(y)/math.gamma(x + y)
    return math.exp(math.lgamma(x) + math.lgamma(y) - math.lgamma(x+y))

# Beta function parameters from mean and variance
def beta_parameters(error,samples):
    return (1.-error)*samples,error*samples

# Stirling number generator (2nd kind)
def stirling(n,k,memory={}):
    if not (n,k) in memory.keys(): 
        n1=n
        k1=k
        if n<=0:
            return 1
        elif k<=0:
            memory[(n,k)] = 0
        elif (n==0 and k==0):
            memory[(n,k)] = -1
        elif n!=0 and n==k:
            memory[(n,k)] = 1
        elif n<k:
            memory[(n,k)] = 0
        else:
            temp1=stirling(n1-1,k1)
            temp1=k1*temp1
            memory[(n,k)] = (k1*(stirling(n1-1,k1)))+stirling(n1-1,k1-1)
    return memory[(n,k)]




'''
Howie Equations
'''

def eq2(c,f,C):
    return nCk(C,c)*(f**c)*((1-f)**(C-c))

def eq3(cx,c,s1,s2):
    return (nCk(c,cx)*beta(cx+s1,c-cx+s2)/beta(s1,s2))

def eq4(w,c,W):
    a = log_factorial_approx(W) + log(stirling(c,w))
    b1 = log(factorial(W-w))
    b2 = c*log(W)
    return np.exp(a - b1 - b2)

def eq6(w_ab,w_a,w_b,w):
    return (nCk(w,w_ab)*nCk(w-w_ab,w_a-w_ab)*nCk(w-w_a,w_b-w_ab))/(nCk(w,w_a)*nCk(w,w_b))


def pairseq(w_ab,w_a,w_b,f,C,W,s1,s2,threshold = 1e-3,silent=True):
    p_total = 0. # prep a summing variable
    
    for c_i in xrange(0,C + 1):
        if not silent: print 'Identifying match probability... {}%\r'.format(100*(c_i+1)/C),
            
        if c_i < max(w_ab,w_a,w_b): # Check on the status of w_i
            continue
        
        p1 = eq2(c_i,f,C)
        
        if p1 < threshold: # Check on the status of p1
            continue # condition for early exit (apply to p)
            
        for w_i in xrange(0,1 + min([c_i,W])):
            if w_i < max(w_a,w_b): # Check on the status of w_i      
                continue
            p2 = eq4(w_i,c_i,W)
            
            if p1*p2 < threshold: # Check on the status of p2
                continue # condition for early exit (apply to p_term)
            
            for c_ia in xrange(0,c_i + 1):
                
                if c_ia < w_a: # Check on the status of c_ia
                    continue # you need atleast c cells to find a in w_a wells
                    
                p3 = eq3(c_ia,c_i,s1,s2)
                    
                if p1*p2*p3 < threshold: # Check on the value of p3
                    continue # condition for early exit (apply to p_term)
                
                for c_ib in xrange(0,c_i + 1):
                    
                    if c_ib < w_b: # Check on the status of c_ib
                        continue # you need atleast c cells to find b in w_b wells
                    
                    p4 = eq3(c_ib,c_i,s1,s2)

                    if p1*p2*p3*p4 < threshold: # Check on the status of p4
                        continue # condition for early exit (apply to p_term)
                    
                    # Normally we would integrate across all values w_a,w_b,w_ab
                    p5 = eq4(w_a,c_ia,W)
                    p6 = eq4(w_b,c_ib,W)

                    if p1*p2*p3*p4*p5*p6 < threshold:
                        continue
                    
                    p7 = eq6(w_ab,w_a,w_b,w_i)
                    p8 = 1
                    
                    if all([p5,p6,p7,p8]):
                        p_total += p1*p2*p3*p4*p5*p6*p7*p8
    
    if not silent: print ''
    return p_total,timer


'''
The Meat-and-Potatoes Function
'''

# Verbose: on range 0 to 9
# TODO: verbose levels

def solve(data,pair_threshold = 0.99,verbose=0):
    
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))
    
    freqs = data.metadata['generated_data']['cell_frequencies']
    cells = data.metadata['cells']
    
    b_freqs = dict([(b,sum([freqs[i] for i,c in enumerate(cells) if b in c[1]])) for b in b_uniques])
    
    # sort in ascending order
    sorted(a_uniques,key=int)
    sorted(b_uniques,key=int)
    
    a_dict = dict([(a,[]) for a in a_uniques]) # create dictionary for each clonotype
    b_dict = dict([(b,[]) for b in b_uniques]) # create dictionary for each clonotype
    
    for w,well in enumerate(data.well_data):
        [a_dict[i].append(ind) for ind,i in enumerate(well[0])]
        [b_dict[i].append(ind) for ind,i in enumerate(well[1])]
    
    # Find generalized parameters for pairSEQ
    W = data.metadata['num_wells']
    C = W*data.metadata['cells_per_well_distribution_params']['cells_per_well']
    error = data.metadata['chain_deletion_prob']
    samples = len(a_uniques) + len(b_uniques)
    s1,s2 = beta_parameters(error,samples)
    
    # Iterate across all possible calculations of pairseq
    for a in a_uniques:
        a_wells = a_dict[a]
        for b in b_uniques:
            b_wells = b_dict[b]
            numel = len(set(a_wells).intersection(b_wells))
            w_ab,w_a,w_b = numel,len(a_wells)-numel,len(b_wells)-numel
            f = b_freqs[b]
            print 'Combination ({},{}): {}'.format(a,b,pairseq(w_ab,w_a,w_b,f,C,W,s1,s2,threshold = 1e-3,silent=False))

    print 'Finished!'

#error,samples = 0.1,1000
#s1,s2 = beta_parameters(error,samples)
#p,timer = pairseq(12,13,13,0.05,500,50,s1,s2,threshold=1e-6)
#def pairseq(w_ab,w_a,w_b,f,C,W,s1,s2,threshold = 1e-3,silent=True):


        
    
    
    
    

