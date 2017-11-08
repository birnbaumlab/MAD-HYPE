
'''

Project: 
MAD-HYPE (Multicell Analytical Deconvolution for High Yield Paired-chain Evaluation)

Update (8/29/2017):
    Modified script to create dictionaries represeting each well, allows circumventation of checking many sequences

Class(s): 
(none)

Function: 
Using analytics based in Bayesian systems, we attempt to identify T-cell clones from low-throughput well sequencing, with each well containing high-throughput numbers of cells

'''

'''
Library Importation
'''

# standard libraries
import math
import operator as op
import time
import os
import pickle
import itertools

# nonstandard libraries
import numpy as np
import scipy.misc
#import matplotlib.pyplot as plt
#from seq_generator import SequencingGenerator as SeqGen
#from seq_data import SequencingData


'''
Factory Methods
'''

# N choose K (OUTDATED, slower by ~25x)
def nCk_old(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

# N choose K
def nCk(n, r):
    return scipy.misc.comb(n,r)
    
'''
Frequency Estimators
'''

# non-match MLE estimator for f_a,f_b,f_ab
def nonmatch_frequency(w_ab,w_a,w_b,w_tot):
    return float(w_ab)/w_tot,float(w_a+w_ab)/w_tot,float(w_b+w_ab)/w_tot

# MLE estimator for f_a,f_b,f_ab
def match_frequency(w_ab,w_a,w_b,w_tot):
    if w_tot-w_ab == 0: f_a = 0
    else: f_a = float(w_a)/(w_tot-w_ab)
    if w_tot-w_ab == 0: f_b = 0
    else: f_b = float(w_b)/(w_tot-w_ab)
    f_ab = max((0,1. - (1.-(float(w_ab)/w_tot))/(1-f_a*f_b)))
    return f_ab,f_a,f_b

'''
Probability Calculators
'''

# prior probability
def nonmatch_probability(w_ab,w_a,w_b,w_tot):
    w_ab,w_a,w_b,w_tot = int(w_ab),int(w_a),int(w_b),int(w_tot)
    f_ab,f_a,f_b = nonmatch_frequency(w_ab,w_a,w_b,w_tot)
    prob = instance_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b)
    if prob == 0.: 
        return float('-inf')
    return math.log10(prob)

def match_probability(w_ab,w_a,w_b,w_tot):
    if w_ab == 0: return float('nan')
    w_ab,w_a,w_b,w_tot = int(w_ab),int(w_a),int(w_b),int(w_tot)
    f_ab,f_a,f_b = match_frequency(w_ab,w_a,w_b,w_tot)
    prob_total = 0
    for w in xrange(0,int(w_ab)+1):
        prob_total += (nCk(w_tot,w)*(f_ab**(w))*((1-f_ab)**(w_tot-(w))))*instance_probability(w_ab-w,w_a,w_b,w_tot-w,f_ab,f_a,f_b)        
    if prob_total == 0.: return float('-inf')
    return math.log10(prob_total)

def instance_probability(w_ab,w_a,w_b,w_tot,f_ab,f_a,f_b):
    a = nCk(w_tot,w_a+w_ab)*(f_a**(w_a+w_ab))*((1-f_a)**(w_tot-(w_a+w_ab)))
    b = nCk(w_a+w_ab,w_ab)*(f_b**w_ab)*((1-f_b)**(w_a+w_ab-w_ab))
    c =  nCk(w_tot-(w_a+w_ab),w_b)*(f_b**w_b)*((1-f_b)**(w_tot-(w_a+w_b+w_ab)))
    return a*b*c

'''
Specialty Calculators
'''
    
def self_occurence(v):
    temp = np.matmul(v,np.transpose(v))
    np.fill_diagonal(temp,0)
    return temp

def match_score(w_ab,w_a,w_b,w_tot, match_prior = 0.5):
    if w_ab == 0: 
        return 0.,0.
    else:
        mp,nmp = match_probability(w_ab,w_a,w_b,w_tot),nonmatch_probability(w_ab,w_a,w_b,w_tot)
        f_ab,_,_ = match_frequency(w_ab,w_a,w_b,w_tot)
        return (10**mp * match_prior)/(10**mp * match_prior + 10**nmp * (1-match_prior)),f_ab

# TODO: what do I even do?
def possible_matches(real_matches,img_ab,a_uniques,b_uniques):
    w_ab = np.sum(img_ab[:,:,:],axis=2)
    real = []
    for rm in real_matches:
        try: 
            i,j = a_uniques.index(rm[0]),b_uniques.index(rm[1])
            if w_ab[i,j] > 0: real.append(rm)
        except ValueError: # if a real match never even occurs
            pass
    return real

def precalculate_match_scores(w_tot,match_prior = 0.5):
    scores,freqs = np.zeros((w_tot+1,w_tot+1,w_tot+1)),np.zeros((w_tot+1,w_tot+1,w_tot+1))
    print ''
    for w_ab in xrange(w_tot+1):
        for w_a in xrange(w_tot+1):
            for w_b in xrange(w_tot+1):
                if w_ab + w_a + w_b > w_tot: 
                    continue
                else:
                    scores[w_ab][w_a][w_b],freqs[w_ab][w_a][w_b] = match_score(w_ab,w_a,w_b,w_tot,match_prior) # effectively take out prior
        #print 'Finished {}/{}...'.format(w_ab,w_tot)
    return scores,freqs

'''
Important Functions
'''

'''

Name: directional_matches
Function: Returns edges connecting matches that occur in group a -> b

Parameters:
    img_ab -
        description: co-occurence matrix for unique members of a and b, in each well
        size: (# unique a) x (# unique b) x (# wells)
    img_a -
        description: occurence matrix for unique members of a, in each well
        size: (# unique a) x (# wells)
    img_b -
        description: occurence matrix for unique members of b, in each well
        size: (# unique b) x (# wells)
    a_uniques -
        description: labels for unique members of a
        size: (# unique a) x (1)
    b_uniques -
        description: labels for unique members of a
        size: (# unique a) x (1)       
    threshold -
        description: crieria for Bayesian probability of match
        size: float in range (0 to 1)
    silent -
        description: whether function prints matches as it finds them
        size: boolean (True/False)   
Returns:
    predicted_ab - 
        description: list of tuples (a label,b label)
        
'''
# TODO: add stochasicity to well dismissal

def directional_matches(a_wells,b_wells,a_uniqs,b_uniqs,w_tot,threshold=0.99,silent=False,distinct=False):

    # important storage variables
    predicted_ab = []
    predicted_frequency = []
    predicted_score = []

    print 'Making C++ data files...'
    print a_wells
    
    with open('chain_data_a.txt','w') as f:
        for w in a_uniqs: f.write('{}'.format(str(a_wells[w]))[1:-1]+'\n')
    with open('chain_data_b.txt','w') as f:
        for w in b_uniqs: f.write('{}'.format(str(b_wells[w]))[1:-1]+'\n')
    with open('uniques_a.txt','w') as f:
        for w in a_uniqs: f.write('{}\n'.format(w))
    with open('uniques_b.txt','w') as f:
        for w in b_uniqs: f.write('{}\n'.format(w))

    raw_input('Hold...')
        
    # create dictionaries that map each well dictionary to the reverse (well # -> chains)
    wells_a = dict([(i,set()) for i in xrange(w_tot)])
    wells_b = dict([(i,set()) for i in xrange(w_tot)])

    # iterate across both dictionaries
    print 'Starting reverse well dictionary formation...'
    for wells_x,x_wells in [(wells_a,a_wells),(wells_b,b_wells)]:
        for k,v in x_wells.items():
            for well in v:
                wells_x[well].add(k)
    print 'Finished!'

    # find a copy of precalulcated scores, build one if it doesn't exist 
    if not os.path.exists('./pickles'): os.makedirs('./pickles')
    if os.path.isfile('./pickles/val{}.p'.format(w_tot)): (scores,freqs) = pickle.load(open('./pickles/val{}.p'.format(w_tot),'r'))
    else: 
        scores,freqs = precalculate_match_scores(w_tot, match_prior=0.5)
        pickle.dump((scores,freqs),open('./pickles/val{}.p'.format(w_tot),'w'))


    # start iterating through well sets
    tested_pairs = set()
    for w in xrange(w_tot):
        print len(wells_a[w]),len(wells_b[w])
        p_tot = len(wells_a[w])*len(wells_b[w])
        pairs = itertools.product(wells_a[w],wells_b[w])
        for ind,p in enumerate(pairs):
            if not distinct and p[0] == p[1]: continue
            if p in tested_pairs:
                continue
            else:
                # calculate well occurances
                #w_ab = len([a for a in a_wells[p[0]] if a in b_wells[p[1]]]) 
                w_ab = len(set(a_wells[p[0]]).intersection(b_wells[p[1]])) # taps in to O(min(a,b)) for set intersection
                w_a,w_b = len(a_wells[p[0]]),len(b_wells[p[1]])

                # produce well counts (cross-over removed)
                n_ab = w_ab
                n_a,n_b = w_a - n_ab,w_b - n_ab
                n_tot = w_tot

                # create some catches for obviously bad cases
                if n_ab <= 3: continue # since atleast 4 occurances, better be higher than 3 matches
                f_ab = float(w_a*w_b)/(w_tot**2)

                # get score and frequency for parameter set
                score,frequency = scores[n_ab,n_a,n_b],freqs[n_ab,n_a,n_b] 
                 
                if score > threshold:
                    predicted_ab.append(p)
                    predicted_frequency.append(frequency)
                    predicted_score.append(score)

                tested_pairs.add(p)

            if ind%1000000 == 0:
                print 'Finished {}/{}'.format(ind+1,p_tot) 
        
    print ''
     
    return predicted_ab,predicted_frequency,predicted_score # returns edges


# NOTE: This is still a relevant function, but I think it makes more sense to make a new fn for compile
# TODO: actually figure out graph solving feature

class CollectResults:
    def __init__(self):
        self.cells = []
        self.cell_frequencies = []
        self.threshold = []
        self.cell_frequencies_CI = []

    def add_results(self,all_edges,all_freqs,all_scores,id):

        # TODO: add aa/bb solver ( this is non-trivial math )

        # results['cells'] = [e[0] for edge in all_edges for e in edge] # eventually, when we start doing multiclones
        if id == 'ab' or id == 'AB':
            self.cells += [((e[0],),(e[1],)) for e in all_edges]
        elif id == 'aa' or id == 'AA':
            self.cells += [((e[0],e[1]),((),)) for e in all_edges]
        elif id == 'aa' or id == 'AA':
            self.cells += [(((),),(e[0],e[1])) for e in all_edges]

        self.cell_frequencies += [f for f in all_freqs]
        self.threshold += [s for s in all_scores]
        self.cell_frequencies_CI += [(-s/2,s/2) for s in all_scores] # this is pseudo, and totally not right

    def export(self):
        results = {
                    'cells':self.cells,
                    'cell_frequencies':self.cell_frequencies,
                    'threshold':self.threshold,
                    'cell_frequencies_CI':self.cell_frequencies_CI
                  }
        return results


'''
The Meat-and-Potatoes Function
'''
# Verbose: on range 0 to 9
# TODO: verbose levels
def solve(data,pair_threshold = 0.99,verbose=0,real_data=False,all_pairs=True):
    
    if verbose >= 5: silent = False
    else: silent = True
    
    # find uniques
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))

    # sort in ascending order
    sorted(a_uniques,key=int)
    sorted(b_uniques,key=int)

    # counts of each unique in wells
    w_tot = len(data.well_data)
    print 'W_tot:',w_tot

    # create dictionaries for well presence
    a_wells = dict([(a_u,[]) for a_u in a_uniques])
    b_wells = dict([(b_u,[]) for b_u in b_uniques])
    
    if verbose >= 1: print 'Starting variable creation...'
    
    # TODO: Replace with well_dictionary method
    # creates all the necessary images of the data
    for w,well in enumerate(data.well_data):
        for a in well[0]: a_wells[a].append(w) # drop well indices in each unique chains entry
        for b in well[1]: b_wells[b].append(w) # drop well indices in each unique chains entry
        print 'Image generation progress... {}%\r'.format(100*(w+1)/w_tot),

    print ''

    if verbose >= 1: print 'Starting edge detection...'
            
    # Setup threshold values (TODO: better way to distinguish ab,ba from aa,bb
    t = pair_threshold
    t_shared = 1 - (1 - pair_threshold)**2
            
    
    # Find each type of available edge
    ab_edges,ab_freqs,ab_scores = directional_matches(
        a_wells,b_wells,a_uniques,b_uniques,w_tot,threshold=t,silent=silent,distinct=True)
    if verbose >= 2: print 'Finished AB edges!'
        
    if all_pairs:
        aa_edges,aa_freqs,aa_scores = directional_matches(
            a_wells,a_wells,a_uniques,a_uniques,w_tot,threshold=t,silent=silent,distinct=False)
        if verbose >= 2: print 'Finished AA edges!'
        bb_edges,bb_freqs,bb_scores = directional_matches(
            b_wells,b_wells,b_uniques,b_uniques,w_tot,threshold=t,silent=silent,distinct=False)
        if verbose >= 2: print 'Finished BB edges!'
    else:
        aa_edges,aa_freqs,aa_scores = [],[],[]
        bb_edges,bb_freqs,bb_scores = [],[],[]
        
    if verbose >= 1: print 'Finished edge detection, analyzing graph...'
    
    if real_data:
        results_ab = {'edges':ab_edges,'freqs':ab_freqs,'scores':ab_scores}
        results_aa = {'edges':aa_edges,'freqs':aa_freqs,'scores':aa_scores}
        results_bb = {'edges':bb_edges,'freqs':bb_freqs,'scores':bb_scores}
        results = {'AB':results_ab,'AA':results_aa,'BB':results_bb}
        return results

    elif not real_data:
        real_matches = data.metadata['cells']
        # checks to see these actually occur in data
        #potential_ab_matches = possible_matches(real_matches,img_ab,a_uniques,b_uniques)
        
        # solves for true edges
        #all_edges = [ab_edges]#,aa_edges,bb_edges]
        #all_freqs = [ab_freqs]#,aa_freqs,bb_freqs]
        #all_scores = [ab_scores]#,aa_scores,bb_scores]
        #all_uniques = [a_uniques,b_uniques]
       
        #results_ab = {'edges':ab_edges,'freqs':ab_freqs,'scores':ab_scores}
        #results_aa = {'edges':aa_edges,'freqs':aa_freqs,'scores':aa_scores}
        #results_bb = {'edges':bb_edges,'freqs':bb_freqs,'scores':bb_scores}

        #results = {'AB':results_ab,'AA':results_aa,'BB':results_bb}
        
        # deposit results in the collect results function 
        compiler = CollectResults()
        compiler.add_results(ab_edges,ab_freqs,ab_scores,'AB')
        #compiler.add_results(aa_edges,aa_freqs,aa_scores,'AA')
        #compiler.add_results(bb_edges,bb_freqs,bb_scores,'BB')
        
        if verbose >= 1: print 'Finished!'
        
        return compiler.export() 
    
    

    '''
    Scanning for alpha -> beta graph edges
    '''




