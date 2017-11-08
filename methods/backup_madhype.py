
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
import multiprocessing
import threading
import collections

# nonstandard libraries
import matplotlib.pyplot as plt # temporary, delete me later
from datetime import datetime
import numpy as np
import scipy.misc



'''
Specialty Calculators
'''
def chunkIt(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

# N choose K
def nCk(n, r):
    return scipy.misc.comb(n,r)

'''
Important Functions
'''

def madhype_thread(args):
    startTime = datetime.now() 
    os.system(os.getcwd() + '/solver/test_backup {} {} {}'.format(*args))
    print 'Thread-{} took {} seconds.\n'.format(args[-1],datetime.now()-startTime)


def determine_core_usage(cores): 
    if cores == 'max':
        print 'Declared max CPU core usage ({})...'.format(multiprocessing.cpu_count())
        return multiprocessing.cpu_count()
    elif multiprocessing.cpu_count() < cores:
        print 'Number of declared cores larger than detected, reducing usage {} -> {}'.format(
                cores,multiprocessing.cpu_count())
        return multiprocessing.cpu_count()
    else:
        return cores

def multithread_madhype(cores,data,args):
    # start a timer
    startTime = datetime.now()

    a_uniques,b_uniques,a_wells,b_wells = data # unpackage data 

    # try to reduce dependencies on data processing here
    for i,a_u in enumerate(chunkIt(a_uniques,cores)):
        with open('./solver/chain_data_a_{}.txt'.format(i+1),'w') as f:
            for w in a_u: f.write('{}'.format(str(a_wells[w]))[1:-1]+'\n')
        with open('./solver/uniques_a_{}.txt'.format(i+1),'w') as f:
            for w in a_u: f.write('{}\n'.format(w))
    with open('./solver/chain_data_b.txt','w') as f:
        for w in b_uniques: f.write('{}'.format(str(b_wells[w]))[1:-1]+'\n')
    with open('./solver/uniques_b.txt','w') as f:
        for w in b_uniques: f.write('{}\n'.format(w))

    # create parameter sets for each core
    arg_lists = [args+[str(i+1)] for i in xrange(cores)] 
    
    # create threads and start the engines
    pool = multiprocessing.Pool(processes = cores)
    pool.map(madhype_thread, arg_lists)

    # let us know how long everything took
    print 'Multithreaded C++ full implementation took {} seconds.\n'.format(datetime.now()-startTime)


"""
Compiles the results from multiple core outputs, returns a dictionary of results
"""
def collect_results(core_count):
    lines = []
    # get all the results into one list
    for i in xrange(core_count):
        with open('results_{}.txt'.format(i+1),'r') as f:
            lines += [l.strip().split('\t') for l in f.readlines()]
    # save in clean form 
    return lines

def export_results(ab_lines,aa_lines,bb_lines):
    # initialize list
    results = []
    # neatly packages all the data from AB,AA,BB
    for lines in [ab_lines,aa_lines,bb_lines]:
        edges = [(int(l[2]),int(l[3])) for l in lines]
        freqs = [float(l[1]) for l in lines]
        scores = [float(l[0]) for l in lines]
        results.append({'edges':edges,'freqs':freqs,'scores':scores})
    # return results
    return {'AB':results[0],'AA':results[1],'BB':results[2]}
    
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
# Adjustment: refers to the prior distribution adjustment account for repertoire density error
def solve(data,pair_threshold = 0.99,verbose=0,real_data=False,all_pairs=False,repertoire_adjustment=False,cores='max'):
    
    if verbose >= 5: silent = False
    else: silent = True

    w_tot = len(data.well_data)

    # Find uniques
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))

    # Generate reverse dictionary
    a_wells = dict([(a_u,[]) for a_u in a_uniques])
    b_wells = dict([(b_u,[]) for b_u in b_uniques])
    
    if verbose >= 1: print 'Starting reverse dictionary creation...'
    
    # creates all the necessary images of the data
    for w,well in enumerate(data.well_data):
        for a in well[0]: a_wells[a].append(w) # drop well indices in each unique chains entry
        for b in well[1]: b_wells[b].append(w) # drop well indices in each unique chains entry
        print 'Reverse dictionary progress... {}%\r'.format(100*(w+1)/w_tot),
    print ''

    # detect the appropriate amount of core usage
    core_count = determine_core_usage(cores) # calls simple function to figure out core counts

    # C++ embedding  
    new_thresh = -math.log10((1./pair_threshold - 1)) if pair_threshold>0 else float('-inf')
    args = [str(w_tot),str(new_thresh)]
    
    startTime = datetime.now()

    ### MAD-HYPE multithreaded C++ implementation ###
    # work horse (which loads everything and sets up threads)

    # First, do the necessary A/B pairs:
    passed_data = [a_uniques,b_uniques,a_wells,b_wells]
    multithread_madhype(core_count,passed_data,args)
    lines_ab = collect_results(core_count)

    if repertoire_adjustment:
        # make histogram dictionaries for multiplicity adjustments
        counter_a = collections.Counter([len(wc) for wc in a_wells.values()])
        counter_b = collections.Counter([len(wc) for wc in b_wells.values()])
        #counter_a = dict([(k,v*math.log10((nCk(w_tot,k)-1)/nCk(w_tot,k))) for k,v in counter_a.items()])
        #counter_b = dict([(k,v*math.log10((nCk(w_tot,k)-1)/nCk(w_tot,k))) for k,v in counter_b.items()])
        # Old set, for later investigation
        counter_a_old = dict([(k,math.log10(v)) for k,v in counter_a.items()]) 
        counter_b_old = dict([(k,math.log10(v)) for k,v in counter_b.items()]) 
        counter_a = dict([(k,math.log10(1./(v*(1./nCk(w_tot,k)) + (1.-(1./nCk(w_tot,k))))))
            for k,v in counter_a.items()]) 
        counter_b = dict([(k,math.log10(1./(v*(1./nCk(w_tot,k)) + (1.-(1./nCk(w_tot,k)))))) 
            for k,v in counter_b.items()]) 

        for k in counter_a.keys():
            print '{}: {} vs. {}'.format(k,counter_a[k],10.**counter_a_old[k])

        

        # make adjustment
        lines_ab = [[float(l[0]) + 1.0*(counter_a[len(a_wells[int(l[2])])] + counter_b[len(b_wells[int(l[3])])])] + l[1:] for l in lines_ab]

    # Next, consider doing A/A or B/B pairs
    if all_pairs:
        passed_data = [a_uniques,a_uniques,a_wells,b_wells]
        multithread_madhype(core_count,passed_data,args)
        lines_aa = collect_results(core_count)

        passed_data = [a_uniques,b_uniques,a_wells,b_wells]
        multithread_madhype(core_count,passed_data,args)
        lines_bb = collect_results(core_count)
    else:
        lines_aa = []
        lines_bb = []

    # real data post-processing section
    if real_data:
        # collect results in each core output file
        results = export_results(lines_ab,lines_aa,lines_bb)
        return results

    # non-so-real data post-processing
    elif not real_data:
        # recalls real matches
        real_matches = data.metadata['cells']

        # parse lines in results file
        ab_edges, ab_freqs, ab_scores = [],[],[]
        for line in lines_ab:
            ab_edges.append((int(line[2]),int(line[3])))
            ab_freqs.append(float(line[1]))
            ab_scores.append(float(line[0]))
        
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




