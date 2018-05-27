
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
    """ This function is what is executed each process """ 
    '''
    Arguments:
        > args: list or tuple, consisting of (w_tot, threshold, process #)
    '''
    startTime = datetime.now() # start a timeer for processing
    os.system('python ' + os.getcwd() + '/solver/test.py ' + ' '.join(args))
    print 'Process-{} took {} seconds.\n'.format(args[-1],datetime.now()-startTime) # update on processing time


def determine_core_usage(cores): 
    """ This function checks that declared core usage is approriate """
    '''
    Arguments:
        > cores: string ('max') or int, number of processes that will be initialized
                 not to be greater than the number of cores on CPU
    '''
    if type(cores) == str:
        if cores.lower() == 'max': # check lower case name
            print 'Declared max CPU core usage ({})...'.format(multiprocessing.cpu_count())
            return multiprocessing.cpu_count()
        else:
            raise KeyError('cores string is not recognized ({})'.format(cores))
    elif not type(cores) == int: 
        raise TypeError('cores is of type <{}>, not type <{}>'.format(type(cores),'int'))
    elif multiprocessing.cpu_count() < cores:
        print 'Number of declared cores larger than detected, reducing usage {} -> {}'.format(
                cores,multiprocessing.cpu_count())
        return multiprocessing.cpu_count()
    else:
        return cores

def try_chain_additions_multithread(init_chainsets, a_uniques, b_uniques, cores, args):
    print args

    # Map alpha/beta chain ids to indices
    alpha_id_to_idx = {a: i for i,a in enumerate(a_uniques)}
    beta_id_to_idx = {b: i for i,b in enumerate(b_uniques)}

    # Output process-specific info to data files
    chunk_size = int(math.ceil(len(init_chainsets) / float(cores)))
    print len(init_chainsets), chunk_size
    for i in range(cores):
        f = open('./solver/initial_{}.txt'.format(i+1), 'w')
        for w in init_chainsets[i*chunk_size:(i+1)*chunk_size]:
            a1 = alpha_id_to_idx[w[0][0]] if len(w[0])>=1 else -1
            a2 = alpha_id_to_idx[w[0][1]] if len(w[0])>=2 else -1
            b1 = beta_id_to_idx[w[1][0]] if len(w[1])>=1 else -1
            b2 = beta_id_to_idx[w[1][1]] if len(w[1])>=2 else -1
            f.write('{},{},{},{}\n'.format(a1,a2,b1,b2))
        f.close()

    # Prep process-specific arguments
    args_list = [args+[str(i+1)] for i in xrange(cores)]

    # Start processes
    pool = multiprocessing.Pool(processes = cores)
    pool.map(madhype_thread, args_list)

    # Collect and parse results
    lines = collect_results(cores)

    cells, scores_dict, freqs_dict = set(), {}, {}
    for line in lines:
        a1,a2,b1,b2 = [int(v) for v in line[2:]]
        a = () if a1==-1 else (a1,) if a2==-1 else (a1,a2)
        b = () if b1==-1 else (b1,) if b2==-1 else (b1,b2)
        c = (a,b)
        cells.add(c)
        scores_dict[c] = float(line[0])
        freqs_dict[c] = float(line[1])
    return cells, scores_dict, freqs_dict

           

def multithread_madhype(cores,data,args_dict):
   # start a timer
    startTime = datetime.now()

    a_uniques,b_uniques,a_wells,b_wells = data # unpackage data 

    # Output general data used by all processes
    with open('./solver/chain_data_a.txt','w') as f:
        for w in a_uniques: f.write('{}'.format(str(a_wells[w]))[1:-1]+'\n')
    with open('./solver/uniques_a.txt','w') as f:
        for w in a_uniques: f.write('{}\n'.format(w))
    with open('./solver/chain_data_b.txt','w') as f:
        for w in b_uniques: f.write('{}'.format(str(b_wells[w]))[1:-1]+'\n')
    with open('./solver/uniques_b.txt','w') as f:
        for w in b_uniques: f.write('{}\n'.format(w))

    # Extract arguments
    w_tot = args_dict['w_tot']
    base_thresh = args_dict['pair_threshold']

    # Run first pass (match alphas to betas)
    pass1_prior = 1.5/math.sqrt(len(a_uniques)*len(b_uniques))
    pass1_thresh = base_thresh - (math.log10(pass1_prior) - math.log10(1-pass1_prior))
    pass1_args = [str(w_tot), str(pass1_thresh), '0', '1']
    pass1_init = [((a,),()) for a in a_uniques]
    p1_cells, p1_scores_dict, p1_freqs_dict = try_chain_additions_multithread(pass1_init, a_uniques, b_uniques, cores, pass1_args)

    cells = set(p1_cells)
    scores_dict = {c: s-pass1_thresh+base_thresh for c,s in p1_scores_dict.iteritems()}
    freqs_dict = p1_freqs_dict.copy()
    print "MADHYPE-PYTHON PASS 1: Found {} cells".format(len(p1_cells))

    # Run second pass (add alpha/beta chains to each chain pair to make duals)
#    pass2_prior = .13/math.sqrt(len(a_uniques)*len(b_uniques))
#    pass2_thresh = base_thresh - (math.log10(pass2_prior) - math.log10(1-pass2_prior))
#    pass2_args = [str(w_tot), str(pass2_thresh), '1', '1']
#    pass2_init = list(cells)
#    p2_cells, p2_scores, p2_freqs = try_chain_additions_multithread(pass2_init, a_uniques, b_uniques, cores, pass2_args)

#    for alist,blist in p2_cells:
#        for c in [((a,),(b,)) for a in alist for b in blist]:  cells.discard(c)
#    print "MADHYPE-PYTHON PASS 2: Discarded {} cells, added {} cells".format(len(set(p1_cells) - cells), len(p2_cells))
#    cells |= p2_cells
#    scores_dict.update({c:s-pass2_thresh+base_thresh for c,s in p2_scores.iteritems()})
#    freqs_dict.update(p2_freqs)

    
     


    

    # let us know how long everything took
    print 'Multithreaded C++ full implementation took {} seconds.\n'.format(datetime.now()-startTime)

    return list(cells), scores_dict, freqs_dict


"""
Compiles the results from multiple core outputs, returns a dictionary of results
"""
def collect_results(core_count,skip_run=False):
    lines = []
    if not skip_run: dirname = 'solver/'
    else: dirname = 'saved_data/'+skip_run+'/'
    # get all the results into one list
    for i in xrange(core_count):
        with open(dirname + 'results_{}.txt'.format(i+1),'r') as f:
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

    def add_results(self,all_edges,all_freqs,all_scores):

        # TODO: add aa/bb solver ( this is non-trivial math )

        # results['cells'] = [e[0] for edge in all_edges for e in edge] # eventually, when we start doing multiclones
        #if id == 'ab' or id == 'AB':
        #    self.cells += [((e[0],),(e[1],)) for e in all_edges]
        #elif id == 'aa' or id == 'AA':
        #    self.cells += [((e[0],e[1]),((),)) for e in all_edges]
        #elif id == 'aa' or id == 'AA':
        #    self.cells += [(((),),(e[0],e[1])) for e in all_edges]
        print all_edges[0]
        self.cells += all_edges

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
# Input parameters:
#   > skip_run: if False, everything will run as normal; if set to string, will dig for a folder in testing for data
def solve(data,custom_params={}):

    # TODO: I would eventually like to remove this declaration, and have the parameter dictionary passed always contain the defaults instead
    # In the interim, this will do
    # default parameters
    default_params = {
            'max_threshold':None,
            'pair_threshold':2,
            'verbose':0,
            'real_data':False, 
            'all_pairs':False,
            'repertoire_adjustment':False,
            'cores':'max', 
            'skip_run':False
            }

    params = Parameters(default_params)

    # included default parameters (TODO: again, I would like to eventually stay away from local namespace pollution, but until we switch all code this will do)
    params.update(custom_params)
    
    if params.verbose >= 5: silent = False
    else: silent = True

    w_tot = len(data.well_data)

    # Find uniques
    a_uniques = list(set([a for well in data.well_data for a in well[0]]))
    b_uniques = list(set([b for well in data.well_data for b in well[1]]))

    # Generate reverse dictionary
    a_wells = dict([(a_u,[]) for a_u in a_uniques])
    b_wells = dict([(b_u,[]) for b_u in b_uniques])
    
    if params.verbose >= 1: print 'Starting reverse dictionary creation...'
    
    # creates all the necessary images of the data
    for w,well in enumerate(data.well_data):
        for a in well[0]: a_wells[a].append(w) # drop well indices in each unique chains entry
        for b in well[1]: b_wells[b].append(w) # drop well indices in each unique chains entry
        print 'Reverse dictionary progress... {}%\r'.format(100*(w+1)/w_tot),
    print ''

    # detect the appropriate amount of core usage
    core_count = determine_core_usage(params.cores) # calls simple function to figure out core counts

    # C++ embedding  
    # old version where we were not using ratios
    #new_thresh = -math.log10(1./pair_threshold - 1) if pair_threshold>0 else float('-inf')
    #args = [str(w_tot),str(new_thresh)]
    args = {'w_tot': w_tot, 'pair_threshold': params.pair_threshold}
    
    startTime = datetime.now()

    ### MAD-HYPE multithreaded C++ implementation ###
    # work horse (which loads everything and sets up threads)

    # First, do the necessary A/B pairs:
    passed_data = [a_uniques,b_uniques,a_wells,b_wells]
    if not params.skip_run:
        cells, scores_dict, freqs_dict = multithread_madhype(core_count,passed_data,args)
    #lines_ab = collect_results(core_count,params.skip_run)

    # TODO: Modify this to work with new results format
    if params.repertoire_adjustment:
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
    # TODO: This is not compatible with the new multithread_madhype()
    if params.all_pairs:
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
    if params.real_data:
        # collect results in each core output file
        results = export_results(lines_ab,lines_aa,lines_bb)
        return results

    # non-so-real data post-processing
    elif not params.real_data:
        # recalls real matches
        real_matches = data.metadata['cells']

        scores = [scores_dict[c] for c in cells]
        freqs = [freqs_dict[c] for c in cells]

        # deposit results in the collect results function 
        compiler = CollectResults()
        compiler.add_results(cells, freqs, scores)
        #compiler.add_results(ab_edges,ab_freqs,ab_scores,'AB')
        #compiler.add_results(aa_edges,aa_freqs,aa_scores,'AA')
        #compiler.add_results(bb_edges,bb_freqs,bb_scores,'BB')
        
        if params.verbose >= 1: print 'Finished!'
        
        return compiler.export() 
    
    


class Parameters:
    """ Creates a default set of analytical parameters, replaces with any input dictionary """ 
    # use a dictionary to update class attributes
    def update(self,params={}):
        # makes params dictionary onto class attributes
        print params
        for key, value in params.items():
            setattr(self, key, value)
        # can apply checks afterward

    # initialize class with default dictionary
    def __init__(self,default_params,custom_params={}):
        self.update(default_params)
        self.model_parameters = default_params.keys()
        self.update(custom_params)

    def dict(self):
        return self.__dict__



