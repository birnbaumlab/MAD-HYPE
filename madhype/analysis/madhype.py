

'''
Solver for MAD-HYPE method
'''


# standard libraries
import time
import os
import pickle
from collections import Counter
from operator import mul
from sys import argv
import random
from multiprocessing import cpu_count, Process, Pipe

# nonstandard libraries
from scipy.misc import comb
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# homegrown libraries
from methods import match_probability
from ..defaults import madhype_options as default_options


def solve(data,**kwargs):

    """ Solver for generated and extracted datasets, using MADHYPE algorithm """

    # Determine options
    options = default_options.copy()
    options.update(kwargs)

    # Pull out values into local namespace
    num_wells = data['options']['num_wells']
    cpw = data['options']['cpw']
    num_cores = options['num_cores']
    fdr = options['fdr']
    silent = options['silent']
    

    # Identify unique a,b sequences
    uniques = {
               'A':frozenset().union(*data['well_data']['A']),
               'B':frozenset().union(*data['well_data']['B'])
              }

    # Initialize well distribution
    well_distribution = {'A':{},'B':{}}

    # Create range markers for number of wells 
    pts = [0]+[sum(num_wells[:i+1]) for i in xrange(len(num_wells))]
    
    """ --- Interpret data --- """

    # Identifying well placements for each unique chain
    for label in ('A','B'):
        if not silent: print 'Starting {} chain identification...'.format(label)
        for index in uniques[label]:
            well_distribution[label][index] = \
                    [set([w for w in xrange(pts[i],pts[i+1]) if index in data['well_data'][label][w]])
                        for i in xrange(len(num_wells))]
        if not silent: print 'Finished {} chain identification.'.format(label)

    """ --- Generate data density filter --- """

    # generate filter dictionary
    filt = BuildFilter(well_distribution,num_wells,fdr)

    # initialize list
    results = [] # initialize list

    # REMOVE THIS
    bypass_filter = False 

    # arguments to pass to worker functions
    args = [well_distribution['B'],
            num_wells,
            cpw,
            filt,
            options['prior_alpha'],
            options['prior_match'],
            options['threshold']]

    if num_cores == 0:
        num_cores = cpu_count()
    elif num_cores > cpu_count():
        print 'Number of cores exceeds CPU count, reducing core usage {}->{}...'.format(
                num_cores,cpu_count())
        num_cores = cpu_count

    alpha_dicts = chunkify_dict(well_distribution['A'],num_cores)

    # multithread solver 
    print 'Starting {} processes...'.format(num_cores)
    results = parmap(create_worker(*args),alpha_dicts)
    print 'Finished!'


    print 'Total potential matches: {}'.format(len(results))
    # flatten list
    #results = [entry for subresults in results for entry in subresults]

    # return results
    return results


#------------------------------------------------------------------------------# 
""" Multiprocessing Methods """
#------------------------------------------------------------------------------# 

def _data_intersect(d1,d2,num_wells):
    """ Inputs two lists of sets of indices, formats for match_probability """
    d1,d2 = (set(d1),),(set(d2),)

    w_ij = tuple(len(s1.intersection(s2)) for s1,s2 in zip(d1,d2))
    w_i  = tuple(len(s1) - w for s1,w in zip(d1,w_ij))
    w_j  = tuple(len(s2) - w for s2,w in zip(d2,w_ij))
    w_o  = tuple(w4 - w2 - w3 - w1 for w1,w2,w3,w4 in zip(w_ij,w_i,w_j,num_wells))
    
    return {
            'w_j':w_j,
            'w_ij':w_ij,
            'w_o':w_o,
            'w_tot':num_wells
            }


def spawn(f,index,barcode):
    """ Attaches arbitrary function to pipe """
    def fun(pipe,x):
        pipe.send(f(x,index,barcode))
        pipe.close()
    return fun

def parmap(f,X):
    barcode = np.random.randint(99999) # random barcode to hold temporary files
    pipe=[Pipe() for x in X]
    proc=[Process(target=spawn(f,index+1,barcode),args=(c,x)) 
            for index,(x,(p,c)) in enumerate(zip(X,pipe))]
#    print 'Starting pipes...'
    [p.start() for p in proc]
#    print 'Joining pipes...'
    [p.join() for p in proc]
#    print 'Processing pipes...'
    #recv = [p.recv() for (p,c) in pipe]
    recv = get_results(barcode,len(X))
#    print 'Returning results...'
    return recv

def get_results(barcode,index_range):

    results = []

    for index in xrange(index_range):
        try:
            fname = ".{}_{}.p".format(barcode,index+1)
            results += pickle.load(open(fname, "rb" ))
            #print 'Loaded:',fname
            os.remove(fname)
        except IOError:
            pass

    return results


def create_worker(betas,num_wells,cpw,filt,prior_alpha,prior_match,threshold):
    """ Creates a multiprocessing worker with preset information """

    def worker(alphas,index,barcode):
        """ Worker function, using betas and num_wells """

        total_alphas = len(alphas)
        step_size = max(1,total_alphas/1000) 

        # initialize list
        results = [] # initialize list

        # REMOVE THIS
        bypass_filter = False 

        # Iterate through combinations!
        for i,(a,a_dist) in enumerate(alphas.items()):

            # give heads up on progress
            if i % step_size == 0 and index == 1: 
                #print 'Process {} - {}% complete...'.format(index,round(100*float(i)/total_alphas,1))
                print '{}% complete...'.format(round(100*float(i)/total_alphas,1))

            # apply filter (before itersection,A)
            if filt.check_dist(a_dist) and not bypass_filter: continue

            for j,(b,b_dist) in enumerate(betas.items()):

                # apply filter (before itersection,B)
                if filt.check_dist(b_dist) and not bypass_filter: continue

                # set up input
                pair_data = _data_intersect(a_dist,b_dist,num_wells)

                # apply filter (on intersection)
                if filt.check_tuple(pair_data['w_ij']) and not bypass_filter: continue

                # add keys to dictionary where needed
                pair_data['alpha'] = prior_alpha
                pair_data['cpw'] = cpw
                pair_data['label'] = ((a,),(b,))

                # calculate match probability
                p,f = match_probability(pair_data,prior_match)

                if p > threshold:
                    if filt.check_tuple(pair_data['w_ij']): continue
                    results.append((((a,),(b,)),p,f[0]))

#        print '\nCore {} finished!'.format(index)

#        if index == 1:
#            print '\nWaiting on other cores to finish...!'

        pickle.dump(results, open( ".{}_{}.p".format(barcode,index), "wb" ))
#        print 'Pickled to: .{}_{}.p'.format(barcode,index)

        #return []#results

    # returns worker function
    return worker 

def chunkify_iterable(iterable,n):
    """ Breaks iterable into n approx. equal chunks """ 
    return list(iterable[i::n] for i in range(n))

def chunkify_dict(mydict,n):
    """ Breaks iterable into n approx. equal chunks """ 
    subdict = [[] for i in range(n)]
    [subdict[i%n].append(item) for i,item in enumerate(mydict.items())]
    return list(dict(s) for s in subdict)

#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _data_intersect(d1,d2,num_wells):
    """ Inputs two lists of sets of indices, formats for match_probability """
    w_ij = tuple(len(s1.intersection(s2)) for s1,s2 in zip(d1,d2))
    w_i  = tuple(len(s1) - w for s1,w in zip(d1,w_ij))
    w_j  = tuple(len(s2) - w for s2,w in zip(d2,w_ij))
    w_o  = tuple(w4 - w2 - w3 - w1 for w1,w2,w3,w4 in zip(w_ij,w_i,w_j,num_wells))
    return {
            'w_i':w_i,
            'w_j':w_j,
            'w_ij':w_ij,
            'w_o':w_o,
            'w_tot':num_wells
           }

#------------------------------------------------------------------------------# 

class BuildFilter(object):

    def __init__(self,well_distribution,num_wells,fdr):

        """ Initializes filter using well distributions and options """

        w_tot = sum(num_wells)
        w_dict = Counter([tuple(len(i) for i in w) for w in well_distribution['A'].values()] +
                         [tuple(len(i) for i in w) for w in well_distribution['B'].values()])
        self.w_filter = dict([(k,1./
                         ((v-1)/reduce(mul,[comb(W,w) for W,w in zip(num_wells,k)]) + 1))
                     for k,v in w_dict.items()])
        self.fdr = fdr

    def check_tuple(self,value):

        """ Checks whether a certain value passes filter """

        try:
            if self.w_filter[value] < 1. - self.fdr: return True
        except KeyError:
            pass

        return False

    def check_dist(self,my_iter):

        """ Checks whether a certain value passes filter """

        value = tuple(len(s) for s in my_iter)

        try:
            if self.w_filter[value] < 1. - self.fdr: return True
        except KeyError:
            pass

        return False

#------------------------------------------------------------------------------# 

