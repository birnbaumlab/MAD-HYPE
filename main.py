
"""
Tests out the new and improved variable solver
"""

# standard libraries
import time
from collections import Counter
from operator import mul

# nonstandard libraries
from scipy.misc import comb

# homegrown libraries
from solver.methods import match_probability
from sequence_generator import *
from post_processing import visualize_results

#------------------------------------------------------------------------------# 

def main(*args,**kwargs):

    options = {
              #'num_wells':(24,24,24,24),
              #'cpw':(1,10,100,1000),
              'num_wells':(24,),
              'cpw':(500,),
              'num_cells':1000,
              'cell_freq_distro':'power-law',
              'cell_freq_constant':1,
              'chain_misplacement_prob':0,
              'chain_deletion_prob':0.1
              }

    # Update settings
    for arg in args: options.update(arg)
    options.update(kwargs)

    # Pull out values into local namespace
    num_wells = options['num_wells']
    cpw = options['cpw']
    threshold = 0.1
    fdr = 0.01
    prior_alpha = 0

    # Generate datasets
    sg = DataGenerator(options)
    sg.generate_cells()
    data = sg.generate_data()

    # Identify unique a,b sequences
    uniques = {
               'A':frozenset().union(*data['well_data']['A']),
               'B':frozenset().union(*data['well_data']['B'])
              }

    # Initialize well distribution
    well_distribution = {'A':{},'B':{}}

    # Create range markers for number of wells 
    pts = [0]+[sum(num_wells[:i+1]) for i in xrange(len(num_wells))]
    
    # Identifying well placements for each unique chain
    for label in ('A','B'):
        print 'Starting {} chain identification...'.format(label)
        for index in uniques[label]:
            well_distribution[label][index] = \
                    [set([w for w in xrange(pts[i],pts[i+1]) if index in data['well_data'][label][w]])
                        for i in xrange(len(num_wells))]
        print 'Finished {} chain identification.'.format(label)
        
    # TODO: include frequency filter 
    # TODO: standardize threshold

    # generate filter dictionary
    w_tot = sum(options['num_wells'])
    w_dict = Counter([tuple(len(i) for i in w) for w in well_distribution['A'].values()] +
                     [tuple(len(i) for i in w) for w in well_distribution['B'].values()])
    w_filter = dict([(k,1./
                     ((v-1)/reduce(mul,[comb(W,w) for W,w in zip(options['num_wells'],k)]) + 1))
                     for k,v in w_dict.items()])

    '''#
    print 'Filter:'
    for k,v in w_filter.items():
        print '{}: p = {}'.format(k,v)
    #'''

    #'''#
    print 'Cells:'
    for i in xrange(options['num_cells']):
        c = ((i,),(i,))
        print '{}: p = {}'.format(c,data['cells'][c])
    #'''

    # initialize list
    results = [] # initialize list
    freqs = [] # initialize list
    pair_datas = [] # initialize list

    # Iterate through combinations!
    for i,a in enumerate(uniques['A']):
        if i % 10 == 0: print 'Starting A-chain {}...'.format(i)
        for j,b in enumerate(uniques['B']):
            # set up input
            pair_data = _data_intersect(
                    well_distribution['A'][a],well_distribution['B'][b],num_wells)

            # apply filter
            try:
                if w_filter[pair_data['w_ij']] < 1. - fdr:
                    continue
            except KeyError:
                pass

            pair_data['alpha'] = prior_alpha
            pair_data['cpw'] = cpw
            pair_data['label'] = ((a,),(b,))
            # calculate match probability
            p,f = match_probability(pair_data)
            
            if p > threshold:
                results.append((((a,),(b,)),p))
                pair_datas.append((pair_data,p))
                freqs.append((f,p))

    results.sort(key=lambda x: -x[1])
    freqs.sort(key=lambda x: -x[1])
    pair_datas.sort(key=lambda x: -x[1])

    for i,r,f,p in zip(xrange(len(results[:1000])),results,freqs,pair_datas):
        if r[0][0] == r[0][1]: continue
        else: print 'WRONG:',(r[0][0][0],r[0][1][0]),r,
        print [len(a) for a in (well_distribution['A'][r[0][0][0]])],
        print [len(b) for b in (well_distribution['B'][r[0][1][0]])],p[0]
        #'''#
        print '> ',[(a) for a in (well_distribution['A'][r[0][0][0]])]
        print '> ',[(b) for b in (well_distribution['B'][r[0][1][0]])]
        print '> ',f[0][0]
        print '> ',f[0][1]
        #'''#

    results.sort(key=lambda x: -x[1])
    freqs.sort(key=lambda x: -x[1])

    visualize_results(results,data)

    '''
    print well_distribution['A'][0]
    print well_distribution['B'][0]
    print data['cells'][((0,),(0,))]
    '''


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

if __name__ == "__main__":
    main()


