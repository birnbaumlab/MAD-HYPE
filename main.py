
"""
Tests out the new and improved variable solver
"""

# standard libraries
import time

# nonstandard libraries

# homegrown libraries
from solver.methods import match_probability
from sequence_generator import *

#------------------------------------------------------------------------------# 

def main(*args,**kwargs):

    options = {
              'num_wells':[96],
              'cpw':[100],
              'num_cells':25,
              'cell_freq_distro':'power-law',
              'cell_freq_constant':-1,
              'chain_misplacement_prob':0,
              'chain_deletion_prob':0.1
              }

    # Update settings
    for arg in args: options.update(arg)
    options.update(kwargs)

    # Pull out values into local namespace
    num_wells = options['num_wells']
    cpw = options['cpw']
    threshold = 0.01
    prior_alpha = 1

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
    pts = [0]+[sum(num_wells[:i]) for i in xrange(len(num_wells))]
    
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

    results = [] # initialize list

    # Iterate through combinations!
    for i,a in enumerate(uniques['A']):
        if i % 1 == 0: print 'Starting A-chain {}...'.format(i)
        for j,b in enumerate(uniques['B']):
            # set up input
            pair_data = _data_intersect(
                    well_distribution['A'][a],well_distribution['B'][b],cpw)
            pair_data['alpha'] = prior_alpha
            pair_data['cpw'] = cpw
            
            # calculate match probability
            p = match_probability(pair_data)
            
            if p > 5:
                results.append((((a,),(b,)),p))

    results.sort(key=lambda x: x[1])

    print results

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
    w_o  = tuple(w4 - w2 - w3 + w1 for w1,w2,w3,w4 in zip(w_ij,w_i,w_j,num_wells))
    print d1,d2,num_wells
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


