
"""
Tests out the new and improved variable solver
"""

# standard libraries
import time
from collections import Counter
from operator import mul

# nonstandard libraries
from scipy.misc import comb
import matplotlib.pyplot as plt
from matplotlib import cm

# homegrown libraries
from solver.methods import match_probability
from sequence_generator import *
from post_processing import visualize_results
from post_processing import analyze_results

#------------------------------------------------------------------------------# 

def test(*args,**kwargs):

    options = {
              # experimental design
              'num_wells':(24,),
              'cpw':(10,),
              # simulated repertoire
              'num_cells':1000,
              'cell_freq_distro':'power-law',
              'cell_freq_constant':1.0,
              'chain_misplacement_prob':0, # TODO: add functionality
              'chain_deletion_prob':0.1,
              # analysis constants
              'threshold':0.1, # minimum ratio accepted by match_probability
              'fdr':0.01, # acceptable fdr (cuts off matches, sets filter)
              'prior_alpha':1.0, # prior for clonal frequency
              'prior_match':1.0, # prior for clonal match ( <= 1.0 )
              # visual cues
              'silent':False,
              'visual':False
              }

    # Update settings
    # TODO: add warning for non-key add
    for arg in args: options.update(arg)
    options.update(kwargs)

    # Pull out values into local namespace
    num_wells = options['num_wells']
    cpw = options['cpw']
    silent = options['silent']

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
    filt = BuildFilter(well_distribution,options)

    '''#
    print 'Filter:'
    for k,v in w_filter.items():
        print '{}: p = {}'.format(k,v)
    #'''

    '''#
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

        # give heads up on progress
        if i % 10 == 0: 
            if not silent: print 'Starting A-chain {}...'.format(i)

        # apply filter (before itersection,A)
        if filt.check_dist(well_distribution['A'][a]): continue

        for j,b in enumerate(uniques['B']):

            # apply filter (before itersection,B)
            if filt.check_dist(well_distribution['B'][b]): continue

            # set up input
            pair_data = _data_intersect(
                    well_distribution['A'][a],well_distribution['B'][b],num_wells)

            # apply filter (on intersection)
            if filt.check_tuple(pair_data['w_ij']): continue

            # add keys to dictionary where needed
            pair_data['alpha'] = options['prior_alpha']
            pair_data['cpw'] = cpw
            pair_data['label'] = ((a,),(b,))

            # calculate match probability
            p,f = match_probability(pair_data,options['prior_match'])
            
            if p > options['threshold']:
                results.append((((a,),(b,)),p))

    '''#
    for i,r,f,p in zip(xrange(len(results[:1000])),results,freqs,pair_datas):
        if r[0][0] == r[0][1]: continue
        else: print 'WRONG:',(r[0][0][0],r[0][1][0]),r,
        print [len(a) for a in (well_distribution['A'][r[0][0][0]])],
        print [len(b) for b in (well_distribution['B'][r[0][1][0]])],p[0]
        print '> ',[(a) for a in (well_distribution['A'][r[0][0][0]])]
        print '> ',[(b) for b in (well_distribution['B'][r[0][1][0]])]
        print '> ',f[0][0]
        print '> ',f[0][1]
    #'''


    # gather results
    compiled_results = analyze_results(results,data,options)

    # visualize results if requested
    if options['visual']:
        visualize_results(results,data,options)

    # return compiled results
    return compiled_results


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

    def __init__(self,well_distribution,options):

        """ Initializes filter using well distributions and options """

        w_tot = sum(options['num_wells'])
        w_dict = Counter([tuple(len(i) for i in w) for w in well_distribution['A'].values()] +
                         [tuple(len(i) for i in w) for w in well_distribution['B'].values()])
        self.w_filter = dict([(k,1./
                         ((v-1)/reduce(mul,[comb(W,w) for W,w in zip(options['num_wells'],k)]) + 1))
                     for k,v in w_dict.items()])
        self.fdr = options['fdr']

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

if __name__ == "__main__":

    mode = 1

    if mode == 1:
        
        options = {
                'num_cells':1000,
                'num_wells':(24,),
                'cpw':(10,200),
                'seed':1,
                # visual cues
                'silent':True,
                'visual':True
                }


        results = test(options)

    if mode == 2:
    
        # system parameters
        repeats = 50
        cpw_range = np.logspace(0,4,41,dtype=int)
        num_wells_range = [12,24,36,48]

        # simulation parameters
        options = {
                'num_cells':1000,
                'num_wells':(24,),
                'seed':1,
                'silent':True
                }

        id_map = np.zeros((len(cpw_range),options['num_cells']))

        for ii,num_wells in enumerate(num_wells_range): 
            for i,cpw in enumerate(cpw_range):
                for seed in xrange(repeats):
                    results = test(options,cpw=(cpw,),seed=seed,num_wells=(num_wells,))
                    id_map[i,:] += results['pattern']
                print 'Finished cpw = {}!'.format(cpw)

            # display graph
            c_labels = [v for v in [1,10,100,1000,10000] if v <= max(cpw_range)]
            c_inds = [min(range(len(cpw_range)), key=lambda i: abs(cpw_range[i]-v)) for v in c_labels]

            fig, ax = plt.subplots()

            cax = ax.imshow(id_map,interpolation='nearest')
            ax.set_aspect(aspect='auto')

            ax.set_yticks(c_inds)
            ax.set_yticklabels(c_labels)

            plt.title('Identification of clones with {} wells'.format(num_wells))
            plt.xlabel('Clonal Index')
            plt.ylabel('Cells per well (#)')
            cbar = fig.colorbar(cax, ticks=[0,repeats])
            cbar.ax.set_yticklabels(['0%','100%'])  # vertically oriented colorbar
            
            fig.savefig('{}_wells.eps'.format(num_wells),format='eps',dpi=1000)

        # show, allow user to exit
        plt.show(block=False)
        raw_input('Press enter to close...')
        plt.close()




