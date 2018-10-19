
"""

Authors: Patrick V. Holec, Joseph Berleant
Created: 4 December 2017

"""

"""

This is a new sequence generator that is compatible with variable cpw

"""

# standard libraries
import itertools

# nonstandard libraries
import numpy as np
from numpy.random import binomial

# homegrown libraries
from ..defaults import general_options as default_options

#------------------------------------------------------------------------------# 

class DataGenerator(object):

    """
    Methods:
        __init__(self,**kwargs)
            input: kwargs
            output: (none)
        update(self,**kwargs)
            input: kwargs
            output: (none)
        generate_data(self, cells, cell_frequencies, seed [op.])
            input: list of cells, list of cell frequencies, int
            output: data dict 
    """ 


    def __init__(self, **kwargs):
        """ Initialize sequence generator object """
        ## Set parameter values
        self.options = default_options.copy()
        
        ## update options 
        self.update(**kwargs)

    def update(self,**kwargs):
        """ Update options using *args and **kwargs """ 
        for k in self.options:
            if k in kwargs:  self.options[k] = kwargs[k]

    def generate_data(self, cells, cell_frequencies, seed = None):
        """ Generate simulated sequencing data from the given list of cells """

        # set random seed
        if seed is not None:
          np.random.seed(seed)

        # transfer options to local namespace
        num_cells = len(cells)
        num_wells = self.options['num_wells']
        cpw = self.options['cpw']
        chain_deletion_prob = self.options['chain_deletion_prob']
        chain_misplacement_prob = self.options['chain_misplacement_prob']

        # Check that num_wells are cpw are consistent with each other
        assert (isinstance(num_wells, int) and isinstance(cpw, int)) or \
                (len(num_wells) == len(cpw)), \
                "num_wells and cpw must both be ints or iterables of equal length"

        # interpret user input for cpw distribution
        if isinstance(cpw,int) and isinstance(num_wells,int):
            self.cpw = [cpw for _ in xrange(num_wells)]
        elif isinstance(cpw,(tuple,list)) and isinstance(num_wells,(tuple,list)):
            assert len(cpw) == len(num_wells),"cpw and num_wells not same length"
            self.cpw = [c for w,c in zip(num_wells,cpw) for _ in xrange(w)]

        # start generating well data
        self.well_data = {'A':[],'B':[]}

        # aggregate misplaced chains
        misplaced_a = []
        misplaced_b = []

        for well_ind,cpw in enumerate(self.cpw):
            # select cell indices for each well
            indices = np.random.choice(num_cells,size=(cpw,),p=cell_frequencies)
            # seperate chains in to a,b
            a = list(itertools.chain(*[cells[i][0] for i in indices]))
            b = list(itertools.chain(*[cells[i][1] for i in indices]))
            # shuffle listed chains
            np.random.shuffle(a)
            np.random.shuffle(b)
            # truncation length, then remove duplicates
            a_ind,b_ind = binomial(len(a),chain_misplacement_prob),binomial(len(b),chain_misplacement_prob)
            a,b = (a[a_ind:]),(b[b_ind:])
            # misplaced chains, aggregated for later distribution
            misplaced_a += a[:a_ind]
            misplaced_b += b[:b_ind]
            # store in wells
            self.well_data['A'].append(a)
            self.well_data['B'].append(b)

        # indices to place misplaced chains
        a_ind = np.random.randint(0,len(self.cpw),(len(misplaced_a),))
        b_ind = np.random.randint(0,len(self.cpw),(len(misplaced_b),))

        for i in xrange(len(self.cpw)):
            self.well_data['A'][i] += [a for a,ind in zip(misplaced_a,a_ind) if ind == i] 
            self.well_data['B'][i] += [b for b,ind in zip(misplaced_b,b_ind) if ind == i] 

        self.well_data['A'] = [list(set(wd[binomial(len(wd),chain_deletion_prob):])) for wd in self.well_data['A']]
        self.well_data['B'] = [list(set(wd[binomial(len(wd),chain_deletion_prob):])) for wd in self.well_data['B']]

        # compile useful information
        data = {
                'well_data':self.well_data,
                'cells':[(((a,),(b,)),f) for c,f in zip(cells,cell_frequencies)
                    for a in c[0] for b in c[1]],
                'complete_cells':[(c,f) for c,f in zip(cells,cell_frequencies)],
                'options':self.options
               }

        return data # return results

def generate_data(cells, freqs, seed = None, **kwargs):
    gen = DataGenerator(**kwargs)
    return gen.generate_data(cells, freqs, seed)
