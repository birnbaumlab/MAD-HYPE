
"""

Author: Patrick V. Holec
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

#------------------------------------------------------------------------------# 

class DataGenerator(object):

    """
    Methods:
        __init__(self,*args,**kwargs)
            input: dicts, kwargs
            output: (none)
        update(self,*args,**kwargs)
            input: dicts, kwargs
            output: (none)
        generate_cells(self)
            input: (none)
            output: (none)
        generate_data(self)
            input: (none)
            output: data dict 
    """ 

    def __init__(self, *args, **kwargs):

        """ Initialize sequence generator object """

        ## Set default parameter values
        self.settings = {
                        'num_wells':96,
                        'cpw':35,
                        'num_cells':1000,
                        'cell_freq_distro':'power-law',
                        'cell_freq_constant':-1,
                        'chain_misplacement_prob':0,
                        'chain_deletion_prob':0.1,
                        'seed':42
                        }


        self.cells = None

        # update settings 
        self.update(*args,**kwargs)

    def update(self,*args,**kwargs):

        """ Update settings using *args and **kwargs """ 

        for arg in args: self.settings.update(arg)
        self.settings.update(kwargs)

    def generate_cells(self):

        """ Generate cells from user settings"""

        # set random seed
        np.random.seed(self.settings['seed'])

        # transfer settings to local namespace
        num_cells = self.settings['num_cells']
        cell_freq_distro = self.settings['cell_freq_distro']
        cell_freq_constant = self.settings['cell_freq_constant']

        # create local cell IDs
        self.cells = [((i,),(i,)) for i in xrange(num_cells)]
        #np.random.shuffle(self.cells) # FIXME

        # create frequencies associations
        if cell_freq_distro == 'constant':
            self.freqs = [1./num_cells for _ in xrange(num_cells)]
        elif cell_freq_distro == 'power-law':
            self.freqs = 10.**(-cell_freq_constant*np.log10(np.arange(1,num_cells+1)))
            self.freqs = self.freqs/sum(self.freqs)

    def generate_data(self):

        """ Generate data from user settings """

        # set random seed
        np.random.seed(self.settings['seed'])

        # transfer settings to local namespace
        num_cells = self.settings['num_cells']
        num_wells = self.settings['num_wells']
        cpw = self.settings['cpw']
        chain_deletion_prob = self.settings['chain_deletion_prob']

        # check if cells have already been generated
        if self.cells == None:
            self.generate_cells()

        # interpret user input for cpw distribution
        if isinstance(cpw,int) and isinstance(num_wells,int):
            self.cpw = [cpw for _ in xrange(num_wells)]
        if isinstance(cpw,(tuple,list)) and isinstance(num_wells,(tuple,list)):
            assert len(cpw) == len(num_wells),"cpw and num_wells not same length"
            self.cpw = [c for w,c in zip(num_wells,cpw) for _ in xrange(w)]

        # start generating well data
        self.well_data = {'A':[],'B':[]}
        for well_ind,cpw in enumerate(self.cpw):
            # select cell indices for each well
            indices = np.random.choice(num_cells,size=(cpw,),p=self.freqs)
            # seperate chains in to a,b
            a = list(itertools.chain(*[self.cells[i][0] for i in indices]))
            b = list(itertools.chain(*[self.cells[i][1] for i in indices]))
            # shuffle listed chains
            np.random.shuffle(a)
            np.random.shuffle(b)
            # truncation length, then remove duplicates
            a_ind,b_ind = binomial(len(a),chain_deletion_prob),binomial(len(b),chain_deletion_prob)
            a,b = set(a[a_ind:]),set(b[b_ind:])
            # store in wells
            self.well_data['A'].append(a)
            self.well_data['B'].append(b)

        # compile useful information
        data = {
                'well_data':self.well_data,
                'cells':dict([(c,f) for c,f in zip(self.cells,self.freqs)]),
                'settings':self.settings
               }

        return data # return results

#------------------------------------------------------------------------------# 
""" Factory Methods """
#------------------------------------------------------------------------------# 

#------------------------------------------------------------------------------# 
""" Namespace Catch """
#------------------------------------------------------------------------------# 
if __name__ == "__main__":
    sg = DataGenerator()
    sg.generate_cells()
    sg.generate_data()

