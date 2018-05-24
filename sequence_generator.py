
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
from scipy.optimize import root

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
                        'num_wells':                 96,
                        'cpw':                       35,
                        'num_cells':               1000,
                        'cell_freq_distro': 'power-law',
                        'cell_freq_constant':        -1,
                        'cell_freq_max':           0.01,
                        'chain_misplacement_prob':    0,
                        'chain_deletion_prob':      0.1,
                        'alpha_dual_prob':           0.,
                        'beta_dual_prob':            0.,
                        'alpha_sharing_probs':     None,
                        'beta_sharing_probs':      None,
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
        cell_freq_max = self.settings['cell_freq_max']
        alpha_dual_prob = self.settings['alpha_dual_prob']
        beta_dual_prob = self.settings['beta_dual_prob']
        alpha_sharing_probs = self.settings['alpha_sharing_probs']
        beta_sharing_probs = self.settings['beta_sharing_probs']

        # interpret sharing chain probs
        if alpha_sharing_probs is None: 
            alpha_sharing_probs = [0.816,0.085,0.021,0.007,0.033,0.005,0.033]
        if beta_sharing_probs is None:
          beta_sharing_probs = [0.859,0.076,0.037,0.019,0.009]
        if isinstance(alpha_sharing_probs,float):
          alpha_sharing_probs = [(1-alpha_sharing_probs),(alpha_sharing_probs)]
        if isinstance(beta_sharing_probs,float):
          beta_sharing_probs = [(1-beta_sharing_probs),(beta_sharing_probs)]


        # Generate the degree for each alpha- and beta-chain from the given distribution
        adegs = np.random.choice(
                range(1,len(alpha_sharing_probs)+1), num_cells, replace=True, p=alpha_sharing_probs)
        bdegs = np.random.choice(
                range(1,len(beta_sharing_probs)+1), num_cells, replace=True, p=beta_sharing_probs)



        # Generate how many alpha- and beta-chains will be in each cell (i.e. how many dual chains)
        adual = [2 if v<=alpha_dual_prob else 1 for v in np.random.uniform(size=num_cells)]
        bdual = [2 if v<=beta_dual_prob else 1 for v in np.random.uniform(size=num_cells)]

        # Cut off at the desired number of cells
        alphas = [(i,) for i,n in enumerate(adegs) for _ in range(n)][:sum(adual)] 
        betas = [(i,) for i,n in enumerate(bdegs) for _ in range(n)][:sum(bdual)]

        # create local cell IDs
        a_inds,b_inds = np.arange(num_cells),np.arange(num_cells)
        np.random.shuffle(a_inds)
        np.random.shuffle(b_inds)


        # Randomly assign alpha- and beta-chains to each other
        np.random.shuffle(alphas)
        np.random.shuffle(betas)

        for i in range(num_cells):
          if adual[i]==2:  alphas[i:i+2] = [tuple(sorted(alphas[i]+alphas[i+1]))]
          if bdual[i]==2:  betas[i:i+2] = [tuple(sorted(betas[i]+betas[i+1]))]

        # Due to random duplicates, there may be slightly less than num_cells cells
        self.cells = list(set(zip(alphas, betas))) 

        for i in xrange(len(self.cells),num_cells):
            #print 'Making adjustment for duplicate at index {}...'.format(i)
            cells.append(((i,),(i,)))

        #self.cells = [((a_inds[i],),(b_inds[i],)) for i in xrange(num_cells)]
        #np.random.shuffle(self.cells) # FIXME

        # create frequencies associations
        if cell_freq_distro == 'constant':
            self.freqs = [1./num_cells for _ in xrange(num_cells)]
        elif cell_freq_distro == 'power-law':
            self.freqs =_power_law_distribution(num_cells,cell_freq_max,cell_freq_constant)

    def generate_data(self):

        """ Generate data from user settings """

        # set random seed
        np.random.seed(self.settings['seed'])

        # transfer settings to local namespace
        num_cells = self.settings['num_cells']
        num_wells = self.settings['num_wells']
        cpw = self.settings['cpw']
        chain_deletion_prob = self.settings['chain_deletion_prob']
        chain_misplacement_prob = self.settings['chain_misplacement_prob']

        # assertion check
        assert type(num_wells) == type(cpw), "num_wells and cpw differ in type"

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

        # aggregate misplaced chains
        misplaced_a = []
        misplaced_b = []

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

        self.well_data['A'] = [set(wd[binomial(len(wd),chain_deletion_prob):]) for wd in self.well_data['A']]
        self.well_data['B'] = [set(wd[binomial(len(wd),chain_deletion_prob):]) for wd in self.well_data['B']]

        # compile useful information
        data = {
                'well_data':self.well_data,
                'cells':dict([(c,f) for c,f in zip(self.cells,self.freqs)]),
                'settings':self.settings
               }

        return data # return results

#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _power_law_distribution(num_cells,max_freq,alpha):
    """ Returns power law distribution using given parameters """ 
    # Lower bound
    if max_freq <= 1./num_cells:
        print 'Max. freq too low! Returning uniform distribution...'
        return [1./num_cells for _ in xrange(num_cells)]
    
    # Upper bound
    if max_freq >= 1.:
        print 'Max. freq too high! Returning delta distribution...'
        return [1] + [0 for _ in xrange(num_cells-1)]
 
    # Find a shift
    shift = root(_get_max_freq_diff, 1.0, args=(num_cells,max_freq,alpha)).x
    
    # Find best
    return _get_freq_distribution(shift,num_cells,max_freq,alpha)

#------------------------------------------------------------------------------# 

def _get_max_freq_diff(shift,num_cells,max_freq,alpha):
    """ Function for finding diff. b/w max_freq and current distribution """
    freqs = _get_freq_distribution(shift,num_cells,max_freq,alpha)
    return max_freq - freqs[0]

#------------------------------------------------------------------------------# 

def _get_freq_distribution(shift,num_cells,max_freq,alpha):
    """ Generate a normalized power-law distribution """
    freqs = np.arange(shift,num_cells+shift) ** -alpha
    return freqs/sum(freqs)

#------------------------------------------------------------------------------# 
""" Namespace Catch """
#------------------------------------------------------------------------------# 
if __name__ == "__main__":
    sg = DataGenerator()
    sg.generate_cells()
    sg.generate_data()

