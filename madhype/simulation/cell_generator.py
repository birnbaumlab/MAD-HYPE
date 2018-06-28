import numpy as np
import scipy.optimize

from ..defaults import general_options as default_options

class CellGenerator(object):
    """ Convenience class to store options for cell generation
    so that multiple cell populations can be generated easily with
    the same statistical properties. """

    def __init__(self, **kwargs):

        ## Store options
        self.options = default_options.copy()

        ## Update options from defaults
        self.options.update(kwargs)

        ## Interpret value of XXX_sharing_probs
        asp = self.options['alpha_sharing_probs']
        bsp = self.options['beta_sharing_probs']
        if asp is None:
            asp = [0.816,0.085,0.021,0.007,0.033,0.005,0.033]
        elif isinstance(asp, float):
            asp = [1-asp, asp]
        if bsp is None:
            bsp = [0.859,0.076,0.037,0.019,0.009]
        elif isinstance(bsp,float):
            bsp = [1-bsp, bsp]
        self.options['alpha_sharing_probs'] = asp
        self.options['beta_sharing_probs'] = bsp

    def generate_cells(self, seed = None):

        """ Generate cells from options"""

        # set random seed
        if seed is not None:
          np.random.seed(seed)

        # transfer options to local namespace
        num_cells = self.options['num_cells']
        cell_freq_distro = self.options['cell_freq_distro']
        alpha_sharing_probs = self.options['alpha_sharing_probs']
        beta_sharing_probs = self.options['beta_sharing_probs']
        alpha_dual_prob = self.options['alpha_dual_prob']
        beta_dual_prob = self.options['beta_dual_prob']

        # Generate the degree for each alpha- and beta-chain from the given distribution
        # PVH: using 2*number cells to ensure we get enough space for dual chains
        adegs = np.random.choice(
                range(1,len(alpha_sharing_probs)+1), 2*num_cells, replace=True, p=alpha_sharing_probs)
        bdegs = np.random.choice(
                range(1,len(beta_sharing_probs)+1), 2*num_cells, replace=True, p=beta_sharing_probs)

        # Generate how many alpha- and beta-chains will be in each cell (i.e. how many dual chains)
        adual = [2 if v<=alpha_dual_prob else 1 for v in np.random.uniform(size=num_cells)]
        bdual = [2 if v<=beta_dual_prob else 1 for v in np.random.uniform(size=num_cells)]

        # Cut off at the desired number of cells
        alphas = [(i,) for i,n in enumerate(adegs) for _ in range(n)][:sum(adual)] 
        betas = [(i,) for i,n in enumerate(bdegs) for _ in range(n)][:sum(bdual)]

        # Randomly assign alpha- and beta-chains to each other
        np.random.shuffle(alphas)
        np.random.shuffle(betas)

        for i in range(num_cells):
            # it's technically possible for alphas[i]==alphas[i+1] in which case
            # this won't work quite right --JB
            if adual[i]==2:  alphas[i:i+2] = [tuple(sorted(set(alphas[i]+alphas[i+1])))]
            if bdual[i]==2:  betas[i:i+2] = [tuple(sorted(set(betas[i]+betas[i+1])))]

        # Due to random duplicates, there may be slightly less than num_cells cells
        cells = list(set(zip(alphas, betas))) 

        if len(cells) < num_cells:
            print 'madhype: WARNING: Cell generation produced {} cells. Adding {} new cells (no chain-sharing or dual clones)'.format(len(cells), num_cells - len(cells))
        for i in xrange(len(cells),num_cells):
            cells.append(((i,),(i,))) # these chain ids will be greater than any already used

        # create frequencies associations
        if cell_freq_distro == 'constant':
            freqs = [1./num_cells for _ in xrange(num_cells)]
        elif cell_freq_distro == 'power-law':
            alpha = self.options['cell_freq_constant']
            freqs = _power_law_distribution(num_cells,alpha)

        return cells, freqs

def generate_cells(seed = None, **kwargs):
    gen = CellGenerator(**kwargs)
    return gen.generate_cells(seed)
        


#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _power_law_distribution(num_cells,alpha):
    """ Returns power law distribution using given parameters """ 

    freqs = (1 - np.arange(0, num_cells, 1.0) / num_cells) ** (1./(1-alpha))
    f_max = 1./freqs.sum()
    return freqs * f_max

