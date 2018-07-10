import numpy as np
import scipy.optimize

from ..defaults import general_options as default_options

# This flag should *almost* always be true, only used rarely when making frequency repertoire plots
shuffle_chains = True 

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
        if shuffle_chains:
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
            freq_max = self.options['cell_freq_max']
            alpha = self.options['cell_freq_constant']
            freqs = _power_law_distribution(num_cells,alpha, freq_max)

        return cells, freqs

def generate_cells(seed = None, **kwargs):
    gen = CellGenerator(**kwargs)
    return gen.generate_cells(seed)
        


#------------------------------------------------------------------------------# 
""" Internal Methods """
#------------------------------------------------------------------------------# 

def _power_law_distribution(num_cells, alpha, freq_max):
    """ Returns power law distribution using given parameters """ 

    # Choose freq_min so that sum(freqs)=1 (valid probability distribution)
    func = lambda x: (_power_law_constrained(num_cells, alpha, x, freq_max).sum() - 1)**2
    freq_min = scipy.optimize.minimize_scalar(func, bounds = (0.0, freq_max), method='bounded', options={'xatol':1e-100}).x
    freqs = _power_law_constrained(num_cells, alpha, freq_min, freq_max)

    if np.fabs(1 - freqs.sum()) / freqs[-1] > 1e-2:
        if freq_max > 1.0:
            print "Unable to generate power law probability distribution with freq_max={} > 0. Assuming delta distribution.".format(freq_max)
            freqs = [1.0] + [0.0]*(num_cells-1)
        elif freq_max < 1.0/num_cells:
            print "Unable to generate power law probability distribution with freq_max*num_cells < 1.0. Using a uniform distribution."
            freqs = [1.0/num_cells] * num_cells
        else:
            print "Unable to generate a power law probability distribution with alpha={} and freq_max={}. Try using default values alpha={}, freq_max={}.".format(alpha, freq_max, default_options['cell_freq_constant'], default_options['cell_freq_max'])
            freqs = [float('nan')] * num_cells

    return freqs

def _power_law_constrained(num_points, alpha, vmin, vmax):
    """ Returns values following a power law constrained between vmin and vmax with exponent alpha """
    vmax_adj = vmax ** (1-alpha)
    vmin_adj = vmin ** (1-alpha)
    data = (vmax_adj - (vmax_adj - vmin_adj) * np.arange(0, num_points, 1.0) / num_points) ** (1. / (1-alpha))
    return data
