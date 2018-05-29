class CellGenerator(object):
    """ Convenience class to store settings for cell generation
    so that multiple cell populations can be generated easily with
    the same statistical properties. """

    default_settings = {
      'num_cells': 1000,
      'alpha_sharing_probs': None,
      'beta_sharing_probs': None,
      'alpha_dual_prob': 0.,
      'beta_dual_prob': 0.,
      'cell_freq_distro': 'power-law',
      'cell_freq_max': 0.01,
      'cell_freq_constant': -1
    }

    def __init__(self, **kwargs):

        ## Store settings
        self.settings = CellGenerator.default_settings.copy()

        ## Update settings from defaults
        self.settings.update(kwargs)

        ## Interpret value of XXX_sharing_probs
        asp, bsp = self.settings['alpha_sharing_probs'], self.settings['beta_sharing_probs']
        if not isinstance(asp, list):
            if asp is None:
                asp = [0.816,0.085,0.021,0.007,0.033,0.005,0.033]
            elif isinstance(asp, float):
                asp = [1-asp, asp]
        if not isinstance(bsp, list):
            if bsp is None:
                bsp = [0.859,0.076,0.037,0.019,0.009]
            if isinstance(bsp,float):
                bsp = [1-bsp, bsp]
        self.settings['alpha_sharing_probs'] = asp
        self.settings['beta_sharing_probs'] = bsp

    def generate_cells(self, seed = None):

        """ Generate cells from settings"""

        # set random seed
        if seed is not None:
          np.random.seed(self.settings['seed'])

        # transfer settings to local namespace
        num_cells = self.settings['num_cells']
        alpha_sharing_probs = self.settings['alpha_sharing_probs']
        beta_sharing_probs = self.settings['beta_sharing_probs']
        alpha_dual_prob = self.settings['alpha_dual_prob']
        beta_dual_prob = self.settings['beta_dual_prob']

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

        # create local cell IDs
        a_inds,b_inds = np.arange(num_cells),np.arange(num_cells)
        np.random.shuffle(a_inds)
        np.random.shuffle(b_inds)

        # Randomly assign alpha- and beta-chains to each other
        np.random.shuffle(alphas)
        np.random.shuffle(betas)

        print 'Length of alphas:',len(alphas)
        print 'Sum of adual:',sum(adual)
        print 'Sum of bdual:',sum(bdual)

        for i in range(num_cells):

          print 'Length of alphas:',len(alphas)

          # it's technically possible for alphas[i]==alphas[i+1] in which case
          # this won't work quite right --JB
          if adual[i]==2:  alphas[i:i+2] = [tuple(sorted(alphas[i]+alphas[i+1]))]
          if bdual[i]==2:  betas[i:i+2] = [tuple(sorted(betas[i]+betas[i+1]))]

        # Due to random duplicates, there may be slightly less than num_cells cells
        cells = list(set(zip(alphas, betas))) 

        for i in xrange(len(cells),num_cells):
            #print 'Making adjustment for duplicate at index {}...'.format(i)
            cells.append(((i,),(i,)))

        #self.cells = [((a_inds[i],),(b_inds[i],)) for i in xrange(num_cells)]
        #np.random.shuffle(self.cells) # FIXME

        # create frequencies associations
        if cell_freq_distro == 'constant':
            freqs = [1./num_cells for _ in xrange(num_cells)]
        elif cell_freq_distro == 'power-law':
            freq_max = self.settings['cell_freq_max']
            alpha = self.settings['cell_freq_constant']
            freqs =_power_law_distribution(num_cells,freq_max, alpha)

        return cells, freqs

def generate_cells(seed = None, **kwargs):
    gen = CellGenerator(**kwargs)
    return gen.generate_cells(seed)
        


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

#
