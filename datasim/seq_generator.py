import random

import itertools as it
import numpy as np

from seq_data import SequencingData

class SequencingGenerator(object):
  """ Encapsulates parameters for generating sequencing data from a simulated experiment.
  Accepts the following parameters:
    Number of wells
      - num_wells [READ/WRITE]
    Distribution of cells/well
      - cells_per_well_distribution [READ ONLY]
      - cells_per_well_distribution_params [READ ONLY]
      - set_cells_per_well()
    Cell information
      - cells [READ/WRITE]
    Cell frequency distribution
      - cell_frequency_distribution [READ ONLY]
      - cell_frequency_distribution_params [READ ONLY]
      - set_cell_frequency_distribution()
    Noise stuff
      - chain_misplacement_prob [READ/WRITE]
  """ 

  def __init__(self, **kwargs):

    ## Set default parameter values
    self.num_wells = 96
    self.set_cells_per_well(distro_type = 'constant', cells_per_well = 35)
    self.num_cells = 1000
    self.cells = []
    self.set_cell_frequency_distribution(distro_type = 'power-law', alpha = -1)
    self.alpha_sharing_probs = None
    self.beta_sharing_probs = None
    self.chain_misplacement_prob = 0
    self.chain_deletion_prob = 0
    self.seed = np.random.randint(999999)

    ## If any parameters are specified in kwargs, modify them from the defaults
    self.set_options(**kwargs)

  ## Well parameters
  
  # Number of wells
  @property
  def num_wells(self):
    return self._num_wells
  @num_wells.setter
  def num_wells(self, val):
    self._num_wells = val

  # Distribution of cells/well
  @property
  def cells_per_well_distribution(self):
    return self._cpw_distro
  @property
  def cells_per_well_distribution_params(self):
    return self._cpw_params
  @property
  def cells_per_well(self):
    return self._cpw
  def set_cells_per_well(self, distro_type, **distro_params):
    self._cpw_distro = distro_type
    self._cpw_params = distro_params


  ## Cell parameters
  
  # Cell info - list of (alpha_id, beta_id) representing distinct cell types
  @property
  def cells(self):
    return self._cells
  @cells.setter
  def cells(self, val):
    self._cells = val[:]

  # Cell frequencies
  @property
  def cell_frequency_distribution(self):
    return self._cell_freq_distro
  @property
  def cell_frequency_distribution_params(self):
    return self._cell_freq_params
  def set_cell_frequency_distribution(self, distro_type, **distro_params):
    self._cell_freq_distro = distro_type
    self._cell_freq_params = distro_params


  ## Noise parameters
  # Probability that a chain will relocate to a random other well
  @property
  def chain_misplacement_prob(self):
    return self._cmp
  @chain_misplacement_prob.setter
  def chain_misplacement_prob(self, prob):
    self._cmp = prob
  # Probability that a chain will not be amplified properly (i.e. deleted)
  @property
  def chain_deletion_prob(self):
    return self._cdp
  @chain_deletion_prob.setter
  def chain_deletion_prob(self, prob):
    self._cdp = prob
  

  def set_options(self, **kwargs):
    if 'num_wells' in kwargs:
      self.num_wells = kwargs['num_wells']
    if 'cells_per_well_distribution' in kwargs and 'cells_per_well_distribution_params' in kwargs:
      self.set_cells_per_well(kwargs['cells_per_well_distribution'], **kwargs['cells_per_well_distribution_params'])
    if 'num_cells' in kwargs:
      self.num_cells = kwargs['num_cells']
    if 'cell_frequency_distribution' in kwargs and 'cell_frequency_distribution_params' in kwargs:
      self.set_cell_frequency_distribution(kwargs['cell_frequency_distribution'], **kwargs['cell_frequency_distribution_params'])
    if 'chain_misplacement_prob' in kwargs:
      self.chain_misplacement_prob = kwargs['chain_misplacement_prob']
    if 'chain_deletion_prob' in kwargs:
      self.chain_deletion_prob = kwargs['chain_deletion_prob']
    if 'alpha_sharing_probs' in kwargs:
      self.alpha_sharing_probs = kwargs['alpha_sharing_probs'] 
    if 'beta_sharing_probs' in kwargs:
      self.beta_sharing_probs = kwargs['beta_sharing_probs'] 
    if 'seed' in kwargs:
      self.seed = kwargs['seed'] 



  @staticmethod
  def generate_cells(num_cells, alpha_sharing_probs = None, beta_sharing_probs = None, alpha_dual_prob = 0.0, beta_dual_prob = 0.0, alpha_start_idx=0, beta_start_idx=0):

    # set random seed
    np.random.seed(self.seed)

    if alpha_sharing_probs is None:
      alpha_sharing_probs = [0.816,0.085,0.021,0.007,0.033,0.005,0.033]
    if beta_sharing_probs is None:
      beta_sharing_probs = [0.859,0.076,0.037,0.019,0.009]
    if isinstance(alpha_sharing_probs,float):
      alpha_sharing_probs = [alpha_sharing_probs,(1-alpha_sharing_probs)]
    if isinstance(beta_sharing_probs,float):
      beta_sharing_probs = [beta_sharing_probs,(1-beta_sharing_probs)]
         

    # Generate the degree for each alpha- and beta-chain from the given distribution
    adegs = np.random.choice(range(1,len(alpha_sharing_probs)+1), num_cells, replace=True, p=alpha_sharing_probs)
    bdegs = np.random.choice(range(1,len(beta_sharing_probs)+1), num_cells, replace=True, p=beta_sharing_probs)

    # Generate how many alpha- and beta-chains will be in each cell (i.e. how many dual chains)
    adual = [2 if v<=alpha_dual_prob else 1 for v in np.random.uniform(size=num_cells)]
    bdual = [2 if v<=beta_dual_prob else 1 for v in np.random.uniform(size=num_cells)]

    # Cut off at the desired number of cells
    alphas = [(i,) for i,n in enumerate(adegs) for _ in range(n)][:sum(adual)] # this truncation will alter the distro somewhat
    betas = [(i,) for i,n in enumerate(bdegs) for _ in range(n)][:sum(bdual)]

    # Randomly assign alpha- and beta-chains to each other
    np.random.shuffle(alphas)
    np.random.shuffle(betas)

    for i in range(num_cells):
      if adual[i]==2:  alphas[i:i+2] = [tuple(sorted(alphas[i]+alphas[i+1]))]
      if bdual[i]==2:  betas[i:i+2] = [tuple(sorted(betas[i]+betas[i+1]))]
    cells = list(set(zip(alphas, betas))) # Due to random duplicates, there may be slightly less than num_cells cells
    # (A slightly more complex method could be used to ensure exactly num_cells cells)
    # NOTE: working on this now (PVH)

    for i in xrange(len(cells),num_cells):
        #print 'Making adjustment for duplicate at index {}...'.format(i)
        cells.append(((i,),(i,)))

    return cells

  @staticmethod
  def generate_cells_lee(num_cells, max_alphas=None, max_betas=None):
    if max_alphas == None:  max_alphas = num_cells
    if max_betas == None:  max_betas = num_cells

    # Generate the degree for each alpha- and beta-chain from a given distribution
    sharing_probs=[0.8375, 0.0805, 0.029, 0.013, 0.021, 0.0025, 0.0165] # Averages from the Lee et al. paper
    adegs = np.random.choice(range(1,8), max_alphas, replace=True, p=sharing_probs)
    bdegs = np.random.choice(range(1,8), max_betas, replace=True, p=sharing_probs)

    # If you want to generate from a power law instead (not sure if this works as expected)
    #adegs = np.floor(np.random.pareto(2.1, size=max_alphas))+1
    #bdegs = np.floor(np.random.pareto(2.1, size=max_alphas))+1

    # Cut off at the desired number of cells
    alphas = sum([[i]*int(n) for i,n in enumerate(adegs)],[])[:num_cells] # this trunc. will skew the distro a bit
    betas = sum([[i]*int(n) for i,n in enumerate(bdegs)], [])[:num_cells]

    # Randomly assign alpha- and beta-chains to each other
    np.random.shuffle(alphas)
    np.random.shuffle(betas)
    cells = list(set(zip(alphas, betas))) # Due to chance dups, there may be slightly less than num_cells cells
    
    # (A slightly more complex method could be used to ensure exactly num_cells cells)
    return cells

  def _sample_cells_per_well(self):
    # Generate a list of number of cells for each well, based on the specified distribution
    # This is called to generate a new list each time generate_sequencing_data() is called.
    distro = self.cells_per_well_distribution
    params = self.cells_per_well_distribution_params
    if distro == 'constant':
      return [params['cells_per_well']] * self.num_wells
    elif distro == 'poisson':
      return list(np.random.poisson(params['lam'], self.num_wells))
    elif distro == 'explicit':
      return params['cells_per_well']
    else:
      assert False, "Unknown distribution of cells/well: {0}".format(distro)
  def _sample_cell_freqs(self):
    # Generate a list of cell frequencies based on the specified distribution
    # This is called each time generate_sequencing_data() is called.
    distro = self.cell_frequency_distribution
    params = self.cell_frequency_distribution_params
    if distro == 'constant':
      freqs = np.array([1]*self.num_cells)
    elif distro == 'power-law-old':
      freqs = np.random.pareto(-params['alpha'],self.num_cells) ## TODO: there's something screwy abt this distro, talk to holec abt it
    elif distro == 'power-law':
        freqs = 10.**(params['alpha']*np.log10(np.arange(1,self.num_cells+1)))
    elif distro == 'Lee':
      p_s = params.get('p_s', 0.5)
      n_s = params.get('n_s', 10)
      freq_min = p_s/(self.num_cells - n_s) # freq of clones in distro tail
      freq_n_s = 1.1*freq_min # lowest clone frequency within top p_s
      r = 2.*(p_s-freq_n_s*n_s)/((n_s-1)*n_s) # freq step size within top p_s
      freqs = [freq_n_s+r*i for i in range(n_s)] + [freq_min]*(self.num_cells-n_s)
      random.shuffle(freqs)
    elif distro == 'explicit':
      freqs = np.array(params['frequencies'])
    else:
      assert False, "Unknown distribution of cell frequencies: {0}".format(distro)

    freqs = freqs / np.sum(freqs) # Normalize freqs so it is a probability distro
    return list(freqs)

  def generate_data(self):
    # Generates sequencing data based on this SequencingGenerator object's parameter values.
    # Results are returned in a SequencingData object, which can be saved to a file with seq_data.save()

    # set random seed
    np.random.seed(self.seed)

    if self.cells == []:
        self.cells = SequencingGenerator.generate_cells(self.num_cells,
                alpha_sharing_probs = self.alpha_sharing_probs,
                beta_sharing_probs = self.beta_sharing_probs)

    cells_per_well = self._sample_cells_per_well()
    cell_freqs = self._sample_cell_freqs()

    misplaced_alphas = []
    misplaced_betas = []
    well_data = []
    cells_per_well_idx = []
    chain_deletions = []
    chain_misplacements = []
    for cpw in cells_per_well:
      # Pick cpw cells based on distro in cell_freqs
      # TODO: use trees to do this in O(log n) time?
      # (tbh: i actually don't know the efficiency of numpy's algorithm)
      cells_idx = np.random.choice(range(len(self.cells)), cpw, replace=True, p=cell_freqs)
      cells = [self.cells[idx] for idx in cells_idx]

      # Extract alpha and beta chains in the well
      alphas, betas = zip(*cells) if len(cells)>0 else ([],[])
      alphas, betas = [a for alist in alphas for a in alist], [b for blist in betas for b in blist]

      # Apply chain deletions and chain misplacements
      alphas_del_flag = [bool(v<self.chain_deletion_prob) for v in np.random.uniform(size=len(alphas))]
      alphas_misp_flag = [bool(v<self.chain_misplacement_prob) for v in np.random.uniform(size=len(alphas))]
      misplaced_alphas.extend([a for a,deleted,misplaced in zip(alphas, alphas_del_flag, alphas_misp_flag) if misplaced and not deleted])
      alphas = [a for a,deleted,misplaced in zip(alphas, alphas_del_flag, alphas_misp_flag) if not deleted and not misplaced]

      betas_del_flag = [bool(v<self.chain_deletion_prob) for v in np.random.uniform(size=len(betas))]
      betas_misp_flag = [bool(v<self.chain_misplacement_prob) for v in np.random.uniform(size=len(betas))]
      misplaced_betas.extend([b for b,deleted,misplaced in zip(betas, betas_del_flag, betas_misp_flag) if misplaced and not deleted])
      betas = [b for b,deleted,misplaced in zip(betas, betas_del_flag, betas_misp_flag) if not deleted and not misplaced] 

      # Remove duplicate chains and add to well_data
      well_data.append([sorted(set(alphas)), sorted(set(betas))])

      # Store actual cells in well for metadata
      cells_per_well_idx.append(list(cells_idx))

      # Store record of chain deletions and misplacements
      chain_deletions.append((alphas_del_flag, betas_del_flag))
      chain_misplacements.append((alphas_misp_flag, betas_misp_flag))

    # Put misplaced chains in random wells
    for a in misplaced_alphas:  well_data[np.random.randint(0, len(well_data))][0].append(a)
    for b in misplaced_betas:  well_data[np.random.randint(0, len(well_data))][1].append(b)
      
    metadata = {
      'num_wells': self.num_wells,
      'cells_per_well_distribution': self.cells_per_well_distribution,
      'cells_per_well_distribution_params': self.cells_per_well_distribution_params,
      'num_cells': self.num_cells,
      'cells': self.cells,
      'cell_frequency_distribution': self.cell_frequency_distribution,
      'cell_frequency_distribution_params': self.cell_frequency_distribution_params,
      'chain_deletion_prob':self.chain_deletion_prob,
      'generated_data': {
        'cells_per_well': cells_per_well,
        'cells_per_well_idx': cells_per_well_idx,
        'chain_deletions_per_well': chain_deletions,
        'chain_misplacements_per_well': chain_misplacements,
        'cell_frequencies': cell_freqs
      }
    }
    seq_data = SequencingData(well_data = well_data, metadata = metadata)
    return seq_data
