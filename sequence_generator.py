
"""

Author: Patrick V. Holec
Created: 4 December 2017

"""

"""

This is a new sequence generator that is compatible with variable cpw

"""

class SequencingGenerator(object):

    # I don't know if I wanna keep this stuff
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

    def __init__(self, *args, **kwargs):

        """ Initialize sequence generator object """

        ## Set default parameter values
        self.settings = {
                        'num_wells':96,
                        'cpw_distro':'constant',
                        'cells_per_well':35,
                        'num_cells':1000,
                        'cell_freq_distro':'power-law',
                        'cell_freq_constant':-1,
                        'chain_misplacement_prob':0,
                        'chain_deletion_prob':0,
                        'alpha_sharing_probs':None,
                        'beta_sharing_probs':None
                        }

        # update settings 
        self.update(*args,**kwargs)

    def update(*args,**kwargs):

        """ Update settings using *args and **kwargs """ 

        for arg in args: self.settings.update(arg)
        self.settings.update(kwargs)








