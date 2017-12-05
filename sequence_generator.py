
"""

Author: Patrick V. Holec
Created: 4 December 2017

"""

"""

This is a new sequence generator that is compatible with variable cpw

"""

#------------------------------------------------------------------------------# 

class DataGenerator(object):

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
                        'chain_deletion_prob':0
                        }

        self.cells = None

        # update settings 
        self.update(*args,**kwargs)

    def update(*args,**kwargs):

        """ Update settings using *args and **kwargs """ 

        for arg in args: self.settings.update(arg)
        self.settings.update(kwargs)

    def generate_cells(self):

        """ Generate cells from user settings"""

        # transfer settings to local namespace
        num_cells = settings['num_cells']
        cell_freq_distro = settings['cell_freq_distro']
        cell_freq_constant = settings['cell_freq_constant']

        # create local cell IDs
        self.cells = [((i,),(i,)) for i in xrange(len(num_cells))]
        np.random.shuffle(self.cells)

        # create frequencies associations
        if cell_freq_distro == 'constant':
            self.freqs = [1./num_cells for _ in xrange(num_cells)]
        elif cell_freq_distro == 'power-law':
            self.freqs = 10.**(-cell_freq_constant*np.log10(np.arange(1,n+1)))
            self.freqs = self.freqs/sum(self.freqs)

    def generate_data(self):

        """ Generate data from user settings """

        # transfer settings to local namespace
        
                        'num_wells':96,
                        'cpw_distro':'constant',
                        'cells_per_well':35,
                        'num_cells':1000,
                        'cell_freq_distro':'power-law',

        # check if cells have already been generated
        if self.cells == None:
            self.generate_cells()



#------------------------------------------------------------------------------# 

if __name__ == "__main__":
    sg = DataGenerator()
    sh


