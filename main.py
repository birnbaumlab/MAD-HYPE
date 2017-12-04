
"""
Tests out the new and improved variable solver
"""

# standard libraries

# nonstandard libraries

# homegrown libraries
from solver.methods import match_probability
from datasim.seq_generator import *

#------------------------------------------------------------------------------# 

def main():
    sg = SequencingGenerator()
    sg.set_cells_per_well(distro_type = 'explicit', cells_per_well = cpw_by_well)
    sg.set_options(**self.options)
    data = sg.generate_data()

def main(*args,**kwargs):

    print 'Starting optimization procedure...'


    '''
    Options: 
        > num_wells
        > cells_per_well_distibution
        > num_cells
        > cells
        > cell_frequency_distribution
        > chain_misplacement_prob
        > chain_deletion_prob
        > alpha_sharing_probs
        > beta_sharing_probs
    '''

    options = {
                'num_cells':100,
                'chain_deletion_prob':0.1,
                'alpha_sharing_probs':0.0,
                'beta_sharing_probs':0.0
              }

    # Convert well partitions to iterperable
    cpw_distro = [(p,w) for p,w in zip(parameters,self.options['well_partitions'])]
    cpw_by_well = [wt[0] for wt in cpw_distro for i in xrange(wt[1])]
    
    # Sequence Generator
    sg = SequencingGenerator()
    sg.set_options(**options)
    sg.set_cells_per_well(distro_type = 'explicit', cells_per_well = cpw_by_well)
    data = sg.generate_data()

    # generate results from solver, from data
    results = solve_madhype(data,pair_threshold=2.0)



#------------------------------------------------------------------------------# 

if __name__ == "__main__":
    main()


