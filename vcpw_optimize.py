
# Testing happens here

#---------------------------------------------------------------------------------#

# standard libraries

# nonstandard libraries
from matplotlib import pyplot as plt
import numpy as np

# homegrown libraries

from methods.madhype import solve as solve_madhype
from datasim.seq_generator import *

from variable import * # functions: test_settings
from analysis import * # functions: process_results

#---------------------------------------------------------------------------------#

def optimize_parameters(*args,**kwargs):

    print 'Starting optimization procedure...'

    options = {
                'step_size':0.1,
                'step_type':'relative', # relative,absolute
                'step_style':'all', # all,single
                'repeats':10, # tests at each step
                'max_iterations':10000,
                'max_parameter_repeat':10,
                'initial_guess':[10,16,48,161]
              }

    # update options with initial settings
    for arg in args: options.update(arg)
    options.update(kwargs)
    
    # Create model object
    fun = ModelFun(**options)
    current_parameters = options['initial_guess']
    parameter_log,coverage_log = [],[]

    # Start the engine
    for iter_num in xrange(options['max_iterations']):

        # Regenerate baseline information
        current_coverage = fun.test_parameters(current_parameters)

        # Generate all possible new steps
        step = options['step_size']
        up_step = [int(p*(1.+step)) if p != int(p*(1.+step)) else p + 1 for p in current_parameters]
        down_step = [int(p*(1.-step)) if p != int(p*(1.-step)) else p - 1 for p in current_parameters]

        if options['step_style'] == 'all':
            testing_parameters =  [[p if i != j else up_step[j] for j,p in enumerate(current_parameters)] 
                                                                for i in xrange(len(current_parameters))]
            testing_parameters += [[p if i != j else down_step[j] for j,p in enumerate(current_parameters)] 
                                                                  for i in xrange(len(current_parameters))]

        # Calculate coverage from available options
        coverage = [fun.test_parameters(parameters) for parameters in testing_parameters]

        # check if better option now exists  
        if max(coverage) > current_coverage:
            current_parameters = testing_parameters[coverage.index(max(coverage))]
            current_coverage = max(coverage)
            print 'New parameters: {} - {}%'.format(current_parameters,100*max(coverage))
        else:
            print 'New parameters: {} - {}%'.format(current_parameters,100*current_parameters)

        # Stores parameters in a log
        parameter_log.append(current_parameters)
        coverage_log.append(current_coverage)

        if parameter_log.count(current_parameters) > options['max_parameter_repeat']: 
            print 'Final discovered parameters:',current_parameters
            print 'Final coverage:',100*current_coverage,'%'

    # Stores logs of everything 
    results = {
                'parameter_log':parameter_log,
                'coverage_log':coverage_log,
                'function_log':fun.log
              }

    # Save results
    pickle.dump(results,open('optimization_results.p','w'))
    
    print 'Finished optimization procedure.'
        
#---------------------------------------------------------------------------------#


class ModelFun(object):

    def __init__(self,*args,**kwargs):

        """ Initialize  model function with default values, user input """

        self.options = {
                       'well_partitions':[24,24,24,24],
                       'num_cells':10000,
                       'chain_deletion_prob':0.1,
                       'alpha_sharing_probs':0.0,
                       'beta_sharing_probs':0.0,
                       'repeats':3
                       }

        self.log = [] # log for where algorithm has traversed

        # update options with initial settings
        for arg in args: self.options.update(arg)
        self.options.update(kwargs)


    def test_parameters(self,parameters):

        """ Tests a parameter set; returns repertoire coverage """
    
        coverage = []

        for i in xrange(self.options['repeats']):
        
            # assertion check
            assert len(parameters) == len(self.options['well_partitions']),'not enough parameters'

            # Create cpw partition list
            cpw_distro = [(p,w) for p,w in zip(parameters,self.options['well_partitions'])]
            cpw_by_well = [wt[0] for wt in cpw_distro for i in xrange(wt[1])]

            # Create data that fits user options
            sg = SequencingGenerator()
            sg.set_cells_per_well(distro_type = 'explicit', cells_per_well = cpw_by_well)
            sg.set_options(**self.options)
            data = sg.generate_data()

            # Start solver with options, dataset
            results = solve_madhype(data,pair_threshold=2.0)

            # Calculate coverage
            coverage.append(repertoire_coverage(results,data))
            
        # Gets avg, stores in log, and returns val
        avg_coverage = np.mean(coverage)
        self.log.append((p,avg_coverage))
        return avg_coverage

    def clear_log(self):
        """ Clears logging variable """ 
        self.log = []


#---------------------------------------------------------------------------------#


# namespace catch
if __name__ == '__main__':
    optimize_parameters()
