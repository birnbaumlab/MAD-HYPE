
"""
Tests out the new and improved variable solver
"""

import heapq

# intrapackage libraries
import simulation
from analysis.madhype import solve as madhype
from analysis.alphabetr import solve as alphabetr
from postprocessing import visualize_results, analyze_results, compare_results
from defaults import general_options as default_options

# create a dictionary called to map strings to function handles
solve = {
        'madhype':madhype,
        'alphabetr':alphabetr
        }

#------------------------------------------------------------------------------# 


def simulate_run(solvers, solver_options = None, **kwargs):
    # Generate datasets
    cells, cell_frequencies = simulation.generate_cells(**kwargs)
    data = simulation.generate_data(cells, cell_frequencies, **kwargs)

    return data, run(data, solvers, solver_options, **kwargs)

def run(data, solvers, solver_options = None, **kwargs):
    # Interpret solver_options to use defaults if None
    print('Starting run...')
    if solver_options is None:  solver_options = [{}]*len(solvers)

    # prepare for storage
    compiled_results = []

    # determine settings
    max_pairs = kwargs.get('max_pairs', default_options['max_pairs'])
    plot_comparison = kwargs.get('plot_comparison', default_options['plot_comparison'])

    # Solve using MAD-HYPE method
    for mode,mode_options in zip(solvers,solver_options):
        results = solve[mode](data,**mode_options)

        if len(results) > max_pairs:
            print 'Reducing number of returned pairs from {}->{} (declared limit)...'.format(
                    len(results), max_pairs)
            results = heapq.nlargest(max_pairs, results, key = lambda x: x[1]) 
        else:
#            print 'Starting results sorting by p-value...'
            results.sort(key=lambda x: -x[1])
#            print 'Finished!'

        # gather results
        compiled_results.append(analyze_results(results,data,**kwargs))

    # comparison if two methods are selected
    if plot_comparison != False:
        compare_results(compiled_results,**kwargs)

    # return compiled results
    return compiled_results


