
"""
Tests out the new and improved variable solver
"""

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


def simulate_run(solvers, solver_options, **kwargs):
    options = default_options.copy()

    # Update options
    options.update(kwargs)

    # Generate datasets
    cells, cell_frequencies = simulation.generate_cells(**options)
    data = simulation.generate_data(cells, cell_frequencies, **options)

    return data, run(data, solvers, solver_options, **options)

def run(data, solvers, solver_options, **kwargs):
    options = default_options.copy()

    # Update options
    options.update(kwargs)

    # prepare for storage
    compiled_results = []

    # Solve using MAD-HYPE method
    for mode,mode_options in zip(solvers,solver_options):
        results = solve[mode](data,**mode_options)

        print 'Starting results sorting by p-value...'
        results.sort(key=lambda x: -x[1])
        print 'Finished!'

        if len(results) > options['max_pairs']:
            print 'Reducing called pairs from {}->{} (declared limit)...'.format(
                    len(results),options['max_pairs'])
            results = results[:options['max_pairs']]

        # gather results
        if options['visual']:
            # visualize results if requested
            print 'Visualizing!'
            compiled_results.append(visualize_results(results,data,**options))
        else:
            # gather results
            compiled_results.append(analyze_results(results,data,**options))

    # comparison if two methods are selected
    if options['compare']:
        compare_results(compiled_results,data,**options)

    # return compiled results
    return compiled_results



