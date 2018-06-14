
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

    return run(data, solvers, solver_options, **options)

def run(data, solvers, solver_options, **kwargs):
    general_options = default_options.copy()

    # Update options
    general_options.update(kwargs)

    # prepare for storage
    compiled_results = []

    # Solve using MAD-HYPE method
    for mode,mode_options in zip(solvers,solver_options):

        options = general_options.copy()
        options.update(mode_options)

        results = solve[mode](data,**options)

        if len(results) > options['max_pairs']:
            print 'Reducing called pairs from {}->{} (declared limit)...'.format(
                    len(results),options['max_pairs'])
            results = results[:options['max_pairs']]

        print 'Starting results sorting by p-value...'
        results.sort(key=lambda x: -x[1])
        print 'Finished!'

        # if there are no references for correct sequences
        if 'cells' in data:
            compiled_results.append(results)
            continue

        # gather results
        if options['visual']:
            # visualize results if requested
            print 'Visualizing!'
            compiled_results.append(visualize_results(results,data,**options))
        else:
            # gather results
            compiled_results.append(analyze_results(results,data,**options))

    # comparison if two methods are selected
    if general_options['compare']:
        compare_results(compiled_results,data,**general_options)

    # return compiled results
    return compiled_results



