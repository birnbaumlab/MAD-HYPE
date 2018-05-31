
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

    run(data, solvers, solver_options, **options)

def run(data, solvers, solver_options, **kwargs):
    options = default_options.copy()

    # Update options
    options.update(kwargs)

    # prepare for storage
    compiled_results = []

    # Solve using MAD-HYPE method
    for mode in solvers:
        print mode
        results = solve[mode](data,**options) ## TODO: change options to **options

        results.sort(key=lambda x: -x[1])

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



