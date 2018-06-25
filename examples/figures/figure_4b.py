
## MAD-HYPE/examples/figures/figure_4b.py

# nonstandard libraries
import matplotlib.pyplot as plt

# Import MAD-HYPE package
import madhype
from madhype.postprocessing.plots import plot_comparison

# Set up run parameters
solvers = ['madhype','alphabetr']
solver_options = [{}, {}] # don't change default parameters

# variants
cpws = [(100,),(30,)]

# Set up parameters that apply to all solvers/simulations
general_options = {
        'num_wells': (96,),
        }

fig,ax = plt.subplots(2,1,figsize = (10,10))

general_options['fig'] = fig

for i,cpw in enumerate(cpws):

    # set the number of cells per well
    general_options['cpw'] = cpw
    general_options['ax']  =  ax[i] 

    # Run MAD-HYPE with default parameters
    data,results = madhype.simulate_run(solvers, solver_options, **general_options)

    plot_comparison(results,**general_options) 

    # Print out results
    for solver, result in zip(solvers, results):

        print "{} Results:".format(solver)
        print "  Total # Cells:", result['total']
        print "  Chain pairs identified:", result['positives']
        print "  Chain pairs not identified:", result['negatives']

     






