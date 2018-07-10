## MAD-HYPE/examples/example1.py
import matplotlib
matplotlib.use('Agg')

# Import MAD-HYPE package
import madhype

# Set up run parameters
solvers = ['madhype']
solver_options = [{}, {}] # don't change default parameters

# Set up parameters that apply to all solvers/simulations
general_options = {
        'cpw':(250,),
        'num_wells':(96,),
        'alpha_sharing_probs': 0.0,
        'num_cells': 100,
        'beta_sharing_probs': 0.0,
        'alpha_dual_prob': 0.5,
        'beta_dual_prob': 0.5,
        'visual':                     True,
        'plot_repertoire':            True,
        'save':                       True,
        #'display':                       True,
        }

# Run MAD-HYPE with default parameters
data,results = madhype.simulate_run(solvers, solver_options, **general_options)

# Print out results
for solver, result in zip(solvers, results):

    print "{} Results:".format(solver)

    print "  Total # Cells:", result['total']
    print "  Chain pairs identified:", result['positives']
    print "  Chain pairs not identified:", result['negatives']

