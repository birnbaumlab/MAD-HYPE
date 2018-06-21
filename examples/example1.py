## MAD-HYPE/examples/example1.py

# Import MAD-HYPE package
import madhype

# Set up run parameters
solvers = ['madhype','alphabetr']
solver_options = [{}, {}] # don't change default parameters

# Set up parameters that apply to all solvers/simulations
general_options = {
        'cpw':(10,),
        'num_wells':(96,),
        }

# Run MAD-HYPE with default parameters
data,results = madhype.simulate_run(solvers, solver_options, **general_options)

# Print out results
for solver, result in zip(solvers, results):

    print "{} Results:".format(solver)

    print "  Total # Cells:", result['total']
    print "  Chain pairs identified:", result['positives']
    print "  Chain pairs not identified:", result['negatives']

