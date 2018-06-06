## MAD-HYPE/examples/example1.py

# Import MAD-HYPE package
import madhype

# Set up run parameters
solvers = ['madhype', 'alphabetr']
solver_options = [{}, {}] # don't change default parameters

# Run MAD-HYPE with default parameters
results = madhype.simulate_run(solvers, solver_options)

# Print out results
for solver, result in zip(solvers, results):
  print "{} Results:".format(solver)
  print "  Total # Cells:", result['total']
  print "  Chain pairs identified:", result['positives']
  print "  Chain pairs not identified:", result['negatives']
