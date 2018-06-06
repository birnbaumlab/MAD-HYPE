## MAD-HYPE/examples/example2.py

# Import JSON parser
# This is used to import the sample sequencing data
import json

# Import MAD-HYPE package
import madhype

# Sample sequencing data is stored in JSON format
sequencing_data = json.load('sample_seq_data.txt')

# Set up solver options
solvers = ['madhype'] # only run MAD-HYPE
solver_optiosn = [{}] # use default parameters

# Run solvers
results = madhype.run(sequencing_data, solvers, solver_options)

# Print out results
for solver, result in zip(solvers, results):
  print "{} Results:".format(solver)
  print "  Total # Cells:", result['total']
  print "  Chain pairs identified:", result['positives']
  print "  Chain pairs not identified:", result['negatives']
