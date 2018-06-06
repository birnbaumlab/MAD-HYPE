# MAD-HYPE

This Python package provides a framework for analyzing multi-well T cell receptor (TCR) sequencing data to extract pairs of co-occurring &alpha; and &beta; chains (i.e. those chain pairs that appear together in a single T cell). Chain pairs are determined through a statistical analysis of the frequency with which each &alpha; and &beta; chain pair appear together in a well, adjusted by the overall frequency of each chain alone. A full description of the methodology used is described in (***CITE***) (in review). Also provided is a Python implementation of an alternative approach, ALPHABETR (***CITE***) for benchmarking purposes.

MAD-HYPE currently works with Python 2.7+ only. Python 3 implementation may occur in a future release.

## Installation

### Quick Install
In most cases, MAD-HYPE dependencies will be automatically installed with these instructions. Only attempt to install the dependencies manually if these instructions fail.
1. Download or clone the git repository locally.
2. In Terminal or Command Prompt, navigate to the base directory of the repository (MAD-HYPE).
3. Run:
   
   `$ python setup.py install`
   
   or
   
   `$ python setup.py install --user`
   
### Full Install
If dependencies are not installed correctly, you will need to manually install the following packages:
* NumPy >= 1.13.3
* SciPy >= 0.19.0
* matplotlib >= 2.0.0

These packages can be installed using the file `requirements.txt` included in the repository. In the base repository directory, run:

`$ pip install -r requirements.txt`

Then repeat Step 3 above.

## Getting Started

### Quick Start
Once installed, MAD-HYPE can be used by importing the package `madhype`. The following script simulates sequencing data for a randomly generated T cell population and analyzes it using both MAD-HYPE and ALPHABETR:

```
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
```

### Running MAD-HYPE on existing sequencing data
Existing sequencing data can be supplied to MAD-HYPE using `run()`, which takes a parameter `data` followed by the same parameters as `simulate_run()`. `data` is a `dict` that must, at a minimum, have the following key-value pairs:
  * `'well_data'`: a `dict` with key-value pairs:
    *  `'A'`: a list of `N` lists, where `N` is the number of wells and each inner list contains an ID for each &alpha; chain sequenced in that well
    *  `'B'`: a list of `N` lists, where each inner list contains an ID for each &beta; chain sequenced in that well
  *  `'options'`: a `dict` with options used to generate the sequencing data. Required key-value pairs are:
    *  `'num_wells'`: a tuple of the number of wells in each partition; for example, if all wells have the same number of cells (i.e. 1 partition), this should be a 1-tuple with the total number of wells.
    *  `'cpw'`: a tuple of the number of cells per well in each partition. Its length should match that of `well_data['options']['num_wells']`

Note that the &alpha; and &beta; chain IDs may be any hashable object (e.g. int, string). Overlap is allowed between &alpha; chain IDs and &beta; chain IDs.

Additionally, the order of num_wells and cpw are not community; that is, if num_wells = (36,60) and cpw = (100,1000), the first 36 elements of well_data['A'] and well_data['B'] will be analyzed as though there are 100 cells/well, and the subsequent 60 elements will be analyzed under 1000 cells/well.

The following example script loads sample sequencing data and runs MAD-HYPE on it, using the `run()` function.

```
## MAD-HYPE/examples/example2.py

# Import JSON parser
# This is used to import the sample sequencing data
import json

# Import MAD-HYPE package
import madhype

# Sample sequencing data is stored in JSON format
sequencing_data = json.load(open('sample_seq_data.txt'))

# Set up solver options
solvers = ['madhype'] # only run MAD-HYPE
solver_options = [{}] # use default parameters

# Run solvers
results = madhype.run(sequencing_data, solvers, solver_options)

# Print out results
for solver, result in zip(solvers, results):
  print "{} Results:".format(solver)
  print "  Total # Cells:", result['total']
  print "  Chain pairs identified:", result['positives']
  print "  Chain pairs not identified:", result['negatives']
```

### Using MAD-HYPE to generated simulated cell populations and sequencing data
Under the hood, `simulate_run()` generates a simulated TCR repertoire and simulated sequencing data. The functions used to generate sample cell populations and sequencing data are exposed in the `madhype.simulation` subpackage.

The following example script generates a sample population of 3000 cells. Simulated sequencing data is then generated from this cell population.

```
## MAD-HYPE/examples/example3.py

# Import MAD-HYPE package
import madhype

cells, cell_frequencies = madhype.simulation.generate_cells(num_cells = 3000)

sequencing_data = madhype.simulation.generate_data(cells, cell_frequencies, num_wells = 96, cpw = 100)
```

Further information about setting options for these functions is found in the section "Configuration Options".


## Configuration Options
The file `madhype/defaults.py` contains default settings for generating simulated cell populations and sequencing data and for analyzing sequencing data with MAD-HYPE or ALPHABETR. This file must be modified _prior to installation_ in order to change the default behavior of this package.

`madhype/defaults.py` contains three `dict`s: `general_options`, `madhype_options`, and `alphabetr_options`. Options can be changed from their defaults at run-time as follows:

* Values in `general_options` can be changed by supplying equivalent keyword arguments to `simulate_run()` and `run()`. For example, the following will perform a simulated run with 1000 cells.

```
madhype.simulate_run(..., num_cells = 1000)
```
   
* Values in `madhype_options` and `alphabetr_options` can be changed by supplying `dict`s with equivalent key-value pairs as the `solver_options` argument to `simulate_run()` and `run()`. For example, the following will run MAD-HYPE with an FDR filter of 20%.

```
solvers = ['madhype']
solver_options = [{
    'fdr': 0.2
}]
madhype.simulate_run(solvers, solver_options)
```
  
  
## Version
0.1.0

## Authors
Patrick V. Holec, Joseph Berleant
