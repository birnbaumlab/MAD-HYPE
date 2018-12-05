## MAD-HYPE/examples/example3.py

# Import MAD-HYPE package
import madhype

cells, cell_frequencies = madhype.simulation.generate_cells(num_cells = 3000)

sequencing_data = madhype.simulation.generate_data(cells, cell_frequencies, num_wells = 96, cpw = 100)
