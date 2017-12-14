
import numpy as np

num_cells = 100
alpha = 2.

freqs1 = 10.**(-alpha*np.log10(np.arange(1,num_cells+1)))
freqs2 = np.arange(1,num_cells+1) ** -alpha

print freqs1
print freqs2
print freqs1 - freqs2
