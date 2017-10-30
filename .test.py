
from scipy.stats import binom
import numpy as np


a = [0.2,0.4,0.4]
b = [0.2,0.8]


print a
print b
print np.convolve(a,b,mode='full')

probs = binom.pmf(xrange(10),10,[0.1,0.2])
print probs

