
"""
Makes the set of preliminary figures
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
from solver.methods import match_probability

def main(*args,**kwargs):
    """ Match probability as a function of w_a,w_b """

    w_i_range = np.linspace(0,37)
    w_j_range = np.linspace(0,37)
    data = np.zeros((len(w_i_range),len(w_j_range)))

    well_data = {'w_i':(0,),
                 'w_j':(0,),
                 'w_ij':(24,),
                 'w_o':(48,),
                 'w_tot':(96,),
                 'cpw':(10,),
                 'alpha':1}

    for i,w_i in enumerate(w_i_range):
        for j,w_j in enumerate(w_j_range):
            well_data['w_i'] = (w_i,)
            well_data['w_j'] = (w_j,)
            #print match_probability(well_data,0.01)[0]
            data[i,j] = np.log10(match_probability(well_data,0.01)[0])


    fig, ax = plt.subplots()

    cax = ax.imshow(data,interpolation='nearest')
    ax.set_aspect(aspect='auto')

    #ax.set_yticks(c_inds)
    #ax.set_yticklabels(c_labels)

    plt.xlabel('$w_{i}$')
    plt.ylabel('$w_{j}$')
    cbar = fig.colorbar(cax)
    #cbar.ax.set_yticklabels(['0%','100%'])  # vertically oriented colorbar

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()
    

    
