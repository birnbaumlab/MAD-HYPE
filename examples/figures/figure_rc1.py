
"""
Makes the set of preliminary figures
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# homegrown libraries
from madhype.analysis.methods import match_probability
import madhype

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 16 


def main():

    settings = {
               }

    # update settings
    #for arg in args: settings.update(arg)
    #settings.update(kwargs)

    """ Match probability as a function of w_a,w_b """

    fig = plt.figure(figsize=(12,10))
    ax1 = plt.gca()

    # Spanning parameters
    w_i_range = (0,24)
    w_j_range = (0,24)
    prior = 1e-5
    
    # Create match built from parameters
    X, Y = np.mgrid[w_i_range[0]:w_i_range[1]+1,w_j_range[0]:w_j_range[1]+1]

    # Preallocate data matrix space
    Z = np.zeros((X.shape[0],Y.shape[1]))

    # Define default well distribution parameters
    well_data = {
                 'w_ij':(24,),
                 'w_o':(48,),
                 'w_tot':(96,),
                 'cpw':(10,),
                 'alpha':1
                }


    for i in xrange(X.shape[0]):
        for j in xrange(Y.shape[1]):
            well_data['w_i'] = (X[i,j],)
            well_data['w_j'] = (Y[i,j],)
            val = match_probability(well_data,0.01)[0]
            Z[i,j] = 1./(1.+prior*val)

    matplotlib.rcParams['xtick.labelsize'] = 16
    matplotlib.rcParams['ytick.labelsize'] = 16

    #cax = ax.imshow(data,interpolation='nearest')
    #ax.set_aspect(aspect='auto')

    pcm = ax1.pcolor(X, Y, Z,
                   norm=colors.Normalize(vmin=0, vmax=1),
                   cmap='PuBu_r')

    cbar = fig.colorbar(pcm, ax=ax1)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_ticks([0, .5, 1])
    cbar.set_ticklabels(['0%', '50%', '100%'])

    ax1.tick_params(width = 3,length=8,labelsize=18)
    plt.xticks([0,10,20])
    plt.yticks([0,10,20])

    plt.xlabel('$w_{i}$',fontsize=20)
    plt.ylabel('$w_{j}$',fontsize=20)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()
    

if __name__ == '__main__':
    main()
    
