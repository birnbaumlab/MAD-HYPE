
"""
Makes the set of preliminary figures
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# homegrown libraries
from madhype.analysis.methods import match_probability
import madhype

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 16 


def main():

    """ Match probability as a function of w_a,w_b """

    fig = plt.figure(figsize=(12,10))

    grid = plt.GridSpec(2,2)

    ax1 = plt.subplot(grid[0, 0])
    ax2 = plt.subplot(grid[0, 1])
    ax3 = plt.subplot(grid[1, :])

    fig = plt.gcf()
    plt.subplots_adjust(wspace=0.45, hspace=0.3)

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
                 'alpha':2
                }

    plt.sca(ax1)

    for i in xrange(X.shape[0]):
        for j in xrange(Y.shape[1]):
            well_data['w_i'] = (X[i,j],)
            well_data['w_j'] = (Y[i,j],)
            val = 10**(-match_probability(well_data,prior)[0])
            Z[i,j] = 1./(1.+val)

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

    #ax.set_yticks(c_inds)
    #ax.set_yticklabels(c_labels)
    ax1.tick_params(width = 3,length=8,labelsize=18)
    plt.xticks([0,10,20])
    plt.yticks([0,10,20])

    plt.xlabel('$w_{i}$',fontsize=20)
    plt.ylabel('$w_{j}$',fontsize=20)
    #cbar.ax.set_yticklabels(['0%','100%'])  # vertically oriented colorbar

    # some default settings
    solvers = ['madhype']
    solver_options = [{}]
    settings = {
        'plot_auroc': True,
        'plot_auroc_options': {
            'ax':ax2,
            'fig':fig,
            },
        'plot_repertoire': True,
        'plot_repertoire_options': {
            'ax':ax3,
            'fig':fig,
            },
        'visual':                True,
        'silent':               False,
        }

    data,results = madhype.simulate_run(solvers, solver_options, **settings)

    plt.savefig('Figure3.png')

#    plt.show(block=False)
#    raw_input('Press enter to close...')
#    plt.close()
    

if __name__ == '__main__':
    main()
    
