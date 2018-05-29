
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
from solver.methods import match_probability
from main import simulate_system 

# library modifications
plt.rcParams["font.family"] = "FreeSerif"


def main(*args,**kwargs):

    settings = {
                'subplots':[2]
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    """ Match probability as a function of w_a,w_b """

    if 1 in settings['subplots']:
        # TODO: Figure out better way to pass arguments
        """ Subplot 1 """

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

        fig, ax = plt.subplots(1,1)

        #cax = ax.imshow(data,interpolation='nearest')
        #ax.set_aspect(aspect='auto')

        pcm = ax.pcolor(X, Y, Z,
                       #norm=colors.LogNorm(vmin=1e-5, vmax=1),
                       norm=colors.Normalize(vmin=0, vmax=1),
                       cmap='PuBu_r')
        cbar = fig.colorbar(pcm, ax=ax)
        cbar.ax.tick_params(labelsize=16)
        cbar.set_ticks([0, .5, 1])
        cbar.set_ticklabels(['0%', '50%', '100%'])

        #ax.set_yticks(c_inds)
        #ax.set_yticklabels(c_labels)

        plt.xlabel('$w_{i}$',fontsize=20)
        plt.ylabel('$w_{j}$',fontsize=20)
        #cbar.ax.set_yticklabels(['0%','100%'])  # vertically oriented colorbar

    if 2 in settings['subplots']:
        #TODO: add functionality
        """  """
        simulate_system(settings)

    if 3 in settings['subplots']:
        #TODO: add functionality
        """ Subplot 3 """
        pass

    if 4 in settings['subplots']:
        #TODO: add functionality
        """ Subplot 4 """
        pass

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()
    

    
