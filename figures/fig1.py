
"""
This program is to make diagrams that are too frustrating to do by hand
"""

import matplotlib.pyplot as plt
import math
import numpy as np

def main(*args,**kwargs):
    """ Makes diagram figures """
    rows,columns = 8,12
    cpw = 10
    colors = ['b','g','r','c','m','y','b']

    well_dev = 0.2
    cell_size = 100
    well_size = 60 
    well_edge_size = 5
    margin = (0.7,2)

    w_i = 5
    w_j = 7
    w_ij = 23
    w_o = 61

    x = np.linspace(1,rows,rows)
    y = np.linspace(1,columns,columns)

    x_mesh,y_mesh = np.meshgrid(x,y)

    for mode in ['cells','dna']:#['cells','dna']:
        np.random.seed(42)

        fig = plt.figure(figsize=(rows+margin[0], 0.9*columns+margin[1]), dpi=200, facecolor='w', edgecolor='k')
        ax = plt.gca()

         
        if mode == 'cells':
            for i,j in zip(x_mesh.flatten(),y_mesh.flatten()):
                    well_colors = [colors[np.random.randint(0,len(colors))] for _ in xrange(cpw)]
                    r, theta = [well_dev*np.sqrt(np.random.rand(cpw,1)), 2*math.pi*np.random.rand(cpw,1)]
                    x_mod,y_mod = r*np.cos(theta),r*np.sin(theta)
                    x_well = i + x_mod
                    y_well = j + y_mod

                    plt.scatter(x_well,y_well,color=well_colors,s=cell_size,zorder=10)
        elif mode == 'dna':
            combos = [(i,j) for i in xrange(1,rows+1) for j in xrange(1,columns+1)]
            np.random.shuffle(combos)
            print combos

            w_i_indices = combos[:w_i+w_ij]
            w_j_indices = combos[w_i:w_i+w_j+w_ij]

            sets = [(w_i_indices,'r',-0.1),(w_j_indices,'b',0.1)]

            for w,c,shift in sets: 
                for x,y in w:
                    dy = np.linspace(-0.25,0.25,101)
                    dx = 0.075*np.sin(12*dy)
                    plt.scatter(x + dy + 0, y + dx + shift,color=c,s=20) 




        plt.plot(x_mesh,y_mesh,'ko',color='black',markerfacecolor='white',markersize=well_size,zorder=0,markeredgewidth=well_edge_size)

        plt.xlim([0.5,rows+0.5])
        plt.ylim([0.5,columns+0.5])
        plt.axis('off')
        #fig.patch.set_facecolor([1.0,1.0,1.0])
        fig.patch.set_facecolor([0.9,0.9,0.9])
        plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

        #plt.savefig('myfig.png',facecolor=fig.get_facecolor(), edgecolor='none',bbox_inches='tight')
        plt.savefig('{}.png'.format(mode),facecolor=fig.get_facecolor(), edgecolor='k',edgewidth=5)

        #plt.show(block=False)
        #raw_input('Press enter to close...')
        #plt.close()


if __name__ == '__main__':
    main()


