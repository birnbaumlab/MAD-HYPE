
"""
This program is to make diagrams that are too frustrating to do by hand
"""

import matplotlib.pyplot as plt
import math
import numpy as np

rows,columns = 8,12
cpw = 10
colors = ['b','g','r','c','m','y','b']

well_dev = 0.2
cell_size = 100
well_size = 60 
well_edge_size = 5
margin = (0.7,2)

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
        w_i = [(x,y) for x,y in zip([4,7,3,5,7,3,1,5,6,8,1,2],[2,2,5,9,12,11,5,1,4,9,9,6])]
        w_j = w_i[:6] + [(x,y) for x,y in zip([3,7,5,5,7,7],[7,7,4,4,9,2])]

        sets = [(w_i,'r',-0.1),(w_j,'b',0.1)]

        for w,c,shift in sets: 
            for x,y in w:
                dy = np.linspace(-0.25,0.25,101)
                dx = 0.075*np.sin(12*dy)
                plt.scatter(x + dy + 0, y + dx + shift,color=c,s=20) 




    plt.plot(x_mesh,y_mesh,'ko',color='black',markerfacecolor='white',markersize=well_size,zorder=0,markeredgewidth=well_edge_size)

    plt.xlim([0.5,rows+0.5])
    plt.ylim([0.5,columns+0.5])
    plt.axis('off')
    fig.patch.set_facecolor([0.9,0.9,0.9])
    plt.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

    #plt.savefig('myfig.png',facecolor=fig.get_facecolor(), edgecolor='none',bbox_inches='tight')
    plt.savefig('{}.png'.format(mode),facecolor=fig.get_facecolor(), edgecolor='k',edgewidth=5)

    #plt.show(block=False)
    #raw_input('Press enter to close...')
    #plt.close()

