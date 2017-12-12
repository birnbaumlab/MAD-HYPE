
"""
Attempt to test variable cells per well
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
#from methods import *


'''
MAIN FUNCTIONS
'''

def visualize_results(results,data,*args,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    # default settings parameters
    settings = {
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    # assertion check
    #assert len(results['cells']) == len(results['threshold']), "differing # cells/thresholds"
    
    #
    cells_with_scores = sorted(results,key=lambda x: -x[1])
    cells_without_scores = [i[0] for i in cells_with_scores]

    # these are actuall cells
    cells_temp = sorted([(a,b) for a,b in data['cells'].items()],key=lambda x: -x[1])
    cells_record = set([c[0] for c in cells_temp])
    cells_label = [c[0] for c in cells_temp]
    cells_freqs = [c[1] for c in cells_temp]
    total_cells = len(cells_temp) 

    # look at error rate (stratified by confidence) [FIGURE 1]
    x1,y1,fdr = [0],[0],0.01
    for c in cells_with_scores:
        try:
            data['cells'][c[0]]
            x1.append(x1[-1])
            y1.append(y1[-1]+1)
        except KeyError:
            x1.append(x1[-1]+1)
            y1.append(y1[-1])
        if x1[-1] > fdr*y1[-1]:
            break

    total_matches_at_fdr = x1[-1] + y1[-1]

    # look at error rate (stratified by confidence) [FIGURE 1]
    freqs_colors,frac_repertoire = [],0.
    x2 = []

    for i,c in enumerate(cells_label):
        if c in cells_without_scores[:total_matches_at_fdr]:
            freqs_colors.append((cells_freqs[i],'green'))
            frac_repertoire += cells_freqs[i]
        else:
            freqs_colors.append((cells_freqs[i],'red'))
        x2.append(i)

    freqs_colors = sorted(freqs_colors,key=lambda x: -x[0])
    y2 = [i[0] for i in freqs_colors]
    colors = [i[1] for i in freqs_colors]

    print 'Total cells:',len(cells_record)
    print 'Matches at FDR = 0.01:',total_matches_at_fdr
    print 'Positives:',colors.count('green')
    print 'Negatives:',colors.count('red')
    print 'Fraction of repertoire:',frac_repertoire

    # 
    plt.figure(1)
    plt.plot(x1,y1)

    plt.figure(2)

    for l,color in zip(('Correct','Incorrect'),('green','red')):
        xs = [i for i,c in zip(x2,colors) if c == color]
        ys = [i for i,c in zip(y2,colors) if c == color]
        print 'xs:',xs
        print 'ys:',ys
        print colors
        if len(xs) > 0: plt.bar(xs,ys,color=color,width=1,log=True,label=l)

    plt.annotate('{}/{} identified'.format(colors.count('green'),len(cells_record)),
            xy=(0.7,0.6),xycoords='axes fraction')
    plt.annotate('{}% repertoire'.format(round(100*frac_repertoire,2)),
            xy=(0.7,0.55),xycoords='axes fraction')
    plt.legend()

    # show plots
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()


if __name__ == '__main__':
    test_settings()
