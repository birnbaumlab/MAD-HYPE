
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

def repertoire_coverage(results,data,*args,**kwargs):

    """ Inputs data,results; outputs percentage repertoire coverage """

    # default settings parameters
    settings = {
            'fdr':0.01
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    # put results in easy to access place
    cells_with_scores = zip(results['cells'],results['threshold'])
    cells_with_scores = sorted(cells_with_scores,key=lambda x: -x[1])
    cells_without_scores = [i[0] for i in cells_with_scores]

    cells_record = set(data.metadata['cells'])
    cells_label = data.metadata['cells']
    cells_freqs = data.metadata['generated_data']['cell_frequencies']
    total_cells = len(cells_record) 

    # look at error rate (stratified by confidence) [FIGURE 1]
    pos,neg = 0,0
    for c in cells_with_scores:
        if c[0] in cells_record: pos += 1
        else: neg += 1
        if neg > settings['fdr']*pos: break

    total_matches_at_fdr = pos + neg

    # look at error rate (stratified by confidence) [FIGURE 1]
    frac_repertoire = 0.

    for i,c in enumerate(cells_record):
        if c in cells_without_scores[:total_matches_at_fdr]:
            frac_repertoire += cells_freqs[i]

    # return fractional repertoire coverage
    return frac_repertoire



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

    cells_temp = sorted([(a,b) for a,b in data['cells'].items()])
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

    for i,c in enumerate(cells_record):
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
    plt.bar(x2,y2,color=colors,width=1,log=True)

    plt.legend()

    # show plots
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()


if __name__ == '__main__':
    test_settings()
