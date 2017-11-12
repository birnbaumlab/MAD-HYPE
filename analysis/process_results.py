
"""
Attempt to test variable cells per well
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
from methods import *

def process_results(results,data,*args,**kwargs):

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
    assert len(results['cells']) == len(results['threshold']), "differing # cells/thresholds"
    
    #
    cells_with_scores = zip(results['cells'],results['threshold'])
    cells_with_scores = sorted(cells_with_scores,key=lambda x: -x[1])
    cells_without_scores = [i[0] for i in cells_with_scores]

    cells_record = set(data.metadata['cells'])
    cells_label = data.metadata['cells']
    cells_freqs = data.metadata['generated_data']['cell_frequencies']
    total_cells = len(cells_record) 

    # look at error rate (stratified by confidence) [FIGURE 1]
    x1,y1,fdr = [0],[0],0.01
    for c in cells_with_scores:
        if c[0] in cells_record:
            x1.append(x1[-1])
            y1.append(y1[-1]+1)
        else:
            x1.append(x1[-1]+1)
            y1.append(y1[-1])
        if x1[-1] > fdr*y1[-1]:
            break

    total_matches_at_fdr = x1[-1] + y1[-1]

    # look at error rate (stratified by confidence) [FIGURE 1]
    p2 = {'x':{'+':[],'-':[]},'y':{'+':[],'-':[]}}
    for i,c in enumerate(cells_record):
        print c
        if c in cells_without_scores[:total_matches_at_fdr]:
            p2['x']['+'].append(i)
            p2['y']['+'].append(cells_freqs[i])
        else:
            p2['x']['-'].append(i)
            p2['y']['-'].append(cells_freqs[i])

    print 'Total cells:',len(cells_record)
    print 'Matches at FDR = 0.01:',total_matches_at_fdr
    print 'Positives:',len(p2['x']['+'])
    print 'Negatives:',len(p2['x']['-'])

    # 
    plt.figure(1)
    plt.plot(x1,y1)

    plt.figure(2)
    plt.yscale('log')
    #plt.scatter(p2['x']['+'],p2['y']['+'],s=1,label='+')
    #plt.scatter(p2['x']['-'],p2['y']['-'],s=1,label='-')

    #plt.scatter(p2['x']['+'],p2['y']['+'],s=1,label='+')
    #plt.scatter(p2['x']['-'],p2['y']['-'],s=1,label='-')

    plt.legend()

    # show plots
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()


if __name__ == '__main__':
    test_settings()
