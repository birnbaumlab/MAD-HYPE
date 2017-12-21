
"""
Attempt to test variable cells per well
"""

# standard libraries
from math import log10,ceil,floor

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
#from methods import *
plt.rcParams["font.family"] = "serif"

'''
MAIN FUNCTIONS
'''

def analyze_results(results,data,*args,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    # default settings parameters
    options = {
              'fdr':0.01
              }

    # update settings
    for arg in args: options.update(arg)
    options.update(kwargs)

    # assertion check
    #assert len(results['cells']) == len(results['threshold']), "differing # cells/thresholds"
    
    # these are guessed cells
    cells_with_scores = sorted(results,key=lambda x: -x[1])
    cells_without_scores = [i[0] for i in cells_with_scores]

    # these are actual cells
    cells_temp = sorted([(a,b) for a,b in data['cells'].items()],key=lambda x: -x[1])
    cells_record = set([c[0] for c in cells_temp])
    cells_label = [c[0] for c in cells_temp]
    cells_freqs = [c[1] for c in cells_temp]
    total_cells = len(cells_temp)

    # look at error rate (stratified by confidence)
    x1,y1,fdr = [0],[0],0.01
    for c in cells_with_scores:
        try:
            data['cells'][c[0]]
            x1.append(x1[-1])
            y1.append(y1[-1]+1)
        except KeyError:
            x1.append(x1[-1]+1)
            y1.append(y1[-1])
        if x1[-1] > options['fdr']*y1[-1]:
            break

    total_matches_at_fdr = x1[-1] + y1[-1]

    freqs = [cells_freqs[i] for i,c in enumerate(cells_label) \
             if c in cells_without_scores[:total_matches_at_fdr]]
    frac_repertoire = sum(freqs)
    pattern = [1 if c in cells_without_scores[:total_matches_at_fdr] else 0 \
               for i,c in enumerate(cells_label)]

    results = {
              'pattern':pattern,
              'frac_repertoire':frac_repertoire,
              'positives':sum(pattern),
              'negatives':len(pattern) - sum(pattern),
              'total':len(pattern),
              'freqs':cells_freqs
              }

    return results

#------------------------------------------------------------------------------# 

def visualize_results(results,data,*args,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    # default settings parameters
    settings = {
            'pos_color':'black',
            'neg_color':'white',
            'legend':True
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
            freqs_colors.append((cells_freqs[i],settings['pos_color']))
            frac_repertoire += cells_freqs[i]
        else:
            freqs_colors.append((cells_freqs[i],settings['neg_color']))
        x2.append(i)

    freqs_colors = sorted(freqs_colors,key=lambda x: -x[0])
    y2 = [i[0] for i in freqs_colors]
    colors = [i[1] for i in freqs_colors]

    print 'Total cells:',len(cells_record)
    print 'Matches at FDR = 0.01:',total_matches_at_fdr
    print 'Positives:',colors.count(settings['pos_color'])
    print 'Negatives:',colors.count(settings['neg_color'])
    print 'Fraction of repertoire:',frac_repertoire

    # 
    plt.plot(x1,y1)

    fig,ax = plt.subplots(figsize=(10,5))

    plt.yscale('log')

    for l,color in zip(('Correct','Incorrect'),(settings['pos_color'],settings['neg_color'])):
        xs = [i for i,c in zip(x2,colors) if c == color]
        ys = [i for i,c in zip(y2,colors) if c == color]
        if len(xs) > 0: plt.bar(xs,ys,color=color,width=1,log=True,label=l,edgecolor='none')

    plt.plot(xrange(len(cells_freqs)),cells_freqs,color='k')

    #plt.annotate('{}/{} identified'.format(colors.count('green'),len(cells_record)),
    #        xy=(0.7,0.6),xycoords='axes fraction')
    #plt.annotate('{}% repertoire'.format(round(100*frac_repertoire,2)),
    #        xy=(0.7,0.55),xycoords='axes fraction')
    if settings['legend'] == True:
        leg = plt.legend()
        leg.get_frame().set_edgecolor('k')


    plt.xlim((0,len(x2)))
    plt.ylim((min(y2),10**ceil(log10(max(y2)-1e-9))))
    ax.set_xticks([0,len(x2)/2,len(x2)])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Clone #',fontweight='bold')
    plt.ylabel('',fontweight='bold',size=12)
    ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(labelsize=20)

    # show plots
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()


if __name__ == '__main__':
    test_settings()
