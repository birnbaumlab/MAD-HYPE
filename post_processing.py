
"""
Attempt to test variable cells per well
"""

# standard libraries
from math import log10,ceil,floor

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
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
              'fdr':0.01,
              'silent':False
              }

    # update settings
    for arg in args: options.update(arg)
    options.update(kwargs)

    # assertion check
    
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
              'freqs':cells_freqs,
              'xy':[x1,y1]
              }

    if not options['silent']:
        # display characteristics of the data
        for k,v in results.items():
            print '{}:{}'.format(k,v)

    return results

#------------------------------------------------------------------------------# 

def visualize_results(results,data,*args,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    # default settings parameters
    settings = {
            'fdr':0.01,
            'fdr_plot':0.01,
            'pos_color':'black',
            'neg_color':'white',
            'legend':True,
            'silent':False
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    # heavy lifting on the data
    cresults = analyze_results(results,data,settings)

    # 

    # 
    ax = plt.figure().gca() # initialize figure

    linewidth = 5

    # get data attributes
    false_limit = cresults['xy'][0][-1]
    true_limit = cresults['xy'][1][-1]
    fdr_limit = min(false_limit,settings['fdr_plot']*true_limit)

    # plot FDR line
    plt.plot((0,fdr_limit),(0,fdr_limit/settings['fdr_plot']),
            linestyle='--',color='k',linewidth=linewidth)

    # plot main data auroc
    plt.plot(*cresults['xy'],linewidth=linewidth) 

    # add labels
    plt.xlabel('False positives (#)')
    plt.ylabel('True positives (#)')

    # change axes
    #plt.xticks(np.arange(0,false_limit+1,1))
    #plt.xlim(-0.1,1.1*false_limit)
    #plt.ylim(-0.1,1.1*true_limit)

    fig,ax = plt.subplots(figsize=(10,5))

    plt.yscale('log')

    for val,color,l in zip(
            (0,1),(settings['neg_color'],settings['pos_color']),('Correct','Incorrect')):
        xs = [i for i,p in enumerate(cresults['pattern']) if p == val]
        ys = [f for f,p in zip(cresults['freqs'],cresults['pattern']) if p == val]
        if len(xs) > 0: plt.bar(xs,ys,color=color,width=1,log=True,label=l,edgecolor='none')

    plt.plot(xrange(len(cresults['freqs'])),cresults['freqs'],color='k')

    if settings['legend'] == True:
        leg = plt.legend()
        leg.get_frame().set_edgecolor('k')


    plt.xlim((0,len(cresults['freqs'])))
    plt.ylim((min(cresults['freqs']),10**ceil(log10(max(cresults['freqs'])-1e-9))))
    ax.set_xticks([0,len(cresults['freqs'])/2,len(cresults['freqs'])])
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
