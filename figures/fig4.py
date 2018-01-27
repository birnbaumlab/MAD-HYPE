
"""
This program is to make diagrams that are too frustrating to do by hand
"""

# standard libraries
import math

# nonstandard libraries
import matplotlib.pyplot as plt
from pylab import setp
import numpy as np

# library modifications
plt.rcParams["font.family"] = "serif"

def main(*args,**kwargs):

    """ Makes diagram figures """

    matches_by_cpw = {
            100:{
                'madhype':[697,664,699,718,740,723,727,693,734,724],
                'alphabetr':[436,430,525,448,450,542,434,496,444,563]
                },
            1000:{
                'madhype':[945,937,941,948,946,941,940,942,946,941],
                'alphabetr':[903,47,914,932,920,89,920,922,899,918]
                }
            }

    coverage_by_cpw = {
            100:{
                'madhype':
                [0.877930962132,0.864676101089,0.881788598911,0.891740493664,0.896221893956,
                 0.889533656519,0.894632759173,0.878798270374,0.899417127207,0.891927331392],
                'alphabetr':
                [0.748020305916,0.752852033505,0.802814157319,0.760043654059,0.760682382822,
                 0.809304374966,0.752895331150, 0.787971286782,0.757627278537,0.822973830961]
                },
            1000:{
                'madhype':
                [0.700650595765,0.67367769407,0.690086260119,0.70870390341,0.703248656572,
                 0.690300620824,0.683764563204,0.691218875063,0.707222467948,0.690284851379],
                'alphabetr':
                [0.739069450512,0.0420614543952,0.748443960839,0.798745740671,0.761087857233,
                 0.0747117077914,0.768919011223,0.76105297092,0.755368621068,0.764853043629],
                }
            }

    # settings
    boxprops = dict(linewidth=3.0,zorder=1)
    meanlineprops = dict(linestyle='-',linewidth=2, color='black', zorder=0)
    fs = 18

    # figure specific properties
    fig,axes = plt.subplots(nrows=2, ncols=1, figsize=(5, 12), sharey=False)
    plt.subplots_adjust(left=0.3,right=0.9,hspace=1.0,wspace=1.0)

    # set border for figure
    for ax in axes:
        [i.set_linewidth(3) for i in ax.spines.itervalues()]

    matches = [matches_by_cpw[100]['madhype'],
               matches_by_cpw[100]['alphabetr'],
               matches_by_cpw[1000]['madhype'],
               matches_by_cpw[1000]['alphabetr']]

    coverage = [coverage_by_cpw[100]['madhype'],
                coverage_by_cpw[100]['alphabetr'],
                coverage_by_cpw[1000]['madhype'],
                coverage_by_cpw[1000]['alphabetr']]

    labels = ['$N$ = 100','$N$ = 100','$N$ = 1000','$N$ = 1000']

    # boxplot matches
    bp = axes[0].boxplot(matches, labels=labels, boxprops=boxprops, meanprops=meanlineprops, widths=0.6, meanline=True, showmeans=True)
    axes[0].plot((2.5,2.5),(0,1000),linestyle='--',color='k')
    setBoxColors(bp)
 
    axes[0].set_xticks((1.5,3.5))
    axes[0].set_xticklabels(('N = 100','N = 1000'),fontsize=fs)
    axes[0].set_xlabel('Cells/well',fontsize=fs)

    axes[0].set_ylim((0,1000))
    axes[0].set_yticks((0,250,500,750,1000))
    axes[0].set_yticklabels((0,250,500,750,1000),fontsize=fs)
    axes[0].set_ylabel('Clonal Matches (#)',fontsize=fs)

    # boxplot coverage 
    bp = axes[1].boxplot(coverage, labels=labels, boxprops=boxprops, meanprops=meanlineprops, widths=0.6, meanline=True, showmeans=True)
    axes[1].plot((2.5,2.5),(0,1),linestyle='--',color='k')
    setBoxColors(bp)

    axes[1].set_xticks((1.5,3.5))
    axes[1].set_xticklabels(('N = 100','N = 1000'),fontsize=fs)
    axes[1].set_xlabel('Cells/well',fontsize=fs)

    axes[1].set_ylim((0.,1.))
    axes[1].set_yticks((0.,.5,1.))
    axes[1].set_yticklabels(('0%','50%','100%'),fontsize=fs)
    axes[1].set_ylabel('Repertoire Coverage',fontsize=fs)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig4C.png', format='png', dpi=200)
    plt.close()

# --- Internal Methods --- #

def setBoxColors(bp):
    """ function for setting the colors of the box plots pairs """
    for i in [0,1]:
        setp(bp['boxes'][2*i+0], color='blue')
        setp(bp['caps'][4*i+0], color='blue')
        setp(bp['caps'][4*i+1], color='blue')
        setp(bp['whiskers'][4*i+0], color='blue')
        setp(bp['whiskers'][4*i+1], color='blue')
        setp(bp['fliers'][2*i+0], color='blue')
        ##setp(bp['fliers'][1], color='blue')
        setp(bp['medians'][2*i+0], color='blue')

        setp(bp['boxes'][2*i+1], color='red')
        setp(bp['caps'][4*i+2], color='red')
        setp(bp['caps'][4*i+3], color='red')
        setp(bp['whiskers'][4*i+2], color='red')
        setp(bp['whiskers'][4*i+3], color='red')
        #setp(bp['fliers'][i+1], color='red')
        #setp(bp['fliers'][3], color='red')
        setp(bp['medians'][2*i+1], color='red')

if __name__ == "__main__":
    main()

