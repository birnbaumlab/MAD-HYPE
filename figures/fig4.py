
"""
This program is to make diagrams that are too frustrating to do by hand
"""

# standard libraries
import math

# nonstandard libraries
import matplotlib.pyplot as plt
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
                'madhype':[800,900],
                'alphabetr':[500,600]
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
                'madhype':[],
                'alphabetr':[]
                }
            }

    # settings

    fig,axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6), sharey=True)

    matches = [matches_by_cpw[100]['madhype'],
               matches_by_cpw[100]['alphabetr'],
               matches_by_cpw[1000]['madhype'],
               matches_by_cpw[1000]['alphabetr']]

    coverage = [coverage_by_cpw[100]['madhype'],
                coverage_by_cpw[100]['alphabetr'],
                coverage_by_cpw[1000]['madhype'],
                coverage_by_cpw[1000]['alphabetr']]

    labels = ['$N$ = 100','$N$ = 1000']*2

    # boxplot matches
    bplt1 = axes[0].boxplot(matches,labels=labels,color= 'black')
    axes[0].set_ylim((250,1000))
    axes[0].set_yticks((250,500,750,1000))

    # boxplot coverage 
    bplt2 = axes[1].boxplot(coverage,labels=labels)
    axes[1].set_ylim((0.,1.))
    axes[1].set_yticks((0.,.5,1.))
    axes[1].set_yticklabels(('0%','50%','100%'))

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()

if __name__ == "__main__":
    main()

