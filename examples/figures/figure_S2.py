
# standard libraries
import copy

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
from madhype import simulate_run 

# library modifications
plt.rcParams["font.family"] = "serif"

def main(*args,**kwargs):

    modifications = {
            'chain_deletion_prob':[.0,.1,.2,.3,.4,.5,.75],
            'chain_misplacement_prob':[.0,.1,.2,.3,.4,.5,.75]
            }

    repeats = 1000

    settings = default_settings()
    settings['cell_freq_max'] = 0.05
    settings['num_cells'] = 1000
    settings['cpw'] = (100,)
    settings['chain_deletion_prob'] = 0.0
    settings['chain_misplacement_prob'] = 0.0

    all_coverage = {}
    all_matches = {}

    solvers = ['madhype']
    solver_options = [{}]

    #
    for mod,values in modifications.items():

        all_results = []
        all_coverage[mod] = []
        all_matches[mod] = []

        for i,v in enumerate(values): 

            # iterate across system
            for r in xrange(repeats):

                specific_settings = copy.copy(settings)

                specific_settings[mod] = v
                specific_settings['seed'] = r

                _,results = simulate_run(solvers,solver_options,**specific_settings)

                all_results += results

            all_coverage[mod].append([results['frac_repertoire'] for results in all_results])
            all_matches[mod].append([results['positives'] for results in all_results])


    # plot/display settings
    fs = 18
    boxprops = dict(linewidth=3.0,zorder=1)
    meanlineprops = dict(linestyle='-',linewidth=2, color='black', zorder=0)
    plt.rcParams['xtick.labelsize'] = fs-4

    # figure specific properties
    fig,axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 12), sharey=False)
    plt.subplots_adjust(left=0.15,right=0.9,hspace=0.3,wspace=0.4)

    # set border for figure
    for axe in axes:
        for ax in axe:
            [i.set_linewidth(3) for i in ax.spines.itervalues()]

    labels = ['0%','10%','20%','30%','40%','50%','75%']
    
    bp = axes[0][0].boxplot(
            all_matches['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[0][0],'Chain Deletion Probability','Clonal Matches (#)',fs=fs)

    bp = axes[0][1].boxplot(
            all_coverage['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[0][1],'Chain Deletion Probability','Repertoire Coverage',fs=fs)

    bp = axes[1][0].boxplot(
            all_matches['chain_misplacement_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[1][0],'Chain Misplacement Probability','Clonal Matches (#)',fs=fs)

    bp = axes[1][1].boxplot(
            all_coverage['chain_misplacement_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[1][1],'Chain Misplacement Probability','Repertoire Coverage',fs=fs)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig2S.png', format='png', dpi=300)
    plt.close()


# --- Internal Methods --- #

def label_figure(ax,xlabel,ylabel,fs=18):

    ax.set_xlabel(xlabel,fontsize=fs)
    ax.set_ylabel(ylabel,fontsize=fs)

    if ylabel == 'Repertoire Coverage':
        ax.set_ylim((0,1))
        ax.set_yticks((0.,.5,1.))
        ax.set_yticklabels(('0%','50%','100%'),fontsize=fs)
    elif ylabel == 'Clonal Matches (#)':
        ax.set_ylim((0,1000))
        ax.set_yticks((0,500,1000))
        ax.set_yticklabels((0,500,1000),fontsize=fs)

def default_settings():
    return {
            # experimental design
            'num_cells':100,
            'num_wells':(96,),
            'analysis':('madhype',),
            # madhype parameters
            'threshold':0.5, # minimum ratio accepted by match_probability
            # alphabetr parameters
            'pair_threshold':0.0001,
            'iters':10,
            # simulation parameters
            'cell_freq_max':0.01,
            'cpw':(100,),
            'seed':1,
            # visual cues
            'silent':False,
            'visual':False,
            'visual_block':False,
            'compare':False
            }

if __name__ == "__main__":
    main()




