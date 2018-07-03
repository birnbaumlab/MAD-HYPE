
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
            'cell_freq_constant':[1.5,2.0,2.5,3.0],
            }

    labels = modifications['cell_freq_constant']
    #labels = ['0%','10%']

    repeats = 10

    settings = default_settings()

    settings['cell_freq_max'] = 0.05
    settings['num_cells'] = 1000
    settings['cpw'] = (50,1000)
    settings['num_wells'] = (48,48)
    settings['chain_deletion_prob'] = 0.1
    settings['chain_misplacement_prob'] = 0.0

    settings['alpha_dual_prob'] = 0.0
    settings['beta_dual_prob'] = 0.0

    settings['cell_freq_distro'] = 'power-law'

    all_coverage = {}
    all_matches = {}

    solvers = ['madhype']
    solver_options = [{}]

    for first_mod,values in modifications.items():

        for chain_sharing in [False,True]:

            if chain_sharing: mod = 'Chain sharing'
            else: mod = 'No chain sharing'

            all_coverage[mod] = []
            all_matches[mod] = []

            for i,v in enumerate(values): 

                all_results = []

                # iterate across system
                for r in xrange(repeats):

                    specific_settings = copy.copy(settings)

                    if first_mod == 'dual_prob':
                        specific_settings['alpha_dual_prob'] = v
                        specific_settings['beta_dual_prob'] = v
                    else:
                        specific_settings[first_mod] = v

                    if chain_sharing == True:
                        specific_settings['alpha_sharing_probs'] = None 
                        specific_settings['beta_sharing_probs'] = None 
                    else:
                        specific_settings['alpha_sharing_probs'] = 0.0 
                        specific_settings['beta_sharing_probs'] = 0.0 

                    specific_settings['seed'] = r

                    _,results = simulate_run(solvers,solver_options,**specific_settings)

                    all_results += results

                all_coverage[mod].append([results['frac_repertoire'] for results in all_results])
                all_matches[mod].append([float(results['positives']) for results in all_results])


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
    
    bp = axes[0][0].boxplot(
            all_matches['No chain sharing'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[0][0].set_title('No chain sharing',fontweight='bold',fontsize=fs)
    label_figure(axes[0][0],'Dual Clone Probability (%)','Clonal Matches (#)',fs=fs)

    bp = axes[1][0].boxplot(
            all_coverage['No chain sharing'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[1][0],'Dual Clone Probability (%)','Repertoire Coverage',fs=fs)

    bp = axes[0][1].boxplot(
            all_matches['Chain sharing'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[0][1].set_title('Chain sharing',fontweight='bold',fontsize=fs)
    label_figure(axes[0][1],'Dual Clone Probability (%)','Clonal Matches (#)',fs=fs)

    bp = axes[1][1].boxplot(
            all_coverage['Chain sharing'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[1][1],'Dual Clone Probability (%)','Repertoire Coverage',fs=fs)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig_S6.png', format='png', dpi=300)
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
    elif ylabel == 'Clonal Matches (%)':
        ax.set_ylim((0,1))
        ax.set_yticks((0.,.5,1.))
        ax.set_yticklabels(('0%','50%','100%'),fontsize=fs)

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




