
# make sure a backend is ok on AWS
import matplotlib
matplotlib.use('Agg')

# standard libraries
import os
import pickle
import copy

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
from madhype import simulate_run 

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 14 
plt.rcParams['ytick.labelsize'] = 14


def _get_fdr_for_match_limit(data,results,match_limit):

    positives,negatives = 0,0

    for i,r in enumerate(results[0]['raw_results'][:match_limit]):

        if r[0] in data['cells']:
            positives += 1
        else:
            negatives += 1

    negatives += match_limit - (positives + negatives)

    return min(0.25,float(negatives)/(positives + negatives))



def main(*args,**kwargs):

    mod_range = [.1,.15,.2,.25,.3,.35,.4]
    labels = ['{}%'.format(int(100*m)) for m in mod_range]


    modifications = {
            'chain_deletion_prob': mod_range,
            }

    repeats = 10
    match_limit_low = 750
    match_limit = 250

    settings = default_settings()
    settings['cell_freq_max'] = 0.01
    settings['num_cells'] = 1000
    settings['cpw'] = (300,)
    settings['chain_deletion_prob'] = 0.1
    settings['chain_misplacement_prob'] = 0.0
    settings['alpha_sharing_probs'] = None
    settings['beta_sharing_probs'] = None


    all_cm_fdr = {}
    all_rep_fdr = {}

        
    for sn in ['madhype','alphabetr']:

        if sn == 'madhype':
            solver_options = [{
                'threshold':10., # minimum ratio accepted by match_probability
                }]

        elif sn == 'alphabetr':
            solver_options = [{}]


        solvers = [sn]
        all_cm_fdr[sn] = {}
        all_rep_fdr[sn] = {}
        
        for mod,values in modifications.items():

            all_cm_fdr[sn][mod] = []
            all_rep_fdr[sn][mod] = []

            for i,v in enumerate(values): 

                all_results = []

                # iterate across system
                for r in xrange(repeats):

                    specific_settings = copy.copy(settings)

                    specific_settings[mod] = v
                    specific_settings['seed'] = r

                    data,results = simulate_run(solvers,solver_options,**specific_settings)

                    results[0]['fdr_for_cm'] = _get_fdr_for_match_limit(data,results,match_limit)
                    results[0]['fdr_for_rep'] = _get_fdr_for_match_limit(data,results,match_limit_low)

                    print 'FDR (cm):',results[0]['fdr_for_cm']
                    print 'FDR (low cm):',results[0]['fdr_for_rep']

                    all_results += results

                all_cm_fdr[sn][mod].append([results['fdr_for_cm'] for results in all_results])
                all_rep_fdr[sn][mod].append([results['fdr_for_rep'] for results in all_results])
        

    # plot/display settings
    fs = 18
    boxprops = dict(linewidth=3.0,zorder=1)
    meanlineprops = dict(linestyle='-',linewidth=2, color='black', zorder=0)
    plt.rcParams['xtick.labelsize'] = fs-4

    # figure specific properties
    fig,axes = plt.subplots(nrows=2, ncols=2, figsize=(14, 12), sharey=False)
    plt.subplots_adjust(left=0.15,right=0.9,hspace=0.3,wspace=0.5)

    # set border for figure
    for axe in axes:
        for ax in axe:
            [i.set_linewidth(3) for i in ax.spines.itervalues()]

    bp = axes[0][0].boxplot(
            all_cm_fdr['madhype']['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[0][0].set_title('MAD-HYPE',fontweight='bold',fontsize=fs)
    label_figure(axes[0][0],'Chain Deletion Probability','FDR (%)',fs=fs)

    bp = axes[0][1].boxplot(
            all_cm_fdr['alphabetr']['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[0][1].set_title('ALPHABETR',fontweight='bold',fontsize=fs)
    label_figure(axes[0][1],'Chain Deletion Probability','FDR (%)',fs=fs)

    bp = axes[1][0].boxplot(
            all_rep_fdr['madhype']['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[1][0].set_title('MAD-HYPE',fontweight='bold',fontsize=fs)
    label_figure(axes[1][0],'Chain Deletion Probability','FDR (%)',fs=fs)

    bp = axes[1][1].boxplot(
            all_rep_fdr['alphabetr']['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[1][1].set_title('ALPHABETR',fontweight='bold',fontsize=fs)
    label_figure(axes[1][1],'Chain Deletion Probability','FDR (%)',fs=fs)

    plt.savefig('Figure5C.png', format='png', dpi=300)


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
    elif ylabel == 'FDR (%)':
        dy = 0.02
        ax.set_ylim((0. - dy,0.25 + dy))
        ax.set_yticks((0.0,0.05,0.1,0.15,0.2,0.25))
        ax.set_yticklabels(('0%','5%','10%','15%','20%','>25%'),fontsize=fs)

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
            }

if __name__ == "__main__":
    main()





