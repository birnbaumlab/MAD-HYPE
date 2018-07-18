


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


def _get_fdr_for_clonal_matches(data,results,clone_match_threshold):

    positives,negatives = 0,0

    for i,r in enumerate(results[0]['raw_results']):

        if r[0] in data['cells']:
            positives += 1
        else:
            negatives += 1

        if clone_match_threshold <= positives:
            break

        if i == len(results[0]['raw_results']) - 1:
            print 'DID NOT CAP'
            return 0.25

    return float(negatives)/(positives + negatives)

def main(*args,**kwargs):

    mod_range = [.0,.1,.2,.3]
    #mod_range = [.0,.05,.1,.15,.2,.25,.3]
    labels = ['{}%'.format(int(100*m)) for m in mod_range]


    modifications = {
            'chain_deletion_prob': mod_range,
            }

    repeats = 1 
    clone_match_threshold = 200

    settings = default_settings()
    settings['cell_freq_max'] = 0.01
    settings['num_cells'] = 400
    settings['cpw'] = (100,)
    settings['chain_deletion_prob'] = 0.1
    settings['chain_misplacement_prob'] = 0.0


    all_fdr = {}

    solver_options = [{
        'threshold':10., # minimum ratio accepted by match_probability
        }]
        
    for sn in ['madhype','alphabetr']:

        solvers = [sn]
        all_fdr[sn] = {}
        
        for mod,values in modifications.items():

            all_fdr[sn][mod] = []

            for i,v in enumerate(values): 

                all_results = []

                # iterate across system
                for r in xrange(repeats):

                    specific_settings = copy.copy(settings)

                    specific_settings[mod] = v
                    specific_settings['seed'] = r

                    data,results = simulate_run(solvers,solver_options,**specific_settings)

                    results[0]['fdr_for_threshold'] = _get_fdr_for_clonal_matches(data,results,clone_match_threshold)

                    print 'FDR:',results[0]['fdr_for_threshold']

                    all_results += results

                all_fdr[sn][mod].append([results['fdr_for_threshold'] for results in all_results])
        

    # plot/display settings
    fs = 18
    boxprops = dict(linewidth=3.0,zorder=1)
    meanlineprops = dict(linestyle='-',linewidth=2, color='black', zorder=0)
    plt.rcParams['xtick.labelsize'] = fs-4

    # figure specific properties
    fig,axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 12), sharey=False)
    plt.subplots_adjust(left=0.15,right=0.9,hspace=0.3,wspace=0.5)

    # set border for figure
    for axe in axes:
        for ax in axe:
            [i.set_linewidth(3) for i in ax.spines.itervalues()]

    bp = axes[0][0].boxplot(
            all_fdr['madhype']['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[0][0].set_title('MAD-HYPE, 50% Clonal Matches',fontweight='bold',fontsize=fs)
    label_figure(axes[0][0],'Chain Deletion Probability','FDR (%)',fs=fs)

    bp = axes[0][1].boxplot(
            all_fdr['alphabetr']['chain_deletion_prob'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    axes[0][1].set_title('MAD-HYPE, 50% Clonal Matches',fontweight='bold',fontsize=fs)
    label_figure(axes[0][1],'Chain Deletion Probability','FDR (%)',fs=fs)


    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig_S5.png', format='png', dpi=300)
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
    elif ylabel == 'FDR (%)':
        ax.set_ylim((0,0.25))
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
            'visual_block':False,
            'compare':False
            }

if __name__ == "__main__":
    main()





