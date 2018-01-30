
# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
from main import simulate_system 

# library modifications
plt.rcParams["font.family"] = "serif"
plt.rcParams['xtick.labelsize'] = 14 
plt.rcParams['xtick.labelsize'] = 14

def main(*args,**kwargs):

    modifications = {
            'prior_alpha':[.1,0.5,1.0,2.0,10]
            }

    repeats =  10

    settings = default_settings()
    settings['cell_freq_max'] = 0.01
    settings['num_cells'] = 1000
    settings['cpw'] = (100,)
    settings['chain_deletion_prob'] = 0.1
    settings['chain_misplacement_prob'] = 0.0

    all_coverage = {}
    all_matches = {}

    #
    for mod,values in modifications.items():

        all_results = []
        all_coverage[mod] = []
        all_matches[mod] = []

        for i,v in enumerate(values): 

            # iterate across system
            for r in xrange(repeats):
                all_results += simulate_system(settings,**{mod:v,'seed':r})

            all_coverage[mod].append([results['frac_repertoire'] for results in all_results])
            all_matches[mod].append([results['positives'] for results in all_results])


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

    labels = [
            '0.1',
            '0.5',
            '1',
            '2',
            '10'
            ]
    
    bp = axes[0][0].boxplot(
            all_matches['prior_alpha'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[0][0],r'Prior $\alpha$','Clonal Matches (#)',fs=fs)

    bp = axes[0][1].boxplot(
            all_coverage['prior_alpha'], 
            labels=labels, 
            boxprops=boxprops, 
            meanprops=meanlineprops, 
            widths=0.6, 
            meanline=True, 
            showmeans=True
            )

    label_figure(axes[0][1],r'Prior $\alpha$','Repertoire Coverage',fs=fs)


    ##### HEATMAPS #####

    print 'Start heatmaps...'

    settings = default_settings()
    settings['chain_deletion_prob'] = 0.1
    settings['chain_misplacement_prob'] = 0.0

    num_wells_range = [(48,),(96,)] 
    cpw_range = np.logspace(0,4,41,dtype=int)
    freq_range = np.logspace(-4,-1,13)
    #cpw_range = np.logspace(0,3,4,dtype=int)
    #freq_range = np.logspace(-2,-1,3)

    repeats =  3
    fs = 18

    id_map = np.zeros((len(freq_range),len(cpw_range)))

    for plot_ind,num_wells in enumerate(num_wells_range):
        for i,freq in enumerate(freq_range):
            for j,cpw in enumerate(cpw_range):
                print 'Starting f = {} / cpw = {}...'.format(freq,cpw)
                val = []
                for r in xrange(repeats):
                    results = simulate_system(
                            settings,
                            num_wells=num_wells,
                            cpw=(cpw,),
                            cell_freq_max=0.0,
                            num_cells=int(1./freq)
                            )
                    val.append(results[0]['frac_repertoire'])
                id_map[i][j] = np.mean(val)

        axes[1][plot_ind].imshow(id_map,interpolation='nearest')
        axes[1][plot_ind].set_aspect(aspect='auto')


        # X axis
        c_labels = [v for v in [1,10,100,1000,10000] if v <= max(cpw_range) and v >= min(cpw_range)]
        c_inds = [min(range(len(cpw_range)), key=lambda i: abs(cpw_range[i]-v)) for v in c_labels]
        axes[1][plot_ind].set_xticks(c_inds)
        axes[1][plot_ind].set_xticklabels(c_labels)

        # Y axis
        c_labels = [v for v in [1e-4,1e-3,1e-2,1e-1] if v <= max(freq_range) and v >= min(freq_range)]
        c_inds = [min(range(len(freq_range)), key=lambda i: abs(freq_range[i]-v)) for v in c_labels]
        axes[1][plot_ind].set_yticks(c_inds)
        axes[1][plot_ind].set_yticklabels(c_labels)


        plt.title('Identification of clones with {} wells'.format(num_wells[0]))
        axes[1][plot_ind].set_xlabel('Cells/well',fontsize=fs)
        axes[1][plot_ind].set_ylabel('Clonal Frequency',fontsize=fs)

    # plot figure
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig3S.png', format='png', dpi=300)
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





