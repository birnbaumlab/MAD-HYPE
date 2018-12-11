general_options = {
    # experimental design
    'num_wells':(96,),
    'cpw':(100,),

    # simulated repertoire
    'num_cells':1000,
#    'seed':1,
    'cell_freq_distro': 'power-law',
    'cell_freq_constant':       2.0,
    'cell_freq_max':            0.01,
    'chain_misplacement_prob':  0.0, # TODO: add functionality
    'chain_deletion_prob':      0.1,
    'alpha_dual_prob':          0.0,
    'beta_dual_prob':           0.0,
    'alpha_sharing_probs':      0.0,
    'beta_sharing_probs':       0.0,

    # postprocessing
    'visual':                    False,
    'save_to_excel':             True,
    'plot_auroc':                False,
    'plot_comparison':            False,
    'plot_repertoire':           False,
    'plot_frequency_estimation': False,

    # general plotting settings 
    'title':            True, # whether title  is displayed
    'legend':           True, # whether legend is displayed
    'save':            False, # whether plots are saved
    'savename': 'img_{}.png', # where plots are saved
    'ax':               None, # an axis to plot on
    'fig':              None, # an axis to plot on
    'figsize':         (8,6), # size of newly generated figure
    'fs':19,
    'linewidth':5,
    'fdr_plot': 0.01,
    'pos_color': 'black',
    'neg_color': 'white',
    'mixed1_color': 'green',
    'mixed2_color': '#FFD870',
    'analysis': ('MAD-HYPE','ALPHABETR'),
    'colorbar': True,
    'xlim': False,
    'ylim': False,
    # plot-specific settings (overrides general settings above)
    'plot_auroc_options': {
        'fs': 18,
        'linewidth': 5,
        'figsize': (8,7),
        'fdr_plot': 0.01
    },
    'plot_comparison_options': {
        'pos_color': 'black',
        'mixed1_color': 'green',
        'mixed2_color': '#FFD870',
        'neg_color': 'white',
        'analysis': ('MAD-HYPE', 'ALPHABETR'),
        'legend': True,
    },
    'plot_frequency_estimation_options': {
        'fs': 18,
        'linewidth': 3,
        'figsize': (6,6),
        'colorbar': False,
        'xlim': False,
        'ylim': False
    },
    'plot_repertoire_options': {
        'fs': 18,
        'linewidth': 5,
        'figsize': (10,5),
        'pos_color': 'black',
        'neg_color': 'white'
    },

    # other postprocessing options
    'max_pairs':               1000000,
    'fdr':                        0.01,
    'reference':                  None,

    # misc
    'silent':False
}

madhype_options = {
    # analysis constants
    'threshold':0.1, # minimum ratio accepted by match_probability
    'fdr':0.01, # acceptable fdr (cuts off matches, sets filter)
    'prior_alpha':2.0, # prior for clonal frequency
    'prior_match':None, # prior for clonal match ( <= 1.0 )
    'num_cores':0, # 0 -> max_core usage

    # visual cues
    'silent':False
}

alphabetr_options = {
    'iters':100,
    'pair_threshold':-1,
    'dual_clones': False,
    'estimate_frequencies': True,

    'num_cores':0, # 0 -> max_core usage

    'silent':False,
}
