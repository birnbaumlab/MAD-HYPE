general_options = {
    # experimental design
    'num_wells':(96,),
    'cpw':(50,),

    # simulated repertoire
    'num_cells':1000,
#    'seed':1,
    'cell_freq_distro': 'power-law',
    'cell_freq_constant':       1.0,
    'cell_freq_max':            0.01,
    'chain_misplacement_prob':  0.0, # TODO: add functionality
    'chain_deletion_prob':      0.1,
    'alpha_dual_prob':          0.0,
    'beta_dual_prob':           0.0,
    'alpha_sharing_probs':      0.0,
    'beta_sharing_probs':       0.0,

    # postprocessing
    'visual':True,
    'compare':False,
    'max_pairs':1000000,
    'fdr': 0.01,
    'reference': None,

    # visualization
    'fdr_plot': 0.01,
    'pos_color': 'black',
    'neg_color': 'white',
    'legend': True,
    'visual_block': True,

    # misc
    'silent':False
}

madhype_options = {
    # analysis constants
    'threshold':100., # minimum ratio accepted by match_probability
    'fdr':0.01, # acceptable fdr (cuts off matches, sets filter)
    'prior_alpha':1.0, # prior for clonal frequency
    'prior_match':1.0, # prior for clonal match ( <= 1.0 )
    'num_cores':0, # 0 -> max_core usage

    # visual cues
    'silent':False
}

alphabetr_options = {
    'iters':100,
    'pair_threshold':0.9,

    'num_cores':0, # 0 -> max_core usage

    'silent':False,
}
