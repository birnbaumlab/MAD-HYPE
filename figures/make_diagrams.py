
"""
This collection of scripts is designed to produce figures that we care about
"""

import fig1
import fig2
import fig4
import fig1s
import fig2s
import fig3s

#------------------------------------------------------------------------------#

options = {
        'subplots':[1,2],
        'fdr':0.03,
        'num_cells':1000,
        'cell_freq_max':0.01,
        'num_wells':(96,),
        'cpw':(100,),
        'seed':1,
        # visual cues
        'silent':False,
        'visual':True
        }

fig1.main()
#fig2.main()
#fig4.main()

#fig1s.main()
#fig2s.main()
#fig3s.main()


#------------------------------------------------------------------------------#
