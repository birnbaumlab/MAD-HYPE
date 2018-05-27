
"""
This collection of scripts is designed to produce figures that we care about
"""

from figures import fig1
from figures import fig2
from figures import fig4
from figures import fig1s
from figures import fig2s
from figures import fig3s
from main import simulate_system

#------------------------------------------------------------------------------#

settings = {
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
