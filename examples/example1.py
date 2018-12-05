
from math import log10

# nonstandard libraries
import numpy as np

import matplotlib
matplotlib.rc('text', usetex=True)

import matplotlib.pyplot as plt

from scipy.stats import norm

# Import MAD-HYPE package
import madhype

# Set up run parameters
solvers = ['madhype']
solver_options = [{}] # don't change default parameters
# Uncomment these lines if you want to compare MADHYPE to ALPHABETR
# solvers = ['madhype', 'alphabetr']
# solver_options = [{}, {}] # don't change default parameters

# Set up parameters that apply to all solvers/simulations
general_options = {
        'cpw':(10,250,),
        'num_wells':(48,48,),
        'chain_deletion_prob':0.1,
        'alpha_sharing_probs': 0.0,
        'num_cells': 1000,
        'cell_freq_max': 0.05,
        'beta_sharing_probs': 0.0,
        'alpha_dual_prob': 0.0,
        'beta_dual_prob': 0.0,
        'visual':                     True,
        'plot_repertoire':            False,
        'plot_frequency_estimation':  True,
        'save':                       True,
        'display':                    True,
        'save_to_excel':              False,
        }

# Run MAD-HYPE
data,results = madhype.simulate_run(solvers, solver_options, **general_options)

all_conf = [r[1] for r in results[0]['raw_results']]
all_correct_conf = [r for r in results[0]['positive_confidence']]

# remove the positives from the negatives
conf_dict = {}
for c in all_conf:
    try: conf_dict[c] += 1
    except KeyError: conf_dict[c] = 1

for c in all_correct_conf:
    conf_dict[c] -= 1

all_conf = [c for c,count in conf_dict.items() for _ in xrange(count)]

# Print out results
for solver, result in zip(solvers, results):

    print "{} Results:".format(solver)

    print "  Total # Cells:", result['total']
    print "  Chain pairs identified:", result['positives']
    print "  Chain pairs not identified:", result['negatives']

bins = np.linspace(-5,max(all_correct_conf),15)

fig,ax = plt.subplots(1,1,figsize = (8,6))

n1,_,_ = plt.hist(all_conf,bins, 50, alpha=0.5, label = 'Incorrect Matches', color = 'red')
n2,_,_ = plt.hist(all_correct_conf,bins, 50, alpha=0.5, label = 'Correct Matches', color = 'green')

xlim,ylim = plt.xlim(),plt.ylim()
all_mu,all_sigma = {},{}

for data,counts,label in zip((all_conf,all_correct_conf),(n1,n2),('Incorrect Matches','Correct Matches')):
    (mu, sigma) = norm.fit(data)
    y = norm.pdf(bins,mu,sigma)
    #l = plt.plot(bins, max(counts)*y/max(y), 'k--', linewidth=2)
    #ax.errorbar(x=(mu,),y=(2,),xerr=((-sigma,sigma),),elinewidth = 10)
    all_mu[label] = mu
    all_sigma[label] = sigma

diff = abs(all_mu['Correct Matches'] - all_mu['Incorrect Matches'])

print 'Diff:',diff

#'''
label_diff(
        all_mu['Incorrect Matches'],all_mu['Correct Matches'],
        max(n1),max(n2),
        text = r'{} $\sigma$'.format(
            np.around(diff/np.sum([v**2 for v in all_sigma.values()])**0.5,2)))
#'''

plt.xlim((xlim[0],xlim[1] + 10))
plt.ylim((1,10*ylim[1]))

plt.xlabel(r'Match ratio ($10^{x}$)',fontsize = 24)
plt.ylabel(r'Number of Matches',fontsize = 24)
ax.tick_params(labelsize = 18,length = 10)
ax.tick_params(which = 'minor',labelsize = 18,length = 5)

fig.tight_layout()

[i.set_linewidth(3.0) for i in ax.spines.itervalues()]


plt.yscale('log', nonposy='clip')

plt.legend(fontsize = 18,frameon = False,loc = 7)

plt.savefig("Distribution Demo.svg")
plt.savefig("Distribution Demo.png")
plt.show(block = False)
raw_input('Press enter to close...')
plt.close()









