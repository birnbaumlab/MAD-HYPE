
## MAD-HYPE/examples/figures/figure_4c.py

# nonstandard libraries
import matplotlib.pyplot as plt
from pylab import setp

# Import MAD-HYPE package
import madhype
from madhype.postprocessing.plots import plot_comparison

# library modifications
plt.rcParams["font.family"] = "serif"

def main():

    # Set up run parameters
    solvers = ['madhype','alphabetr']
    solver_options = [{}, {}] # don't change default parameters

    # variants
    cpws = [(10,),(30,)]
    num_simulations = 5

    # Set up parameters that apply to all solvers/simulations
    general_options = {
            'num_wells': (96,),
            }

    fig,ax = plt.subplots(2,1,figsize = (10,10))

    general_options['fig'] = fig

    matches_by_cpw = {
            }

    coverage_by_cpw = {
            }

    for _,cpw in enumerate(cpws):

        print 'Running simulations with {} cpw...'.format(cpw)

        # set the number of cells per well
        general_options['cpw'] = cpw
        matches_by_cpw[cpw] = dict([(s,[]) for s in solvers])
        coverage_by_cpw[cpw] = dict([(s,[]) for s in solvers])

        for index in xrange(num_simulations):

            print 'Starting simulation {}/{}...'.format(index+1,num_simulations)

            # Run MAD-HYPE with default parameters
            data,results = madhype.simulate_run(solvers, solver_options, **general_options)

            # iterate across data
            for method_name,result in zip(solvers,results):
                matches_by_cpw[cpw][method_name].append(result['positives'])
                coverage_by_cpw[cpw][method_name].append(result['frac_repertoire'])

    ### START FIGURE ###

    # settings
    boxprops = dict(linewidth=3.0,zorder=1)
    meanlineprops = dict(linestyle='-',linewidth=2, color='black', zorder=0)
    fs = 18

    # figure specific properties
    fig,axes = plt.subplots(nrows=len(cpws), ncols=1, figsize=(2*len(cpws)+1, 12), sharey=False)
    plt.subplots_adjust(left=0.3,right=0.9,hspace=1.0,wspace=1.0)

    # set border for figure
    for ax in axes:
        [i.set_linewidth(3) for i in ax.spines.itervalues()]

    matches  = [matches_by_cpw[c][s] for s in solvers for c in cpws]
    coverage = [coverage_by_cpw[c][s] for s in solvers for c in cpws]

    labels = ['$N$ = {}'.format(c) for _ in xrange(2) for c in cpws]

    # boxplot matches
    bp = axes[0].boxplot(matches, labels=labels, boxprops=boxprops, meanprops=meanlineprops, widths=0.6, meanline=True, showmeans=True)
    axes[0].plot((2.5,2.5),(0,1000),linestyle='--',color='k')
    setBoxColors(bp)

    axes[0].set_xticks((1.5,3.5))
    axes[0].set_xticklabels(('100','1000'),fontsize=fs)
    axes[0].set_xlabel('Cells/well (#)',fontsize=fs)

    axes[0].set_ylim((0,1000))
    axes[0].set_yticks((0,250,500,750,1000))
    axes[0].set_yticklabels((0,250,500,750,1000),fontsize=fs)
    axes[0].set_ylabel('Clonal Matches (#)',fontsize=fs)

    # boxplot coverage 
    bp = axes[1].boxplot(coverage, labels=labels, boxprops=boxprops, meanprops=meanlineprops, widths=0.6, meanline=True, showmeans=True)
    axes[1].plot((2.5,2.5),(0,1),linestyle='--',color='k')
    setBoxColors(bp)

    axes[1].set_xticks((1.5,3.5))
    axes[1].set_xticklabels(('N = 100','N = 1000'),fontsize=fs)
    axes[1].set_xlabel('Cells/well',fontsize=fs)

    axes[1].set_ylim((0.,1.))
    axes[1].set_yticks((0.,.5,1.))
    axes[1].set_yticklabels(('0%','50%','100%'),fontsize=fs)
    axes[1].set_ylabel('Repertoire Coverage',fontsize=fs)

    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.savefig('fig4C.png', format='png', dpi=200)
    plt.close()


def setBoxColors(bp):
    """ function for setting the colors of the box plots pairs """
    for i in [0,1]:
        setp(bp['boxes'][2*i+0], color='green')
        setp(bp['caps'][4*i+0], color='green')
        setp(bp['caps'][4*i+1], color='green')
        setp(bp['whiskers'][4*i+0], color='green')
        setp(bp['whiskers'][4*i+1], color='green')
        setp(bp['fliers'][2*i+0], color='green')
        setp(bp['medians'][2*i+0], color='green')

        setp(bp['boxes'][2*i+1], color='#FFD870')
        setp(bp['caps'][4*i+2], color='#FFD870')
        setp(bp['caps'][4*i+3], color='#FFD870')
        setp(bp['whiskers'][4*i+2], color='#FFD870')
        setp(bp['whiskers'][4*i+3], color='#FFD870')
        setp(bp['medians'][2*i+1], color='#FFD870')


if __name__ == "__main__":
    main()
