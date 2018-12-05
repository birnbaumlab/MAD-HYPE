
# standard libraries
import os
from math import ceil,floor,log10

# nonstandard libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

"""

Functions:

    plot_auroc
    plot_comparison
    plot_frequency_estimation
    plot_repertoire

"""


def plot_auroc(cresults,**kwargs):

    label = 'auroc'
    options = _default_plot_options()

    # default options parameters
    new_options = {
            'fs':         18,
            'linewidth':   5,
            'figsize': (8,7),
                  }

    # update options
    options.update(new_options)
    options.update(kwargs)
    _update_options(options,label)

    # get references for fig,ax
    fig,ax = _get_axis(options)

    # local namespace
    fs,linewidth = options['fs'],options['linewidth']

    ### START FIGURE ###

    # get data attributes
    false_limit = cresults['xy'][0][-1]
    true_limit = cresults['xy'][1][-1]
    fdr_limit = min(false_limit,options['fdr_plot']*true_limit)

    # plot FDR line
    plt.plot((0,fdr_limit),(0,fdr_limit/options['fdr_plot']),
            linestyle='--',color='k',linewidth=linewidth)

    # plot main data auroc
    plt.plot(*cresults['xy'],linewidth=linewidth) 

    # add labels
    plt.xlabel('False positives (#)',fontsize=fs)
    plt.ylabel('True positives (#)', fontsize=fs)

    # show plots
    _output(options)

def plot_comparison(cresults,**kwargs):

    label = 'comparison'
    options = _default_plot_options()

    # default options parameters
    new_options = {
                  'pos_color':'black',
                  'mixed1_color':'green',
                  'mixed2_color':'#FFD870',
                  'neg_color':'white',
                  'analysis':('MAD-HYPE','ALPHABETR'),
                  'legend':True,
                  }

    # update options
    options.update(new_options)
    options.update(kwargs)
    _update_options(options,label)

    # get references for fig,ax
    fig,ax = _get_axis(options)

    if len(cresults) != 2:
        print 'Results the wrong length ({})'.format(len(cresults))
        return None 

    # heavy lifting on the data
    assert cresults[0]['freqs'] == cresults[1]['freqs'], \
            'Frequencies between results not identical' 

    ### START FIGURE ###  

    plt.yscale('log') # get y-scale

    total = cresults[0]['total'] # 

    for val,color,l in zip(
            ((0,0),(1,0),(0,1),(1,1)),
            (options['neg_color'],options['mixed1_color'],
                options['mixed2_color'],options['pos_color']),
                ('Neither Correct','MAD-HYPE Correct','ALPHABETR Correct','Both Correct')):
        xs = [i for i,p1,p2 in zip(xrange(total),cresults[0]['pattern'],cresults[1]['pattern']) 
                if p1 == val[0] and p2 == val[1]]
        ys = [f for f,p1,p2 in zip(cresults[0]['freqs'],cresults[0]['pattern'],cresults[1]['pattern'])
                if p1 == val[0] and p2 == val[1]]
        if len(xs) > 0: plt.bar(xs,ys,color=color,width=1,log=True,label=l,edgecolor='none')

    plt.plot(xrange(len(cresults[0]['freqs'])),cresults[0]['freqs'],color='k')

    plt.xlim((0,len(cresults[0]['freqs'])))
    plt.ylim((min(cresults[0]['freqs']),10**ceil(log10(max(cresults[0]['freqs'])-1e-9))))
    ax.set_xticks([0,len(cresults[0]['freqs'])/2,len(cresults[0]['freqs'])])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Clone #',fontweight='bold')
    plt.ylabel('',fontweight='bold',size=12)
    ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(labelsize=20)


    # show plots
    _set_legend(options)
    _output(options)

#------------------------------------------------------------------------------#

def plot_frequency_estimation(cresults,**kwargs):

    label = 'frequency_estimation'
    options = _default_plot_options()

    # default options parameters
    new_options = {
            'fs':          18,
            'linewidth':    3,
            'figsize':  (6,6),
            'colorbar': False,
            'xlim': False,
            'ylim': False,
                  }

    # update options
    options.update(new_options)
    options.update(kwargs)
    _update_options(options,label)

    # get references for fig,ax
    fig,ax = _get_axis(options)

    # local namespace
    fs,linewidth = options['fs'],options['linewidth']

    ### START FIGURE ###

    [i.set_linewidth(linewidth) for i in ax.spines.itervalues()]

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.tick_params(width = 3,length=8,labelsize=18)

    xv = [c[0] for c in cresults['positive_matched_freqs']]
    yv = [c[1] for c in cresults['positive_matched_freqs']]

    positive_confidence = cresults['positive_confidence']
    
    plt.rcParams['image.cmap'] = 'pink'

    if len(cresults['positive_matched_freqs']) != 0:

        min_conf = min(positive_confidence)
        max_conf = max(positive_confidence)
        
        if min_conf == max_conf:
            colors = 'black'
        else:
            colors = [(p - min_conf)/(max_conf - min_conf) for p in positive_confidence]

        if options['xlim']: ax.set_xlim(*options['xlim'])
        else: ax.set_xlim((10**floor(log10(min(xv))),10**ceil(log10(max(xv)))))

        if options['ylim']: ax.set_ylim(*options['ylim'])
        else: ax.set_ylim((10**floor(log10(min(yv))),10**ceil(log10(max(yv)))))

        # plot available species
        sc = plt.scatter(*zip(*cresults['positive_matched_freqs']),
                c=colors,linewidth=0.0,edgecolor='black',s=50) 

    if cresults['negative_matched_freqs']: 
        plt.scatter(*zip(*cresults['negative_matched_freqs']),c='r', marker='x')

    else:
        print 'No matches made, nothing plotted in frequency esimation!'

    # label axes
    plt.xlabel('Clonal Frequency',fontsize=fs,fontweight='bold')
    plt.ylabel('Predicted Frequency',fontsize=fs,fontweight='bold')

    # create colorbar
    if options['colorbar']:
        cbar = plt.colorbar(sc, ticks=[])

    # show plots
    _output(options)

#------------------------------------------------------------------------------#

def plot_repertoire(cresults,**kwargs):

    label = 'repertoire'
    options = _default_plot_options()

    print 'Here'
    # default options parameters
    new_options = {
            'fs':           18,
            'linewidth':     5,
            'figsize':  (10,5),
            'legend':     True,
                  }

    # update options
    options.update(new_options)
    options.update(kwargs)
    _update_options(options,label)

    # get references for fig,ax
    fig,ax = _get_axis(options)

    # local namespace
    fs = options['fs']

    ### START FIGURE ###

    plt.yscale('log')

    for val,color,l in zip(
            (0,1),(options['neg_color'],options['pos_color']),('Incorrect','Correct')):
        xs = [i for i,p in enumerate(cresults['pattern']) if p == val]
        ys = [f for f,p in zip(cresults['freqs'],cresults['pattern']) if p == val]
        if len(xs) > 0: plt.bar(xs,ys,color=color,width=1,log=True,label=l,edgecolor='none')

    plt.plot(xrange(len(cresults['freqs'])),cresults['freqs'],color='k')

    if options['legend'] == True:
        leg = plt.legend()
        leg.get_frame().set_edgecolor('k')

    plt.xlim((0,len(cresults['freqs'])))
    plt.ylim((min(cresults['freqs']),10**ceil(log10(max(cresults['freqs'])-1e-9))))
    ax.set_xticks([0,len(cresults['freqs'])/2,len(cresults['freqs'])])
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Clone #',fontweight='bold',fontsize=fs)
    plt.ylabel('',fontweight='bold',size=12)
    ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(labelsize=20)

    # show plots
    _set_legend(options)
    _output(options)

#------------------------------------------------------------------------------#

def _update_options(options,label):
    if isinstance(options['plot_{}'.format(label)],dict):
        options.update(options['plot_{}'.format(label)])

#------------------------------------------------------------------------------#

def _output(options):
    """ Output """
    # save figure
    if options['save']:
        if '{}' in options['savename']:
            # find a unique filename
            for i in xrange(1000):
                if not os.path.isfile(options['savename'].format(i)):
                    options['fig'].savefig(options['savename'].format(i))
                    print 'Saved figure as:',options['savename'].format(i)
                    break
        else:
            # overwrite
            options['fig'].savefig(options['savename'])


#------------------------------------------------------------------------------#

def _get_axis(options):
    """ Get axes, prioritizing user inputs """
    if not options['fig'] and not options['ax']:
        fig,ax = plt.subplots(figsize = options['figsize'])
        try:
            options['fig'],options['ax'] = plt.gcf(),plt.gca()
        except: # if running on AWS
            pass
    elif not options['ax']:
        fig,ax = options['fig'],plt.gca()
        options['ax'] = plt.gca()
    elif not options['fig']:
        fig,ax = plt.gcf(),options['ax']
        options['fig'] = plt.gcf()
    else:
        fig,ax = options['fig'],options['ax']

    plt.sca(ax)

    return fig,ax

#------------------------------------------------------------------------------#

def _set_legend(options):
    """ Set legend for current figure """
    if options['legend'] == True:
        leg = plt.legend(
                prop={'size': 18},
                edgecolor = 'black',
                frameon = False,
                )
        leg.get_frame().set_edgecolor('k')

#------------------------------------------------------------------------------#

def _default_plot_options():
    """ Options common to all plots """
    return {
            'title':            True, # whether title  is displayed
            'legend':           True, # whether legend is displayed
            'save':             True, # whether plots are saved
            'savename': 'img_{}.png', # whether plots are saved
            'ax':               None, # an axis to plot on
            'fig':              None, # an axis to plot on
            'figsize':         (6,12), # size of newly generated figure
            }

#------------------------------------------------------------------------------#

