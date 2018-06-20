
import os
from math import ceil,log10
import matplotlib.pyplot as plt

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

    # display figure (with potential hold)
    if options['display'] and options['hold']: 
        plt.show(block=False)
        raw_input('Press enter to close...')
        plt.close()
    elif options['display']:
        plt.show()

def _get_axis(options):
    """ Get axes, prioritizing user inputs """
    if not options['fig'] and not options['ax']:
        fig,ax = plt.subplots(figsize = options['figsize'])
        options['fig'],options['ax'] = plt.gcf(),plt.gca()
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

def _set_legend(options):
    """ Set legend for current figure """
    if options['legend'] == True:
        leg = plt.legend()
        leg.get_frame().set_edgecolor('k')

def _default_plot_options():
    """ Options common to all plots """
    return {
            'title':            True, # whether title  is displayed
            'legend':           True, # whether legend is displayed
            'save':             True, # whether plots are saved
            'savename': 'img_{}.png', # whether plots are saved
            'display':         False, # whether plots are displayed
            'hold':             True, # if plots are displayed, pause for user input
            'ax':               None, # an axis to plot on
            'fig':              None, # an axis to plot on
            'figsize':         (10,5), # size of newly generated figure
            }

def plot_comparison(cresults,**kwargs):

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

    _set_legend(options)

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
    _output(options)


