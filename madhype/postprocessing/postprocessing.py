
"""
Attempt to test variable cells per well
"""

# standard libraries
from math import log10,ceil,floor
import os

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import openpyxl

# homegrown libraries
plt.rcParams["font.family"] = "serif"

# library setup
plotly.tools.set_credentials_file(username='Pandyr', api_key='AVy42TUJYGQm0TxLEPMl')

'''
MAIN FUNCTIONS
'''

#------------------------------------------------------------------------------# 

def analyze_results(results,data,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    print 'Analyzing results..'

    # default options parameters
    options = {
              'fdr':0.05,
              'silent':False,
              'reference':None,
              }

    # update options
    options.update(kwargs)

    # if this is not a simulated dataset with a hard set of 
    if not 'cells' in data:

        if options['reference'] == None:
            print 'No cell or subject reference provided, cannot produce meaningful results!'
            return {}

        if not is_reference_valid(options):
            print 'Reference not valid, cannot perform analysis!'
            return {}
        
        print 'Starting analysis using subject reference'

        return get_results_from_subject_reference(results,data,options)

    else:

        return get_results_from_cell_reference(results,data,options)

def get_results_from_subject_reference(results,data,options):

    reference = options['reference']

    cells_with_scores = sorted(results,key=lambda x: -x[1])
    cells_without_scores = [i[0] for i in cells_with_scores]

    reported_matches = [('TCRA','TCRB','p-ratio','TCRA origin','TCRB origin','f_ij','f_i','f_j')]
    total = 0
    x,y,n = 0,0,0
    x1,y1 = [0],[0]
    pattern = []

    for index,c in enumerate(cells_with_scores):

        if index < 1000:
            print index,' : ',c 

        a,b = '',''

        i = c[0][0][0]
        j = c[0][1][0]

        if i in reference['A']['X']: a += 'X'
        if i in reference['A']['Y']: a += 'Y'
        if j in reference['B']['X']: b += 'X'
        if j in reference['B']['Y']: b += 'Y'

        reported_matches.append((i,j,c[1],a,b,c[2]['i'],c[2]['j'],c[2]['ij']))

        total += 1

        if (a == 'X' and b == 'X') or (a == 'Y' and b == 'Y'):
            pattern.append(1)
            y += 1
        elif (a == 'X' and b == 'Y') or (a == 'Y' and b == 'X'):
            pattern.append(-1)
            x1.append(x)
            y1.append(y)
            x += 1
            x1.append(x)
            y1.append(y)
        else:
            pattern.append(0)
            n += 1
            continue 

        if x1[-1] > 2*options['fdr']*y1[-1]:
            print 'FDR met!'
            break

    results = {
              'positives':y1[-1],
              'negatives':x1[-1],
              'neutrals':n,
              'pattern':pattern,
              'total':total,
              'x':x1,
              'y':y1,
              'options':options,
              'matches':reported_matches,
              }

    write_results_to_xslx(results)

    trace0 = go.Scatter(x=x1,y=y1) 
    py.plot([trace0,],filename='Howie FDR')

    print 'Finished plotly stuff!'

    return results

def write_results_to_xslx(results):

    # pull options and matches from results
    options = results.pop('options',{})
    matches = results.pop('matches',{})
    options.pop('reference',{})
    results.pop('reference',{})

    # create new workbook
    wb = openpyxl.Workbook() # open new workbook

    # settings sheet
    _write_dict_to_wb(wb,options,'OPTIONS')
    _write_dict_to_wb(wb,results,'SUMMARY')

    # data sheet
    ws = wb.create_sheet(title='MATCHES')
    for match in matches:
        ws.append(match)

    # save file to unique name
    _unique_wb_save(wb,'results','./results')

def _write_dict_to_wb(wb,settings,sheet_name):
    """ Settings writing """
    ws = wb.create_sheet(title=sheet_name)

    for k,v in settings.items():
        if isinstance(v,dict):
            ws.append((k,'dict->'))
            for k2,v2 in v.items():
                if isinstance(v2,(tuple,list)):
                    ws.append(('->',k2)+tuple(v2))
                else:
                    ws.append(('->',k2,v2))
        elif isinstance(v,(tuple,list)):
            ws.append([k] + list(v))
        else:
            ws.append((k,str(v)))

def _unique_wb_save(wb,name,directory):

    if not os.path.isdir(directory):
        os.mkdir(directory)

    # reference object
    for i in range(1,10000):

        reference_fname = os.path.join(directory,'{}_{}.xlsx'.format(name,str(i).zfill(4)))

        if os.path.isfile(reference_fname): continue

        # try removing starting blank sheet
        try: del wb['Sheet']
        except: pass

        print 'Saving to excel:',reference_fname

        wb.save(reference_fname) # save final workbook

        return

def get_results_from_cell_reference(results,data,options):

    # these are guessed cells
    cells_with_scores = sorted(results,key=lambda x: -x[1])
    cells_without_scores = [i[0] for i in cells_with_scores]

    # these are actual cells
    cells_temp = sorted([(a,b) for a,b in data['cells'].items()],key=lambda x: -x[1])
    cells_record = set([c[0] for c in cells_temp])
    cells_label = [c[0] for c in cells_temp]
    cells_freqs = [c[1] for c in cells_temp]
    total_cells = len(cells_temp)

    # look at error rate (stratified by confidence)
    x1,y1,fdr = [0],[0],0.01
    warning=True # warning if you didn't pass enough points

    for c in cells_with_scores:
        try:
            data['cells'][c[0]]
            x1.append(x1[-1])
            y1.append(y1[-1]+1)
        except KeyError:
            x1.append(x1[-1]+1)
            y1.append(y1[-1])
        if x1[-1] > options['fdr']*y1[-1]:
            warning=False
            break

    if warning:
        print 'WARNING: Number of passed guesses did not meet FDR!'
        print 'Lower the match threshold and you will see better results!'

    total_matches_at_fdr = x1[-1] + y1[-1]

    freqs = [cells_freqs[i] for i,c in enumerate(cells_label) \
             if c in cells_without_scores[:total_matches_at_fdr]]
    frac_repertoire = sum(freqs)
    pattern = [1 if c in cells_without_scores[:total_matches_at_fdr] else 0 \
               for i,c in enumerate(cells_label)]

    # compare actual frequencies to those that are predicted
    positive_matched_freqs = []
    negative_matched_freqs = []
    non_matched_freqs = []
    positive_confidence = []

    pos_freq_dict = dict([(c[0],(c[1],c[2]))
        for c in cells_with_scores[:total_matches_at_fdr]])
    neg_freq_dict = dict([(c[0],(c[1],c[2]))
        for c in cells_with_scores[total_matches_at_fdr:]])
    
    # assign positive and negative matches
    for cell,true_freq in zip(cells_label,cells_freqs):
        try:
            positive_matched_freqs.append((true_freq,pos_freq_dict[cell][1]['ij']))
            positive_confidence.append(np.log10(pos_freq_dict[cell][0]))
        except KeyError:
            try:
                negative_matched_freqs.append((true_freq,neg_freq_dict[cell][1]['ij']))
            except KeyError:
                non_matched_freqs.append((true_freq,true_freq))
                #print '{} did not show up!'.format(cell)

    positive_confidence = [(p - min(positive_confidence))/
            (max(positive_confidence) - min(positive_confidence)) 
            for p in positive_confidence]

    results = {
              'positive_matched_freqs':positive_matched_freqs,
              'negative_matched_freqs':negative_matched_freqs,
              'non_matched_freqs':non_matched_freqs,
              'positive_confidence':positive_confidence,
              'pattern':pattern,
              'frac_repertoire':frac_repertoire,
              'positives':sum(pattern),
              'negatives':len(pattern) - sum(pattern),
              'total':len(pattern),
              'freqs':cells_freqs,
              'xy':[x1,y1]
              }

    if True:#not options['silent']:
        # display characteristics of the data

        print 'Positives:',results['positives']
        print 'Negatives:',results['negatives']
        print 'Fraction of repertoire:',results['frac_repertoire']
        print 'Total clones:',results['total']

    return results

#------------------------------------------------------------------------------# 

def visualize_results(results,data,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    # default options parameters
    options = {
            'fdr':0.01,
            'fdr_plot':0.01,
            'pos_color':'black',
            'neg_color':'white',
            'legend':True,
            'silent':False,
            'visual_block':True
               }

    # update options
    options.update(kwargs)

    # heavy lifting on the data
    cresults = analyze_results(results,data,**options)

    # 

    ### AUROC FIGURE ###

    ax = plt.figure().gca() # initialize figure

    linewidth = 5
    fs = 18

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
    plt.xlabel('False positives (#)')
    plt.ylabel('True positives (#)')
    plt.savefig('routine1.png', format='png', dpi=300)

    ### FREQUENCY ESTIMATION FIGURE ###

    ax = plt.figure(figsize=(9,6)).gca()
    plt.subplots_adjust(left=0.3,right=0.9,top=0.9,bottom=0.2,hspace=1.0,wspace=1.0)
    [i.set_linewidth(3) for i in ax.spines.itervalues()]

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim((7e-5,1.5e-2))

    plt.rcParams['image.cmap'] = 'pink'

    # plot available species
    sc = plt.scatter(*zip(*cresults['positive_matched_freqs']),c=cresults['positive_confidence'],linewidth=0.0,edgecolor='black') 

    if cresults['negative_matched_freqs']: plt.scatter(*zip(*cresults['negative_matched_freqs']),c='r', marker='x')
    #if cresults['non_matched_freqs']: plt.scatter(*zip(*cresults['non_matched_freqs']),c='k', marker='x')

    # label axes
    plt.xlabel('Clonal Frequency',fontsize=fs,fontweight='bold')
    plt.ylabel('Predicted Frequency',fontsize=fs,fontweight='bold')

    # add colorbar
    #plt.colorbar(sc)
    cbar = plt.colorbar(sc, ticks=[])
    plt.savefig('routine2.png', format='png', dpi=300)

    ### REPERTOIRE DISPLAY FIGURE ###

    fig,ax = plt.subplots(figsize=(10,5))

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
    plt.xlabel('Clone #',fontweight='bold')
    plt.ylabel('',fontweight='bold',size=12)
    ax.axes.get_xaxis().set_visible(False)
    plt.tick_params(labelsize=20)
    plt.savefig('routine3.png', format='png', dpi=300)

    # show plots
    plt.show(block=False)
    if options['visual_block']:
        raw_input('Press enter to close...')
        plt.close()

    return cresults

#------------------------------------------------------------------------------# 

def is_reference_valid(options):

    reference = options['reference'] # transfer to local namespace

    passing = True

    # if this is not a simulated dataset with a hard set of 
    if reference == None:
        print 'No cell or subject reference provided, cannot produce meaningful results!'
        return False  

    # check to make sure the right dict keys exist
    for s in 'AB':

        if not s in reference:
            print '{} not in reference dict!'.format(s)
            passing = False
            continue

        for c in 'XY':
            if not c in reference[s]:
                print '{} not in reference[{}] subdict!'.format(c,s)
                passing = False

    return passing
    




    check_subject_reference(options)


def compare_results(cresults,**kwargs):

    # default options parameters
    options = {
            'fdr':0.01,
            'fdr_plot':0.01,
            'pos_color':'black',
            'mixed1_color':'green',
            'mixed2_color':'#FFD870',
            'neg_color':'white',
            'analysis':('MAD-HYPE','ALPHABETR'),
            'legend':True,
            'silent':False
               }

    # update options
    options.update(kwargs)

    if len(cresults) != 2:
        print 'Results the wrong length ({})'.format(len(cresults))
        return None 

    # heavy lifting on the data
    assert cresults[0]['freqs'] == cresults[1]['freqs'],'Frequencies between results not identical' 

    ### REPERTOIRE DISPLAY FIGURE ###

    fig,ax = plt.subplots(figsize=(10,5))

    plt.yscale('log')

    total = cresults[0]['total']    

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

    if options['legend'] == True:
        leg = plt.legend()
        leg.get_frame().set_edgecolor('k')


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
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()
