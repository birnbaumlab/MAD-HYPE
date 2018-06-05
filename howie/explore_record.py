

"""
        ### analysis on presented data
        filesX,filesY = subjectXYdata(dirnameX,dirnameY) # returns dictionaries
        catalog_repertoire(filesX,filesY,overwrite=False) 
        data = data_assignment(dirname_exp,threshold=(5,91),overwrite=False,silent=False) # no save due to memory
"""


# standard libraries
from copy import deepcopy
import os
import cPickle as pickle
import gzip 
import csv
import dbm 
from collections import Counter

# nonstandard libraries
import gzip
import matplotlib.pyplot as plt
import plotly
import plotly.plotly as py
import plotly.graph_objs as go

# homegrown libraries


# library setup

plotly.tools.set_credentials_file(username='Pandyr', api_key='AVy42TUJYGQm0TxLEPMl')

def main():
    """ Runs collection of main function """

    repertoire = None

    dirnameDATA = 'pairseq-data'
    dirnameX    = dirnameDATA + '/subjectX' 
    dirnameY    = dirnameDATA + '/subjectY' 

    dirnameEXP = {
            'Experiment 1': dirnameDATA + '/experiment1',
            'Experiment 2': dirnameDATA + '/experiment2',
            }

    #check_unzip_data(dirnameDATA)

    files = subjectXYdata(dirnameX,dirnameY)

    repertoire = catalog_repertoire(*files,overwrite = False)

    for i in range(1,7):
        reference = make_reference(repertoire,well_threshold = i)
    
    results = load_howie_data(**dirnameEXP)

    combine(results,reference)

    print 'Finished!'

    return repertoire

def combine(results,repertoire):
    
    for label,seqs in results.items():
        print 'Starting analysis on {}...'.format(label)
        seqs.sort(key=lambda x: x[2])
        x,y = [0],[0]
        xys = []

        fdr_limit = False
        fdr_threshold = 0.01
        fdr_count = {
                'X':    0,
                'Y':    0,
                'XY':   0,
                }

        for i,seq in enumerate(seqs):

            xy = find_origin(seq,repertoire)
            state = check_origin(xy)
            xys.append(xy)

            if state in ('X','Y'):
                if state == 'X' and not fdr_limit: fdr_count['X'] += 1
                if state == 'Y' and not fdr_limit: fdr_count['Y'] += 1
                x.append(x[-1])
                y.append(y[-1]+1)
            elif state == False:
                if not fdr_limit: fdr_count['XY'] += 1
                if fdr_threshold*(sum(fdr_count.values())) < 2*fdr_count['XY']: fdr_limit = True
                x.append(x[-1]+1)
                y.append(y[-1])
            elif state == None:
                pass

        xy_frequency = Counter(xys)
        
        print 'States:'
        for xy,frequency in xy_frequency.items():
            print '{}: {}'.format(xy,frequency)
        print 'FDR cap:'
        for xy,frequency in fdr_count.items():
            print '{}: {}'.format(xy,frequency)
        print 'Total pairs called:',sum(xy_frequency.values())
        

        trace0 = go.Scatter(x=x,y=y) 
        py.plot([trace0,],filename='Howie FDR')

        raw_input('Press enter to close..')

def find_origin(sequence,repertoire):
    """ Find patient origin of sequenece, returns tuple """
    print sequence[0]
    alpha,beta = '',''
    print type(repertoire['A'])
    print len(repertoire['A'])
    if sequence[0] in repertoire['A']['X']: alpha += 'X'
    if sequence[0] in repertoire['A']['Y']: alpha += 'Y'
    if sequence[1] in repertoire['B']['X']: beta  += 'X'
    if sequence[1] in repertoire['B']['Y']: beta  += 'Y'
    return (_order(alpha),_order(beta)) 

def _order(state):
    if state == 'YX': return 'XY'
    return state

def check_origin(origin):
    """ Sorting function for origin """
    if origin == ('X','X'):
        return 'X' 
    if origin == ('Y','Y'):
        return 'Y' 
    elif '' in origin or 'XY' in origin:
        return None 
    else:
        return False

def check_origin_howie(origin):
    """ Sorting function for origin """
    if origin == ('X','X') or origin == ('XY','X'):
        return 'X' 
    if origin == ('Y','Y') or origin == ('XY','Y'):
        return 'Y' 
    elif '' in origin or 'XY' in origin:
        return None 
    else:
        return False

def load_howie_data(**dirnames):

    howie_results = {}

    for label,dirname in dirnames.items():
        # prepare howie results
        with open(dirname + '/tcr_pairseq_fdr1pct.pairs','r') as f: 
            howie_results[label] = [(a[0],a[2],float(a[5])) for i,a in enumerate(csv.reader(f, dialect="excel-tab")) if i != 0]
    
    return howie_results

def check_unzip_data(zip_fname):
    """ Checked for unzipped data from Howie """

    dirnameX,dirnameY = './data/howie/subjectX','./data/howie/subjectY'

    if os.path.isdir(zip_fname):
        return 

    if os.path.isfile(zip_fname+'.tgz'):
        raise OSError('Pairseq file found but needs to be extracted!')
                     
    else:
        raise OSError('Pairseq file not found as: pairseq-data.tgz') 


def subjectXYdata(dirnameX,dirnameY):
    """ Pulls out file names corresponding the directories submitted, subjects X/Y """
    # initialize dictionaries
    subject_x_files,subject_y_files = {'gdna':{},'cdna':{}},{'gdna':{},'cdna':{}}
    
    # iterate across dirnames and file dictionaries
    for dirname,files in zip([dirnameX,dirnameY],[subject_x_files,subject_y_files]):
        for data_type in ['cdna']:
            # story directory contents
            dirfiles = os.listdir(dirname+'/{}_data'.format(data_type))
            # check across files for file in dirfiles:
            for file in dirfiles:
                # check if starting file is gzipped
                if file.endswith('.gz') and not any([file[:-3] in d for d in dirfiles if not d.endswith('.gz')]):
                    with gzip.open(dirname+'/{}_data/'.format(data_type)+file,'rb') as f:
                        with open(dirname+'/{}_data/'.format(data_type)+file[:-3],'wb') as f_new:
                            f_new.write(f.read(file[:-3]))
                    print 'Unzipped {}.'.format(dirname+'/{}_data/'.format(data_type)+file)
                    files[data_type][file[file.index('TCR')+3] + file[file.index('.well')+5:file.index('.results')]] = \
                            dirname+'/{}_data/'.format(data_type)+file[:-3]

                # otherwise store file location
                elif file.endswith('.tsv'):
                    files[data_type][file[file.index('TCR')+3] + file[file.index('.well')+5:file.index('.results')]] = \
                            dirname+'/{}_data/'.format(data_type)+file
            
    return subject_x_files,subject_y_files

def make_reference(repertoire=None,count_threshold = 1,well_threshold = 5):

    fname = '{}n-{}m.p'.format(count_threshold,well_threshold)

    if not file_check(fname):

        if repertoire == None:

            print 'Need repertoire to construct new references!'
            return None

        print 'Copying repertoire...'
        reference = deepcopy(repertoire)
        print 'Finished copy!'

        # change reference
        change(reference,count_threshold,well_threshold)
        # dump contents as pickle
        pickle_save(reference,fname)

    else:
        print 'Loading reference...'
        reference = pickle_load(fname)

    print 'Done!'

    return reference
    

def change(d,n,m):
    """ Recursively modify dictionary entries """
    if isinstance(d,dict):
        for k,v in d.items():
            d[k] = change(v,n,m) # recursively recall 
            if d[k] == None: d.pop(k,None) # remove key if none referenced
        return d
    else:
        hits = len([i for i in d if i >= n])
        if hits >= m: return hits
        return None

def catalog_repertoire(filesX,filesY,overwrite=False):
    """ Creates a dict that takes each unique sequence and assigns ownership between subject X/Y """
    """ This analysis is very slow, so I'm pulling as many tricks out as I can """

    fnames = ('origin_dict.p',)

    final_dict = {'A':{},'B':{}}

    existing_datasets = file_check(*fnames)

    if overwrite == True or existing_datasets == False:

        origin_dict = {} 
        final_dict = {
                'A':{
                    'X':{},
                    'Y':{},
                    }, 
                'B':{
                    'X':{},
                    'Y':{},
                    },
                }

        # iterate across repertoires to make independent dictionaries
        for patient_id,files in zip(['X','Y'],[filesX,filesY]):
            origin_dict[patient_id] = {}
            for chain_id in ['A','B']: 
                found,keep = set(),[] # create some storage variables
                add,app = found.add,keep.append
                # return 
                for file_id,file in files['cdna'].items(): 
                    if not chain_id in file_id: continue # catch to weed out nonmatching chains
                    print 'Processing sample {} for patient {}'.format(file_id,patient_id) 
                    # go through file and pull sequences
                    seq_freq = {} 
                    with open(file,'rb') as f:
                        for line in csv.reader(f, dialect="excel-tab"):
                            if line[0] not in seq_freq:
                                seq_freq[line[0]] = 1
                            else:
                                seq_freq[line[0]] += 1
                    for seq,freq in seq_freq.items():
                        if seq in final_dict[chain_id][patient_id]:
                            final_dict[chain_id][patient_id][seq].append(freq)
                        else:    
                            final_dict[chain_id][patient_id][seq] = [freq]
                        

        print 'Writing...'
        pickle.dump(final_dict,open('./database/origin_dict.p','wb'),protocol=pickle.HIGHEST_PROTOCOL) 
        print 'Finished writing!'
        # DBM save
        #db_save(final_dict['A'],'origin_dict_A') 
        #db_save(final_dict['B'],'origin_dict_B') 
    
    # if there already exists a dictionary that isn't to be overwritten 
    else:
        # TODO: This really isn't necessary
        print 'Loading existing origin dictionary...'
        final_dict = pickle.load(open('./database/origin_dict.p','rb'))
        print 'Finished loading!'
        #for c in 'AB':
        #    final_dict[c] = db_load(fname_dict[c])

    # final counts on chain sequence popularity 
    #for chain_id,seq_dict in final_dict.items():
    #    print 'Total sequences for chain {}: {}'.format(chain_id,len(seq_dict))

    return final_dict

""" File Functions """

def file_check(*fnames):
    """ Checks for presence of pag and dir files, returns boolean """
    db_files = [os.path.isfile('./database/'+f) for f in fnames]
    val = all(db_files)
    return val

""" Pickle Functions """

def pickle_load(fname):
    """ Loads pickle by filename """
    return pickle.load(open('./database/' + fname,'rb'))

def pickle_save(fname):
    """ Saves a pickled version of files """
    if not os.path.isdir('database'):
        os.mkdir('database')
    pickle.dump(final_dict,open('./database/' + fname,'wb'))


""" DBM Functions """

def db_load(fname):
    return dbm.open('./database/'+os.path.splitext(fname)[0],'c')

def db_save(my_dict,fname):
    if not os.path.isdir('database'):
        os.mkdir('database')
    db = dbm.open('database/'+fname,'c')
    for k,v in my_dict.items():
        db[str(k)] = str(v) 
    db.close()


if __name__ == "__main__":
    rep = main()


