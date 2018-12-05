# standard libraries
from copy import deepcopy
import os
import cPickle as pickle
import gzip 
import csv
import dbm 
from collections import Counter
from datetime import datetime

# nonstandard libraries
import gzip
import matplotlib.pyplot as plt

# homegrown libraries
import madhype
from madhype.defaults import madhype_options as default_options
from madhype.postprocessing.postprocessing import analyze_results

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

    files = subjectXYdata(dirnameX,dirnameY)
    repertoire = catalog_repertoire(*files,overwrite = False)


    data = {
            'options':
            {
                'cpw':     (160000,),
                'num_wells': (96,),
            }
           }

    processing_options = {
            'data_well_threshold':      3,  # minimum number of wells for a chain to appear to be added to dataset
            'data_well_cap':           96,  # maximum number of wells for a chain to appear to be added to dataset
            'overwrite':            False,  # surpresses some inte
            'silent':               False,  # supresses some intermediate input/output
            'subject_well_threshold':   3,  # minimum number of wells for a chain to appear to be added to subject reference 
            'subject_well_cap':        96,  # maximum number of wells for a chain to appear to be added to subject reference 
            'appears_in_reference': False,  # only add a sequence from data if it appears in a reference
            }

    # create a subject reference
    reference = make_reference(repertoire,**processing_options)
    #reference = None

    results = load_howie_data(**dirnameEXP)
    analyze_results(results['Experiment 1'],data,reference=reference)
    analyze_results(results['Experiment 2'],data,reference=reference)
    #raw_input('Waiting..')
    
    options = {
            'reference': reference,
            'visual':        False,
            'max_pairs':    500000,
            'num_wells':     (96,),
            'cpw':         (160000,),
            }

    options.update(processing_options)

    '''
    # run madhype on datasets
    for label,dirname in dirnameEXP.items():
        
        data['well_data'] = data_assignment(label,dirname,reference,**processing_options)

        solvers = ['madhype'] # only run MAD-HYPE
        solver_options =   [{'num_cores':0}] # use default parameters

        # results
        startTime = datetime.now()
        results_madhype = madhype.run(
                data,
                solvers,
                solver_options,
                **options
                )

        print 'MAD-HYPE took {} seconds.\n'.format(datetime.now()-startTime)
    '''

    #combine(results,reference)

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
        
        '''
        print 'States:'
        for xy,frequency in xy_frequency.items():
            print '{}: {}'.format(xy,frequency)
        '''
        print 'FDR cap:'
        for xy,frequency in fdr_count.items():
            print '{}: {}'.format(xy,frequency)
        print 'Total pairs called:',sum(xy_frequency.values())
        
def find_origin(sequence,repertoire):
    """ Find patient origin of sequenece, returns tuple """
    alpha,beta = '',''
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

def load_howie_data(**dirnames):

    howie_results = {}

    for label,dirname in dirnames.items():
        # prepare howie results
        with open(dirname + '/tcr_pairseq_fdr1pct.pairs','r') as f: 
            howie_results[label] = [(((a[0],),(a[2],)),1.-float(a[5]),{'ij':0.0}) for i,a in enumerate(csv.reader(f, dialect="excel-tab")) if i != 0]
    
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

def make_reference(repertoire=None,**kwargs):

    options = {
            'subject_well_threshold': 5,
            'subject_well_cap':      94,
            }

    options.update(kwargs)

    well_threshold = options['subject_well_threshold']
    well_cap =       options['subject_well_cap']

    fname = '{}n-{}m.p'.format(well_threshold,well_cap)

    if not file_check(fname):

        if repertoire == None:

            print 'Need repertoire to construct new references!'
            return None

        print 'Copying repertoire...'
        reference = deepcopy(repertoire)
        print 'Finished copy!'

        # change reference
        change(reference,well_threshold,well_cap)
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
        hits = len([i for i in d])
        if hits <= m and hits >= n: return hits
        return None

def catalog_repertoire(filesX,filesY,overwrite=False,**kwargs):
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
    
    # if there already exists a dictionary that isn't to be overwritten 
    else:
        # TODO: This really isn't necessary
        print 'Loading existing origin dictionary...'
        final_dict = pickle.load(open('./database/origin_dict.p','rb'))
        print 'Finished loading!'
        #for c in 'AB':

    # final counts on chain sequence popularity 
    #for chain_id,seq_dict in final_dict.items():
    #    print 'Total sequences for chain {}: {}'.format(chain_id,len(seq_dict))

    return final_dict

""" File Functions """

def file_check(*fnames):
    """ Checks for presence of pag and dir files, returns boolean """
    files = [os.path.isfile('./database/'+f) for f in fnames]
    val = all(files)
    return val

""" Pickle Functions """

def pickle_load_unique(data_label,my_id,id_label):
    """ Loads pickle by filename """
    if not os.path.isdir('database'):
        os.mkdir('database')
    for file in os.listdir("./database"):
        if file.startswith(id_label) and file.endswith(".p"):
            loaded_id = pickle.load(open('./database/' + file,'rb'))
            if my_id == loaded_id:
                print 'Found matching record: {}'.format(file)
                index = file[file.rindex('_')+1:-2]
                return pickle.load(open('./database/' + data_label + '_{}.p'.format(index),'rb'))
    print 'No matching datasets found!'            
    return None

def pickle_save_unique(my_data,data_label,my_id,id_label):
    """ Saves a pickled version of files """
    if not os.path.isdir('database'):
        os.mkdir('database')
    for i in range(1000):
        if not os.path.isfile('./database/' + data_label + '_{}_.p'.format(i)):
            pickle.dump(my_data,open('./database/' + data_label + '_{}.p'.format(i),'wb'))
            pickle.dump(my_id,open('./database/' + id_label + '_{}.p'.format(i),'wb'))
            return None
    print 'No available names, please delete some files!'

def pickle_load(fname):
    """ Loads pickle by filename """
    return pickle.load(open('./database/' + fname,'rb'))

def pickle_save(my_data,fname):
    """ Saves a pickled version of files """
    if not os.path.isdir('database'):
        os.mkdir('database')
    pickle.dump(my_data,open('./database/' + fname,'wb'))

""" Great """

def data_assignment(label,dirname,reference,**kwargs):
    """ Pulls out data from experiment directory and assigns TCR wells """
    """ Only valid for two subject testing """
    """ Note: threshold mimicked from Howie """

    options = {
            'silent':False,
            'data_well_threshold':     4,
            'data_well_cap':          96,
            'appears_in_reference': True,
            'dirname':           dirname,
            }

    options.update(kwargs) # update dictionary

    # local namespace
    silent = options['silent']
    threshold = options['data_well_threshold']
    cap = options['data_well_cap']
    appears_in_reference = options['appears_in_reference']

    # create DBM dictionaries
    files = {'A':{},'B':{}}

    well_data = pickle_load_unique('well_data',options,'options')

    if well_data:
        print 'Using existing well data savepoint!'
        return well_data

    """ Start unpackaging data, create the four dictionaries """ 
    # initialize dictionaries
    well_data = {'A':[],'B':[]} # will hold lists of lists containing chain indices

    # directory name
    dirfiles = os.listdir(dirname+'/cdna_data/')

    """ Find the file names of all the data, unzip where needed """
    # check across all files, unzip as needed
    for file in dirfiles:
        # check if starting file is gzipped
        if file.endswith('.gz') and not any([file[:-3] in d for d in dirfiles if not d.endswith('.gz')]):
            with gzip.open(dirname+'/cdna_data/'+file,'rb') as f:
                with open(dirname+'/cdna_data/'+file[:-3],'wb') as f_new:
                    f_new.write(f.read(file[:-3]))
            print 'Unzipped {}.'.format(dirname+'/cdna_data/'+file)
            files[file[file.index('TCR')+3]][int(file[file.index('.well')+5: \
                    file.index('.results')])] = dirname+'/cdna_data/'+file[:-3]

        # otherwise store file location
        elif file.endswith('.tsv'):
            files[file[file.index('TCR')+3]][int(file[file.index('.well')+5: \
                    file.index('.results')])] = dirname+'/cdna_data/'+file

    for chain_id,chain_files in files.items(): # iterate across file locations 

        print 'Chain:',chain_id
        print 'Number of chain files:',len(files.items())


        # iterate across repertoires to make independent dictionaries
        for well_id in sorted(chain_files.keys(),key=int):

            # go through file and pull sequences
            with open(files[chain_id][well_id],'rb') as f:

                if not silent: print 'Analyzing well {} for chain {}...'.format(well_id,chain_id)
                well_data[chain_id].append([])

                for i,line in enumerate(csv.reader(f, dialect="excel-tab")):

                    if i == 0: continue # skip header line 

                    #if line[0] == 'TCTCTGCACATTGTGCCCTCCCAGCCTGGAGACTCTGCAGTGTACTTCTGTGCAGCATTAGGTGGTTCTGCAAGGCAACTGACCTTT':

                    '''
                    print 'Found a match in well_id {}'.format(well_id)
                    raw_input('hold')
                    for k in reference[chain_id]['X'].keys():
                        print k
                    '''

                    if not appears_in_reference:
                        well_data[chain_id][-1].append(line[0])
                    elif line[0] in reference[chain_id]['X'] or line[0] in reference[chain_id]['Y']:
                        well_data[chain_id][-1].append(line[0])


                   
    # remove chains that occur less than threshold 
    print 'Adjusting well data for well occurance threshold...'

    chain_keys = {'A':set(k for k,v in Counter(flatten(well_data['A'])).items() if v >= threshold and v <= cap),
                  'B':set(k for k,v in Counter(flatten(well_data['B'])).items() if v >= threshold and v <= cap)}

    print 'Starting well A adjustment...'
    well_data['A'] = [set(compare_lists(j,chain_keys['A'])) for j in well_data['A']]
    print 'Finished well A adjustment!'
    print 'Starting well B adjustment...'
    well_data['B'] = [set(compare_lists(j,chain_keys['B'])) for j in well_data['B']]
    print 'Finished well B adjustment!'

    print 'Key size:'
    print ' > A:',len(chain_keys['A'])
    print ' > B:',len(chain_keys['B'])
    print 'Uniques:' 
    print ' > A:',len(set([i for w in well_data['A'] for i in w]))
    print ' > B:',len(set([i for w in well_data['B'] for i in w]))

    pickle_save_unique(well_data,'well_data',options,'options')
    print 'Finished saving!'

    return well_data


    
class Results:
    """ Quick class to package output data """
    def __init__(self,well_dict,chain_origin,chain_seq):
        # save data to object
        self.well_data = [[a,b] for a,b in zip(well_dict['A'],well_dict['B'])]
        self.chain_origin = chain_origin
        self.chain_seq = chain_seq
        
        # also, pick up some slack and make some files for c++ embedding
        with open('well_data_a.txt','w') as f:
            for w in well_dict['A']: f.write('{}\n'.format(str(w)[1:-1]))
        with open('well_data_b.txt','w') as f:
            for w in well_dict['B']: f.write('{}\n'.format(str(w)[1:-1]))
        with open('uniques_a.txt','w') as f:
            for w in chain_origin['A'].keys(): f.write('{}\n'.format(w))
        with open('uniques_b.txt','w') as f:
            for w in chain_origin['B'].keys(): f.write('{}\n'.format(w))


def flatten(l):
    """ Flatten a list into a 1D list """
    return [item for sublist in l for item in sublist]

def compare_lists(my_list,my_set):
    """ compares two lists of numerals looking for a match """
    return [l for l in my_list if l in my_set]



if __name__ == "__main__":
    rep = main()


