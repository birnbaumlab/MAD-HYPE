
"""
Attempt to test variable cells per well
"""

# standard libraries

# nonstandard libraries
import numpy as np
import matplotlib.pyplot as plt

# homegrown libraries
from methods import *

def process_results(results,data,*args,**kwargs):

    """
    Processes results dictionary using data object as reference
    """

    # default settings parameters
    settings = {
               }

    # update settings
    for arg in args: settings.update(arg)
    settings.update(kwargs)

    # assertion check
    assert len(results['cells']) == len(results['threshold']), "differing # cells/thresholds"
    
    #
    cells_with_scores = zip(results['cells'],results['threshold'])
    cells_with_scores = sorted(cells_with_scores,key=lambda x: -x[1])

    cells_record = set(data.metadata['cells'])
    total_cells = len(cells_record) 

    # 
    x,y = [0],[0]
    for c in cells_with_scores:
        if c[0] in cells_record:
            x.append(x[-1])
            y.append(y[-1]+1)
        else:
            x.append(x[-1]+1)
            y.append(y[-1])
        if x[-1] == total_cells:
            break

    print 'Total cells:',len(cells_record)

    # 
    plt.plot(x,y)
    plt.show(block=False)
    raw_input('Press enter to close...')
    plt.close()


if __name__ == '__main__':
    test_settings()
