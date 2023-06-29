'''
Created on Apr 14, 2020

@author: sjjamsa
'''

import copy
import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid, prefix):
    """
    Read particle marker input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "marker_shined/"+prefix+"_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            # Make all read-in datasets 1-d arrays,
            # regardless of their original dimensionality
            d=f[path][key][:]
            out[key] = np.reshape(d,newshape=(d.size,))

    out["ids"] = out["id"]
    del out["id"]
    return out



class mrk(AscotData):
    '''
    A class acting as a superclass for all marker types.
    '''

    def prune(self,keepNmarkers, data=None, probabilities=None):
        '''
        Keep a subset of markers.
        
        Args:
            keepNmarkers : int <br>
                How many markers to keep.
                
            data=None         : dict <br>
                The data from .read() -method.
                If None, will be automatically read.
                
            probabilities=None
                What is the probability to include each marker
                If None, use the weights, (sum normalized to 1.0).
               
        Returns:
            A dictionary as if from .read()
            
            
            
        
        '''
        
        
        if data is None:
            data = self.read()
        else:
            data = copy.deepcopy(data)
        
        n = data['n'][0]
        totalWeight = np.sum(data["weight"])
        
        if probabilities is None:
            probabilities = data["weight"] / totalWeight
        
        
        
        fortune = np.random.choice(np.arange(n),size=(keepNmarkers,),replace=False,p=probabilities)
        
        for k in data.keys():
            if k=='n':
                data[k][0]=keepNmarkers
                continue
            data[k] = data[k][fortune]

        newWeight = np.sum(data["weight"])
        
        data['weight'] *= totalWeight / newWeight

        return data
