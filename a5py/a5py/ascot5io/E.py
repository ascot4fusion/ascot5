'''
Created on Apr 16, 2020

@author: sjjamsa
'''
from ._iohelpers.treedata import DataGroup

class E(DataGroup):
    '''
    A parent class for all Electric field classes.
    '''

    def write_dummy(self,fn,desc=None):
        from a5py.ascot5io.E_TC import write_hdf5_dummy
        
        if desc is not None:
            return write_hdf5_dummy(fn=fn, desc=desc)
        
        return write_hdf5_dummy(fn=fn)