"""
Initialize groups for storing ASCOT5 data.
"""
import numpy as np
import h5py
import random
import datetime

def setmetadata(group):
    """
    Set metadata to HDF5 group
    
    Metdata consists of unique (not quaranteed but almost 
    certainly) id and date when field was created.

    Parameters
    ----------

    group : HDF5 group
        Group where metadata will be written.
    """
    qid = random.getrandbits(32)
    group.attrs["qid"]  = np.int32(qid)
    group.attrs["date"] = np.string_(datetime.datetime.now())

    
def replacegroup(hdf5file, path):
    """
    Create a new group or replace old with new.

    Parameters
    ----------
    
    hdf5file : HDF5 file
        Target HDF5 file.
    path : str
        Path to group e.g. "bfield/B_2D".
    """
    
    # Remove group if one is already present.
    if path in hdf5file:
        del hdf5file[path]
    hdf5file.create_group(path)

def setgrouptype(hdf5file, group, type_):
    """
    Create a new master group and set its type.

    Parameters
    ----------
    
    hdf5file : HDF5 file
        Target HDF5 file.
    path : str
        Path to group e.g. "bfield".
    type : str
        Which subfield should be used e.g. "B_2D".
    """

    # Create group and set the type to this one.
    if not group in hdf5file:
        o = hdf5file.create_group(group)
        o.attrs["type"] = np.string_(type_)
    else:
        o = hdf5file[group]
        del o.attrs["type"]
        o.attrs["type"] = np.string_(type_)
