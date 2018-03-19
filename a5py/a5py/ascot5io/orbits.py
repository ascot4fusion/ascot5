"""
Orbits HDF5 IO module.
"""
import numpy as np
import h5py
import random
import datetime

def read_hdf5(fn, qid):
    """
    Read orbits.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    qid : str
        qid of the run these orbit data correspond to.

    Returns
    -------

    Dictionary storing the orbits that were read.
    """

    f = h5py.File(fn,"r")
    orbits = f["/results/run-"+qid+"/orbits"]

    out = {}

    # Read data from file.
    for orbgroup in orbits:
        out[orbgroup] = {}
        for field in orbits[orbgroup]:
            out[orbgroup][field]           = orbits[orbgroup][field][:]
            out[orbgroup][field + "_unit"] = orbits[orbgroup][field].attrs["unit"]

    # Find how many markers we have and their ids.
    for orbgroup in orbits:
        out[orbgroup]["N"]        = (np.unique(orbits[orbgroup]["id"][:])).size
        out[orbgroup]["uniqueId"] = np.unique(orbits[orbgroup]["id"][:])

    # Sort fields by id (major) and time (minor), both ascending.
    
    for orbgroup in orbits:
        if out[orbgroup]["N"] > 0:
            ind = np.lexsort((out[orbgroup]["id"], out[orbgroup]["time"]))
            for field in orbits[orbgroup]:
                if field[-4:] != "unit":
                    out[orbgroup][field] = out[orbgroup][field][ind]
    
    f.close()
    
    return out


def write_hdf5(fn, orbits, qid):
    """
    Write orbits.

    Unlike most other "write" functions, this one takes dictionary
    as an argument. The dictionary should have exactly the same format
    as given by the "read" function in this module. The reason for this
    is that this function is intended to be used only when combining 
    different HDF5 files into one.

    TODO not compatible with new format

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    orbits : dictionary
        Orbit data to be written in dictionary format.
    qid : int
        Run id these orbits correspond to.
    """

    group = "orbits"
    # Create group
    f = h5py.File(fn, "a")
    if not group in f:
        o = f.create_group(group)
    else:
        o = f[group]

    for orbgroup in orbits:
        path = "orbits/" + orbgroup

        # Remove group if one is already present.
        if path in f:
            del f[path]
        f.create_group(path)

        # TODO Check that inputs are consistent.
        
        
        # Write data to file.
        for field in orbits[orbgroup]:
            if field[-4:] != "unit" and field != "N" and field != "uniqueId":
                d = f.create_dataset(path + "/" + field, data=orbits[orbgroup][field])
                d.attrs["unit"] = orbits[orbgroup][field + "_unit"]

    f.close()
