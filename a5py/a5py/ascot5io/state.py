"""
State HDF5 IO module.

File: state.py
"""

import numpy as np
import h5py
import random
import datetime

def read_hdf5(fn, qid, read="all"):
    """
    Read all or specified states.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    read : string list, optional
        Which states are read e.g. "inistate". Default is all.

    Returns
    -------

    Dictionary storing the states that were read.
    """

    path = "results/run-"+qid

    with h5py.File(fn,"r") as f:
        out = {}
        if read == "all":
            read = ["inistate", "endstate"]

        # Read data from file
        for statename in read:
            out[statename] = {}
            if statename in f[path]:
                for field in f[path][statename]:
                    out[statename][field]           = f[path][statename][field][:]
                    out[statename][field + "_unit"] = f[path][statename][field].attrs["unit"]
            else:
                print("Warning: State " + statename + " does not exists.")

        # TODO Parse endconditions.

        # Find number of markers and check that no markers share same id
        # (which they shouldn't).
        for state in out:
            out[state]["N"] = np.unique(out[state]["id"]).size
            if out[state]["N"] != out[state]["id"].size:
                print("Warning: Markers don't have unique Id.")

    return out

def write_hdf5(fn, states, qid, state=["inistate","endstate"]):
    """
    Write states.

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
    states : dictionary
        State data to be written in dictionary format.
    qid : int
        Run id these states correspond to.
    """

    with h5py.File(fn, "a") as f:
        # Read data from file
        for statename in [state]:

            path = "results/run-" + qid + '/' + statename

            # Remove group if one is already present.
            if path in f:
                del f[path]
            f.create_group(path)

            # TODO Check that inputs are consistent.

            # Read data from file.
            for field in states[state]:
                if field[-4:] != "unit" and field != "N" and field != "uniqueId":
                    d = f.create_dataset(path + "/" + field, data=states[state][field])
                    d.attrs["unit"] = states[state][field + "_unit"]
