"""
State HDF5 IO module.

File: state.py
"""

import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotData
import a5py.postprocessing.markereval as meval

def read_hdf5(fn, qid, name):
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

    path = "results/run-" + qid + "/" + name

    with h5py.File(fn,"r") as f:
        out = {}

        for field in f[path]:
            out[field]           = f[path][field][:]
            out[field + "_unit"] = f[path][field].attrs["unit"]

        # TODO Parse endconditions.

        # Find number of markers and check that no markers share same id
        # (which they shouldn't).
        out["N"] = np.unique(out["id"]).size
        if out["N"] != out["id"].size:
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

class State(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid(), self.get_type())

    def __getitem__(self, key):
        mode = "gc"
        h5   = self._open()
        item = meval.evaluate(h5, key, mode)
        self._close()
        return item
