"""
Options IO.

File: options.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, options, desc=None):
    """
    Write options.

    Unlike most other "write" functions, this one takes dictionary
    as an argument. The dictionary should have exactly the same format
    as given by the "read" function in this module.

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    options : dictionary
        Options to be written in dictionary format.
    """

    parent = "options"
    group  = "opt"

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        for opt in options:
            if opt != "qid" and opt != "date" and opt != "description":
                g.create_dataset(opt, (options[opt].size,), data=options[opt])


def read_hdf5(fn, qid):
    """
    Read options from HDF5 file.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file.
    qid : str
        qid of the options to be read.

    Returns
    -------

    Dictionary containing options.
    """

    path = "options" + "/opt-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for opt in f[path]:
            out[opt] = f[path][opt][:]

    return out

class Opt(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
