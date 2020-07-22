"""
Transport coefficient data IO module.

File: transcoef.py
"""
import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid):
    """
    Read transport coefficient data from HDF5 file to dictionary.

    Args:
        fn : str <br>
            Full path to HDF5 file.
        qid : str <br>
            QID of the run whose transport coefficient data is read.

    Returns:
        Dictionary storing the coefficients and ids (unsorted).
    """

    with h5py.File(fn,"r") as f:
        transcoef = f["/results/run_"+qid+"/transcoef"]

        out = {}

        # Read data from file.
        for field in transcoef:
            out[field]           = data[field][:]

    return out


class Transcoef(AscotData):
    """
    Object representing transport coefficient data.
    """

    def __init__(self, hdf5, runnode):
        """
        Initialize transcoef object from given HDF5 file to given RunNode.
        """
        self._runnode = runnode
        super().__init__(hdf5)


    def read(self):
        """
        Read coefficient data to dictionary.
        """
        return read_hdf5(self._file, self.get_qid())


    def __getitem__(self, key):
        with self as h5:
            h5keys = list(h5.keys())
            item = None
            if key in h5keys:
                item = h5[key][:]

            return item
