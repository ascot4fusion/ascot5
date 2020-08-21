"""
Transport coefficient data IO module.

File: transcoef.py
"""
import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotData
from a5py.marker.losstime import eval_coefficients

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
        data = f["/results/run_"+qid+"/transcoef"]

        out = {}

        # Read data from file.
        for field in data:
            out[field] = data[field][:]

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

            # Sort by id
            idx  = h5["ids"][:].argsort()
            return item[idx]


    def get_matrix(self, dim=(1,), lossfrac=1.1, method="radindep",
                   userho=0):
        """
        Sort the data into a histogram and return the mean values in each bin.

        This function assumes markers were sorted during the creation.
        """
        ids = self["ids"]
        k   = self["k"]
        d   = self["d"]

        rho = self._runnode.inistate["r"]
        if userho:
            rho = self._runnode.inistate["rho"]
        energy = self._runnode.inistate["energy"]
        pitch  = self._runnode.inistate["pitch"]
        ec     = self._runnode.inistate["endcond"]

        nsize = 1
        for i in np.arange( len(dim) ):
            nsize *= dim[i]

        ntot = self._runnode.inistate["ids"].size
        if ntot % nsize > 0:
            print("The number of markers must be a multiple of dimensions.")
            return

        nslot = int(ntot / nsize)
        K = np.zeros(dim).ravel()
        D = np.zeros(dim).ravel()
        A1 = np.zeros(dim).ravel()
        A2 = np.zeros(dim).ravel()
        A3 = np.zeros(dim).ravel()

        i = 0
        n = 0
        while n < ntot:
            idx = np.logical_and( ids > n, ids <= n + nslot )
            K[i] = np.mean(k[idx])
            D[i] = np.mean(d[idx])

            # Evaluate transport from losses instead
            lf = np.sum(np.logical_or.reduce([ec==32, ec==128, ec==1])[ids[idx]-1]) / nslot
            if lf >= lossfrac:
                coefs = eval_coefficients(
                    [self._runnode], method=method, endrho=1,
                    ids=[ids[idx]-1])
                K[i] = coefs[1]
                D[i] = coefs[2]

            A1[i] = np.mean(rho[ids[idx]-1])
            A2[i] = np.mean(energy[ids[idx]-1])
            A3[i] = np.mean(pitch[ids[idx]-1])

            n = n + nslot
            i = i + 1

        K = np.reshape(K, dim)
        D = np.reshape(D, dim)
        A1 = np.reshape(A1, dim)
        A2 = np.reshape(A2, dim)
        A3 = np.reshape(A3, dim)
        return (A1, A2, A3, K, D)
