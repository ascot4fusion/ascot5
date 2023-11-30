"""Transport coefficient data IO module.
"""
import numpy as np
import h5py

from .coreio.treedata import DataContainer

class Transcoef(DataContainer):
    """Object representing transport coefficient data.
    """

    def __init__(self, root, hdf5):
        """
        Initialize transcoef object from given HDF5 file to given RunNode.
        """
        super().__init__(root, hdf5)

    def read(self):
        """Read raw coefficient data to dictionary.
        """
        out = {}
        with self as f:
            for field in f:
                out[field] = f[field][:]

        return out

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
        """Sort the data into a histogram and return the mean values in each
        bin.

        This function assumes markers were sorted during the creation.
        """
        ids = self["ids"]
        k   = self["k"]
        d   = self["d"]

        ids0 = self._runnode.inistate["id"]
        rho  = self._runnode.inistate["r"]
        if userho:
            rho = self._runnode.inistate["rho"]
        energy = self._runnode.inistate["energy"]
        pitch  = self._runnode.inistate["pitch"]
        ec     = self._runnode.endstate["endcond"]

        lost   = np.logical_or.reduce([ec==32, ec==128, ec==1])

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
            idx  = np.logical_and(  ids > n,  ids <= n + nslot )
            idx0 = np.logical_and( ids0 > n, ids0 <= n + nslot )

            # Evaluate transport from losses instead
            lf = np.sum(lost[ids0[idx0]-1]) / nslot
            if lf >= lossfrac:
                coefs = eval_coefficients(
                    [self._runnode], method=method, endrho=1,
                    ids=[ids0[idx0]-1])
                K[i] = coefs[1]
                D[i] = coefs[2]

                A1[i] = np.mean(rho[ids0[idx0]-1])
                A2[i] = np.mean(energy[ids0[idx0]-1])
                A3[i] = np.mean(pitch[ids0[idx0]-1])

            else:
                K[i] = np.mean(k[idx])
                D[i] = np.mean(d[idx])

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
