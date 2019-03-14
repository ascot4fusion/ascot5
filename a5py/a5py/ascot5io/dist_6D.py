"""
Distribution HDF5 IO module.

File: dist_6D.py
"""
import numpy as np
import h5py

import a5py.dist as distmod
from a5py.marker.alias import get as alias

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid):
    """
    Read 6D distribution.

    Args:
        fn : str <br>
            Full path to HDF5 file.
        qid : str <br>
            QID of the run from where the distribution is read.
    Returns:
        Dictionary storing the distributions that were read.
    """

    with h5py.File(fn,"r") as f:

        path = "/results/run_"+qid+"/dist6d/"
        dist = f[path]
        out = {}

        # A Short helper function to calculate grid points from grid edges.
        def edges2grid(edges):
            return np.linspace(0.5*(edges[0]+edges[1]),
                               0.5*(edges[-2]+edges[-1]), num=edges.size-1)
        # Abscissa info
        abscissae = ["r", "phi", "z", "vr", "vphi", "vz", "time", "charge"]
        abscissae_units = ["m", "deg", "m", "m/s", "m/s", "m/s", "s", "e"]

        for i in range(0,len(abscissae)):
            name = abscissae[i]
            out[name + "_edges"] = dist["abscissa_vec_0"+str(i+1)][:]
            out[name]            = edges2grid(out[name + "_edges"])
            out[name + "_unit"]  = abscissae_units[i]
            out["n" + name]      = out[name].size

        out["abscissae"] = abscissae
        out["histogram"] = dist["ordinate"][0,:,:,:,:,:,:,:,:]

    return out


class Dist_6D(AscotData):

    def __init__(self, hdf5, runnode):
        """
        Object representing orbit data.
        """
        self._runnode = runnode
        super().__init__(hdf5)


    def read(self):
        """
        Read distribution data from HDF5 file to a dictionary.

        Returns:
            Distribution dictionary.
        """
        return read_hdf5(self._file, self.get_qid())


    def get_dist(self, dist=None, **kwargs):
        """
        Return distribution dictionary.

        The ordinate is density which is integrated over the dimensions which
        are given in kwargs.

        Args:
            dist : dict_like, optional <br>
               Give input distribution explicitly instead of reading one from
               HDF5 file.
            kwargs : Name(s) (R, phi, z, vpa, vpe, time, charge) of those
               dimensions along which the distribution is either sliced or
               integrated over. Names are given as keyword value pairs where
               a tuple value (a, b) are indices for slicing and array [a, b]
               are indices for integration. A scalar zero means whole
               dimension is integrated.

        Returns:
            Distribution dictionary.
        """
        if not dist:
            dist = distmod.histogram2distribution(self.read())
        distmod.squeeze(dist, **kwargs)

        return dist


    def plot_dist(self, *args, logscale=False, equal=False, axes=None,
                  dist=None):
        """
        Plot distribution.

        Either makes a 1D or a 2D plot depending on number of input arguments.

        Args:
            args : str, str <br>
                Name of the x-coordinate, and optionally y-coordinate if a 2D
                plot is desired.
            equal : bool, optional <br>
                Make axes equal.
            axes : Axes, optional <br>
                Axes where plotting happens. If None, a new figure will be
                created.
            dist : dict_like, optional <br>
               Give input distribution explicitly instead of reading one from
               HDF5 file. Dimensions that are not x or y are integrated over.
        """
        abscissae = {"r" : 0, "phi" : 0, "z" : 0, "vr" : 0,
                     "vphi" : 0, "vz" : 0, "time" : 0, "charge" : 0}

        x = alias(args[0])
        del abscissae[x]
        y = None
        if len(args) > 1:
            y = alias(args[1])
            del abscissae[y]

        if not dist:
            dist = self.get_dist()

        for k in list(abscissae.keys()):
            if k not in dist["abscissae"]:
                del abscissae[k]

        distmod.squeeze(dist, **abscissae)

        if not y:
            distmod.plot_dist_1D(dist, logscale=logscale, axes=axes)
        else:
            distmod.plot_dist_2D(dist, x, y, logscale=logscale, equal=equal,
                                 axes=axes)
