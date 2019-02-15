"""
5D distribution module.

File: dist_5D.py
"""
import numpy as np
import h5py

import a5py.dist as distmod
import a5py.marker.interpret as interpret

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid):
    """
    Read 5D distribution from a HDF5 file to a dictionary.

    Args:
        fn : str <br>
            HDF5 file filename.
        qid : str <br>
            QID of the run whose distribution is read.
    Returns:
        Distribution dictionary.
    """

    with h5py.File(fn,"r") as f:

        path = "/results/run-"+qid+"/dists/R_phi_z_vpa_vpe_t_q/"
        dist = f[path]
        out = {}

        # A Short helper function to calculate grid points from grid edges.
        def edges2grid(edges):
            return np.linspace(0.5*(edges[0]+edges[1]),
                               0.5*(edges[-2]+edges[-1]), num=edges.size-1)

        out["abscissae"] = ["R", "phi", "z", "vpa", "vpe", "time", "charge"]
        for i in range(0,len(out["abscissae"])):
            name = out["abscissae"][i]
            out[name + '_edges'] = dist['abscissa_vec_00000'+str(i+1)][:]
            out[name]            = edges2grid(out[name + '_edges'])
            out['n_' + name]     = out[name].size

        out['histogram']      = np.squeeze(dist['ordinate'][:], axis=0)

    return out


class Dist_5D(AscotData):
    """
    Object representing 5D distribution data.
    """

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
            dist = distmod.histogram2density(self.read())
        distmod.squeeze(dist, **kwargs)

        return dist

    def get_E_xi_dist(self, E_edges=None, xi_edges=None, dist=None,
                      **kwargs):
        """
        Return distribution dictionary where (vpa, vpe) is converted to (E, xi).

        This is otherwise same function as Dist_5D.get_dist() except velocity
        components (vpa, vpe) are converted to Energy and pitch (E, xi). Energy
        is in electronvolts and pitch is vpa/(vpa^2 + vpe^2)^0.5.

        Args:
            E_edges : array  (optional) <br>
                Energy grid edges in the new distribution.
            xi_edges : (optional) <br>
                Pitch grid edges in the new distribution.
            dist :(optional) <br>
                Use this distribution instead of reading one from HDF5 file.
            kwargs : <br>
                Name(s) (R, phi, z, vpa, vpe, time, charge) of those dimensions
                along which the distribution is either sliced or integrated
                over. Names are given as keyword value pairs where a tuple value
                (a, b) are indices for slicing and array [a, b]are indices for
                integration. A scalar zero means whole dimension is integrated.

        Returns:
            Distribution dictionary.
        """

        if not dist:
            dist = distmod.histogram2density(self.read())

        # Mass from inistate.
        masskg = interpret.mass_kg(self._runnode.inistate["mass"][0])

        Exidist = distmod.convert_vpavpe_to_Exi(dist, masskg, E_edges, xi_edges)
        distmod.squeeze(Exidist, **kwargs)

        return Exidist

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
        abscissae = {"R" : 0, "phi" : 0, "z" : 0, "vpa" : 0,
                     "vpe" : 0, "time" : 0, "charge" : 0}

        x = args[0]
        del abscissae[x]
        y = None
        if len(args) > 1:
            y = args[1]
            del abscissae[y]

        if not dist:
            dist = self.get_dist()

        for k in abscissae.keys():
            if k not in dist["abscissae"]:
                del abscissae[k]

        distmod.squeeze(dist, **abscissae)

        if not y:
            distmod.plot_dist_1D(dist, logscale=logscale, axes=axes)
        else:
            distmod.plot_dist_2D(dist, x, y, logscale=logscale, equal=equal,
                                 axes=axes)

    def plot_E_xi_dist(self, *args, E_edges=None, xi_edges=None,
                       logscale=False, equal=False, axes=None, dist=None):
        """
        Convert (vpa, vpe) to (E, xi) and plot the distribution.

        Either makes a 1D or a 2D plot depending on number of input arguments.

        Args:
            args : str, str <br>
                Name of the x-coordinate, and optionally y-coordinate if a 2D
                plot is desired.
            axes : Axes, optional <br>
                Axes where plotting happens. If None, a new figure will be
                created.
            equal : bool, optional <br>
                Make axes equal.
            dist : dict_like, optional <br>
               Give input distribution explicitly instead of reading one from
               HDF5 file. Dimensions that are not x or y are integrated over.
        """
        abscissae = {"R" : 0, "phi" : 0, "z" : 0, "E" : 0,
                     "xi" : 0, "time" : 0, "charge" : 0}

        x = args[0]
        del abscissae[x]
        y = None
        if len(args) > 1:
            y = args[1]
            del abscissae[y]

        if not dist:
            dist = self.get_E_xi_dist(E_edges=E_edges, xi_edges=xi_edges)

        for k in list(abscissae.keys()):
            if k not in dist["abscissae"]:
                del abscissae[k]

        distmod.squeeze(dist, **abscissae)

        if not y:
            distmod.plot_dist_1D(dist, logscale=logscale, axes=axes)
        else:
            distmod.plot_dist_2D(dist, x, y, logscale=logscale, equal=equal, axes=axes)
