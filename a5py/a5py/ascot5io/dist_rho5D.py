"""
Distribution HDF5 IO module.

File: dist_rho5D.py
"""
import numpy as np
import h5py

import a5py.dist as distmod
from a5py.marker.alias import get as alias
import a5py.marker.interpret as interpret

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid):
    """
    Read rho 5D distribution.

    Args:
        fn : str <br>
            Full path to HDF5 file.
        qid : str <br>
            QID of the run where distribution is read.

    Returns:
        Dictionary storing the distributions that were read.
    """

    with h5py.File(fn,"r") as f:

        path = "/results/run_"+qid+"/distrho5d/"
        dist = f[path]
        out = {}

        # A Short helper function to calculate grid points from grid edges.
        def edges2grid(edges):
            return np.linspace(0.5*(edges[0]+edges[1]),
                               0.5*(edges[-2]+edges[-1]), num=edges.size-1)

        # These could be read directly from HDF5 file, but for clarity
        # we list them here
        abscissae = ["rho", "theta", "phi", "vpar", "vperp", "time", "charge"]
        abscissae_units = ["", "deg", "deg", "m/s", "m/s", "s", "e"]

        for i in range(0,len(abscissae)):
            name = abscissae[i]
            out[name + "_edges"] = dist["abscissa_vec_0"+str(i+1)][:]
            out[name]            = edges2grid(out[name + "_edges"])
            out[name + "_unit"]  = abscissae_units[i]
            out["n" + name]      = out[name].size

        out["abscissae"] = abscissae
        out["histogram"] = dist["ordinate"][0,:,:,:,:,:,:,:]

    return out


class Dist_rho5D(AscotData):
    """
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
            dist = distmod.histogram2distribution(self.read())
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
            dist = distmod.histogram2distribution(self.read())

        # Mass from inistate.
        masskg = interpret.mass_kg(self._runnode.inistate["mass"][0])

        Exidist = distmod.convert_vpavpe_to_Exi(dist, masskg, E_edges, xi_edges)
        distmod.squeeze(Exidist, **kwargs)

        return Exidist

    def plot_dist(self, *args, logscale=False, equal=False, axes=None,
                  dist=None,label=None):
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
        abscissae = {"rho" : 0, "theta" : 0, "phi" : 0, "vpar" : 0,
                     "vperp" : 0, "time" : 0, "charge" : 0}

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
            distmod.plot_dist_1D(dist, logscale=logscale, axes=axes,label=label)
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
        abscissae = {"rho" : 0, "theta" : 0, "phi" : 0, "energy" : 0,
                     "pitch" : 0, "time" : 0, "charge" : 0}

        x = alias(args[0])
        del abscissae[x]
        y = None
        if len(args) > 1:
            y = alias(args[1])
            del abscissae[y]

        if not dist:
            dist = self.get_E_xi_dist(E_edges=E_edges, xi_edges=xi_edges)

        for k in abscissae.keys():
            if k not in dist["abscissae"]:
                del abscissae[k]

        distmod.squeeze(dist, **abscissae)

        if not y:
            distmod.plot_dist_1D(dist, logscale=logscale, axes=axes)
        else:
            distmod.plot_dist_2D(dist, x, y, logscale=logscale, equal=equal,
                                 axes=axes)



    def eval_1d_dist(self, quantity,rmin, rmax, zmin, zmax):
        '''
            The R/z min/max are needed for numerical calculation of the rho slot volumes.
            They should cover the whole plasma.
        '''
        from a5py.ascotpy import Ascotpy
        a5 = Ascotpy(self._file)
        a5.init(bfield=self._runnode.bfield.get_qid(),
                plasma=self._runnode.plasma.get_qid())

        ma = self._runnode.inistate["mass"][0]
        print('mass {} kg'.format(ma)) 
        # ma *= const.physical_constants["atomic mass constant"][0]
        qa = self._runnode.inistate["charge"][0]
        # qa *=  * const.e
        print('charge {} C'.format(qa)) 

        

        dist = distmod.eval1d_fromRhodist(a5, self.get_dist(), quantity, 
                                          rmin, rmax, zmin, zmax, ma=ma, qa=qa)


        a5.free(bfield=True, plasma=True)
        
        return dist
        