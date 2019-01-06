"""
5D distribution module.

File: dist_5D.py
"""
import numpy as np
import h5py

import a5py.dist as distmod

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


def write_hdf5(fn, dists, qid):
    """
    Write 5D distribution to a HDF5 file.

    TODO not compatible with new format

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    orbits : dictionary
        Distributions data to be written in dictionary format.
    qid : int
        Run id these distributions correspond to.
    """

    with h5py.File(fn, "a") as f:

        path = "results/run-" + qid + "/dists/R_phi_z_vpa_vpe_t_q/"

        # Remove group if one is already present.
        if path in f:
            del f[path]
            f.create_group(path)

        # TODO Check that inputs are consistent.

        # Write data to file.

        f.create_group(path)

        abscissae_names = ["R", "phi", "z", "vpa", "vpe", "time", "charge"]
        f.create_dataset(path + "abscissa_ndim", data=len(abscissae_names))
        for i, name in enumerate(abscissae_names):
            f.create_dataset(path + "abscissa_name_00000" + str(i+1),
                             data=name)
            f.create_dataset(path +  "abscissa_vec_00000" + str(i+1),
                             data=dist[name + '_edges'])
            f.create_dataset(path + "abscissa_unit_00000" + str(i+1),
                             data=dist[name + '_unit'])
            f.create_dataset(path + "abscissa_nslot_00000" + str(i+1),
                             data=dist['n_' + name])

        f.create_dataset(path + "ordinate",
                         data=np.expand_dims(dists["ordinate"],0))
        f.create_dataset(path + "ordinate_name_000001",
                         data=dist["ordinate_name"])
        f.create_dataset(path + "ordinate_ndim",data=1)
        f.create_dataset(path + "ordinate_unit_000001",
                         data=dists["ordinate_unit"])

class Dist_5D(AscotData):

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

    def get_E_xi_dist(self, masskg, E_edges=None, xi_edges=None, dist=None,
                      **kwargs):
        """
        Return distribution dictionary where (vpa, vpe) is converted to (E, xi).

        This is otherwise same function as Dist_5D.get_dist() except velocity
        components (vpa, vpe) are converted to Energy and pitch (E, xi). Energy
        is in electronvolts and pitch is vpa/(vpa^2 + vpe^2)^0.5.

        Args:
            masskg : float <br>
                Mass of the species (required for energy conversion) in kg. Note
                that distribution is assumed to consist of markers with equal
                mass.
            E_edges : array  (optional) Energy grid edges in the new distribution.
            xi_edges (optional) Pitch grid edges in the new distribution.
            dist     (optional) Use this distribution instead of reading one
                     from HDF5 file.
            kwargs   Name(s) (R, phi, z, vpa, vpe, time, charge) of those
                     dimensions along which the distribution is either sliced or
                     integrated over. Names are given as keyword value pairs
                     where a tuple value (a, b) are indices for slicing and
                     array [a, b]are indices for integration. A scalar zero
                     means whole dimension is integrated.

        Returns:
            Distribution dictionary.
        """

        if not dist:
            dist = distmod.histogram2density(self.read())

        Exidist = distmod.convert_vpavpe_to_Exi(dist, masskg, E_edges, xi_edges)
        distmod.squeeze(Exidist, **kwargs)

        return Exidist

    def plot_dist(self, *args, axes=None, dist=None):
        """
        Plot distribution.

        Either makes a 1D or a 2D plot depending on number of input arguments.

        Args:
            args : str, str <br>
                Name of the x-coordinate, and optionally y-coordinate if a 2D
                plot is desired.
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
            distmod.plot_dist_1D(dist, axes=axes)
        else:
            distmod.plot_dist_2D(dist, x, y, axes=axes)

    def plot_E_xi_dist(self, masskg, *args, E_edges=None, xi_edges=None,
                       axes=None, dist=None):
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
            dist = self.get_E_xi_dist(masskg,
                                      E_edges=E_edges, xi_edges=xi_edges)

        for k in abscissae.keys():
            if k not in dist["abscissae"]:
                del abscissae[k]

        distmod.squeeze(dist, **abscissae)

        if not y:
            distmod.plot_dist_1D(dist, axes=axes)
        else:
            distmod.plot_dist_2D(dist, x, y, axes=axes)
