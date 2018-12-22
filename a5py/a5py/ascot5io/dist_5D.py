"""
Distribution HDF5 IO module.

File: dist_5D.py
"""
import numpy as np
import h5py
import copy

from a5py.ascot5io.ascot5data import AscotData
import scipy.constants as constants
from scipy.interpolate import griddata, RectBivariateSpline

from a5py.postprocessing.disttrasnformations import vpavpe2Epitch

def read_hdf5(fn, qid):
    """
    Read distributions.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    qid : str
        qid of the run where distribution is read.

    Returns
    -------

    Dictionary storing the distributions that were read.
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
            #out[name + '_unit']  = abscissae_units[i]
            out['n_' + name]     = out[name].size

        out['histogram']      = np.squeeze(dist['ordinate'][:], axis=0)

    return out


def write_hdf5(fn, dists, qid):
    """
    Write distributions.

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
        return read_hdf5(self._file, self.get_qid())

    def histogram2dist(self, dist):
        vol = 1
        for coord in dist["abscissae"]:
            edges = dist[coord + '_edges']
            #dv = (edges[0:-1] + edges[1:]) / 2
            dv = (edges[1:] - edges[0:-1])
            vol = np.multiply.outer(vol, dv)

        dist["density"] = dist["histogram"] / vol
        del dist["histogram"]
        return dist

    def integrate(self, dist, coord, slices=None):
        dim = dist["abscissae"].index(coord)
        edges = dist[coord + '_edges']
        weights = (edges[1:] - edges[0:-1])

        if slices is not None:
            mask = np.zeros(weights.shape)
            mask[slices] = 1
            weights = weights*mask

        dist["density"],s = np.average(dist["density"],
                                       axis=dim, weights=weights,
                                       returned=True)
        dist["density"] *= s
        del dist[coord]
        del dist[coord + '_edges']
        #del dist[coord + '_unit']
        del dist['n_' + coord]
        del dist["abscissae"][dim]
        return dist

    def get_E_xi_dist(self, E_edges, xi_edges, masskg):

        ## Create E-xi distribution ##
        dist = self.histogram2dist(self.read())
        Exidist = copy.deepcopy(dist)

        # Remove vpa and vpe components
        del Exidist["density"]
        Exidist["abscissae"].remove("vpa")
        Exidist["abscissae"].remove("vpe")
        for k in list(Exidist):
            if "vpa" in k or "vpe" in k:
                del Exidist[k]

        # Add E and xi abscissae and initialize a new density
        Exidist["abscissae"].insert(3, "E")
        Exidist["abscissae"].insert(4, "xi")

        Exidist["E"]        = (E_edges[0:-1] + E_edges[1:]) / 2
        Exidist["E_edges"]  = E_edges
        Exidist["n_E"]      = Exidist["E"].size

        Exidist["xi"]       = (xi_edges[0:-1] + xi_edges[1:]) / 2
        Exidist["xi_edges"] = xi_edges
        Exidist["n_xi"]     = Exidist["xi"].size

        Exidist["density"]  = np.zeros((Exidist['n_R'],
                                        Exidist['n_phi'],
                                        Exidist['n_z'],
                                        Exidist["n_E"],
                                        Exidist["n_xi"],
                                        Exidist['n_time'],
                                        Exidist['n_charge']))

        # Transform E-xi grid to points in (vpa,vpa) space that are used in
        # interpolation.
        xig, Eg = np.meshgrid(Exidist["xi"], Exidist["E"])
        vpag = ( xig * np.sqrt( 2*Eg*constants.e/masskg ) ).ravel()
        vpeg = (np.sqrt(1 - xig*xig) * np.sqrt(2*Eg*constants.e/masskg)).ravel()

        # Coordinate transform Jacobian: dvpa dvpe = |jac| dE dxi
        # Jacobian for transform (vpa, vpe) -> (v, xi) is v / sqrt(1-xi^2)
        # because jac = dvpa / dv  = xi, dvpe / dv  = sqrt(1-xi^2)
        #               dvpa / dxi = v,  dvpe / dxi = -xi v / sqrt(1-xi^2),
        # and the Jacobian for (v, xi) -> (E, xi) is e / (mass*v) when
        # E is in electronvolts. Therefore the combined Jacobian is
        # (e/mass) / sqrt(1-xi*xi).
        jac = (constants.e/masskg) / np.sqrt(1 - xig*xig)

        # Interpolate.
        for iR in range(0, dist['n_R']):
            for iphi in range(0, dist['n_phi']):
                for iz in range(0, dist['n_z']):
                    for itime in range(0, dist['n_time']):
                        for iq in range(0, dist['n_charge']):
                            f = RectBivariateSpline(
                                dist['vpa'], dist['vpe'],
                                np.squeeze(dist['density'][iR,iphi,iz,:,:,itime,iq]),
                                kx=1, ky=1)
                            Exidist["density"][iR,iphi,iz,:,:,itime,iq] = \
                                np.reshape(f.ev(vpag, vpeg), (Exidist["n_E"], Exidist["n_xi"])) * jac

        return Exidist

    def get_E_dist(self, E_edges, xi_edges, masskg):
        dist = self.get_E_xi_dist(E_edges, xi_edges, masskg)
        return np.sum(np.squeeze(dist["density"]), 1) \
            * (xi_edges[1] - xi_edges[0])

    def get_xi_dist(self, E_edges, xi_edges, masskg):
        dist = self.get_E_xi_dist(E_edges, xi_edges, masskg)
        return np.sum(np.squeeze(dist["density"]), 0) \
            * (E_edges[1] - E_edges[0])
