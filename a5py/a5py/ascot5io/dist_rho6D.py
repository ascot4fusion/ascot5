"""
Distribution HDF5 IO module.

File: dist_rho6D.py
"""
import numpy as np
import h5py

from a5py.ascot5io.ascot5data import AscotOutput

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

        path = "/results/run-"+qid+"/dists/rho_pol_phi_vr_vphi_vz_t_q"
        dist = f[path]

        # A Short helper function to calculate grid points from grid edges.
        def edges2grid(edges):
            return np.linspace(0.5*(edges[0]+edges[1]),
                               0.5*(edges[-2]+edges[-1]), num=edges.size-1)

        out = {}

        # These could be read directly from HDF5 file, but for clarity
        # we list them here
        abscissae_names = ["rho", "pol", "phi", "vr", "vphi", "vz", "time", "charge"]
        abscissae_units = ["", "deg", "deg", "m/s", "m/s", "m/s", "s", "e"]
        abscissae_realnames = ["Radial coordinate", "Poloidal angle",
                               "Toroidal angle",
                               "Velocity R component", "Velocity phi component",
                               "Velocity z component",
                               "Time", "Charge"]

        for i in range(0,len(abscissae_names)):
            name = abscissae_names[i]
            out[name + '_edges'] = dist['abscissa_vec_00000'+str(i+1)][:]
            out[name]            = edges2grid(out[name + '_edges'])
            out[name + '_unit']  = abscissae_units[i]
            out['n_' + name]     = out[name].size

            out['ordinate']       = dist['ordinate'][0,:,:,:,:,:,:,:]
            out['ordinate_name'] = 'density'
            out['ordinate_unit'] = 's/m^2*deg^2*e'

    return out


def write_hdf5(fn, dist, qid):
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

        path = "results/run-" + qid + "/dists/rho_pol_phi_vpa_vpe_t_q/"

        # Remove group if one is already present.
        if path in f:
            del f[path]
            f.create_group(path)

        # TODO Check that inputs are consistent.

        # Write data to file.
        f.create_group(path + "rho_pol_phi_vr_vphi_vz_t_q")

        abscissae_names = ["rho", "pol", "phi", "vr", "vphi", "vz", "time", "charge"]
        f.create_dataset(path + "abscissa_ndim", data=len(abscissae_names))
        for i, name in enumerate(abscissae_names):
            f.create_dataset(path + "abscissa_name_00000" + str(i+1),
                             data=name)
            f.create_dataset(path + "abscissa_vec_00000" + str(i+1),
                             data=dists["rho_pol_phi_vpa_vpe_t_q"][name + '_edges'])
            f.create_dataset(path + "abscissa_unit_00000" + str(i+1),
                             data=dist[name + '_unit'])
            f.create_dataset(path + "abscissa_nslot_00000" + str(i+1),
                             data=dist['n_' + name])

        f.create_dataset(path + "ordinate",
                         data=np.expand_dims(dist["ordinate"],0))
        f.create_dataset(path + "ordinate_name_000001",
                         data=dist["ordinate_name"])
        f.create_dataset(path + "ordinate_ndim",data=1)
        f.create_dataset(path + "ordinate_unit_000001",
                         data=dist["ordinate_unit"])

class Dist_rho6D(AscotOutput):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
