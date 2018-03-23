"""
Distribution HDF5 IO module.
"""
import numpy as np
import h5py
import random
import datetime

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

    f = h5py.File(fn,"r")

    path = "/results/run-"+qid+"/dists"
    dists = f[path]
    out = {}

    # A Short helper function to calculate grid points from grid edges.
    def edges2grid(edges):
        return np.linspace(0.5*(edges[0]+edges[1]), 0.5*(edges[-2]+edges[-1]), num=edges.size-1)

    if "R_phi_z_vpa_vpe_t_q" in dists:
        out["R_phi_z_vpa_vpe_t_q"] = {}
        temp = out["R_phi_z_vpa_vpe_t_q"]
        disttemp = dists["R_phi_z_vpa_vpe_t_q"]

        # These could be read directly from HDF5 file, but for clarity we list them here
        abscissae_names = ["R", "phi", "z", "vpa", "vpe", "time", "charge"]
        abscissae_units = ["m", "deg", "m", "m/s", "m/s", "s", "e"]
        abscissae_realnames = ["Major radius", "Toroidal angle", "Height", "Velocity parallel to magnetic field", "Velocity perpendicular to magnetic field", "Time", "Charge"]

        for i in range(0,len(abscissae_names)):
            name = abscissae_names[i]
            temp[name + '_edges'] = disttemp['abscissa_vec_00000'+str(i+1)][:]
            temp[name]            = edges2grid(temp[name + '_edges'])
            temp[name + '_unit']  = abscissae_units[i]
            temp['n_' + name]     = temp[name].size

        temp['ordinate']       = disttemp['ordinate'][0,:,:,:,:,:,:,:]
        temp['ordinate_name'] = 'density'
        temp['ordinate_unit'] = 's/m^5*deg*e'

    f.close()

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

    f = h5py.File(fn, "a")

    path = "results/run-" + qid + "/dists/"

    # Remove group if one is already present.
    if path in f:
        del f[path]
        f.create_group(path)

    # TODO Check that inputs are consistent.

    # Write data to file.
    if "rzVDist" in dists:
        f.create_group(path + "rzVDist")
        f.create_group(path + "rzVDist/abscissae")

        f.create_dataset(path + "rzVDist/abscissae/dim1", data=dists["rzVDist"]["R_edges"])
        f.create_dataset(path + "rzVDist/abscissae/dim2", data=dists["rzVDist"]["z_edges"])
        f.create_dataset(path + "rzVDist/abscissae/dim3", data=dists["rzVDist"]["vpa_edges"])
        f.create_dataset(path + "rzVDist/abscissae/dim4", data=dists["rzVDist"]["vpe_edges"])
        f.create_dataset(path + "rzVDist/abscissae/dim5", data=dists["rzVDist"]["time_edges"])

        f.create_dataset(path + "rzVDist/ordinate", data=dists["rzVDist"]["ordinate"])

    if "R_phi_z_vpa_vpe_t_q" in dists:
        f.create_group(path + "R_phi_z_vpa_vpe_t_q")

        abscissae_names = ["R", "phi", "z", "vpa", "vpe", "time", "charge"]
        f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/" + "abscissa_ndim", data=len(abscissae_names))
        for i, name in enumerate(abscissae_names):
            f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/" + "abscissa_name_00000" + str(i+1),
                             data=name)
            f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/" + "abscissa_vec_00000" + str(i+1),
                             data=dists["R_phi_z_vpa_vpe_t_q"][name + '_edges'])
            f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/" + "abscissa_unit_00000" + str(i+1),
                             data=dists["R_phi_z_vpa_vpe_t_q"][name + '_unit'])
            f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/" + "abscissa_nslot_00000" + str(i+1),
                             data=dists["R_phi_z_vpa_vpe_t_q"]['n_' + name])

        f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/ordinate",
                         data=np.expand_dims(dists["R_phi_z_vpa_vpe_t_q"]["ordinate"],0))
        f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/ordinate_name_000001",
                         data=dists["R_phi_z_vpa_vpe_t_q"]["ordinate_name"])
        f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/ordinate_ndim",data=1)
        f.create_dataset(path + "R_phi_z_vpa_vpe_t_q/ordinate_unit_000001",
                         data=dists["R_phi_z_vpa_vpe_t_q"]["ordinate_unit"])

    f.close()
