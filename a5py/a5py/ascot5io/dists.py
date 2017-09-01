"""
Distribution HDF5 IO module.
"""
import numpy as np
import h5py
import random
import datetime

def read_hdf5(fn):
    """
    Read distributions.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.

    Returns
    -------

    Dictionary storing the distributions that were read.
    """

    f = h5py.File(fn,"r")
    dists = f["distributions"]
    out = {}

    if "rzVDist" in dists:
        out["rzVDist"]    = {}
        temp = out
        out  = temp["rzVDist"]

        out['R_edges']    = dists['rzVDist/abscissae/dim1'][:]
        out['R']          = np.linspace(0.5*(out['R_edges'][0]+out['R_edges'][1]), 0.5*(out['R_edges'][-2]+out['R_edges'][-1]), num=out['R_edges'].size-1)
        out['R_unit']     = 'm'
        out['nR']         = out['R'].size

        out['z_edges']    = dists['rzVDist/abscissae/dim2'][:]
        out['z']          = np.linspace(0.5*(out['z_edges'][0]+out['z_edges'][1]), 0.5*(out['z_edges'][-2]+out['z_edges'][-1]), num=out['z_edges'].size-1)
        out['z_unit']     = 'm'
        out['nz']         = out['z'].size

        out['vpa_edges']  = dists['rzVDist/abscissae/dim3'][:]
        out['vpa']        = np.linspace(0.5*(out['vpa_edges'][0]+out['vpa_edges'][1]), 0.5*(out['vpa_edges'][-2]+out['vpa_edges'][-1]), num=out['vpa_edges'].size-1)
        out['vpa_unit']   = 'm/s'
        out['nvpa']       = out['vpa'].size

        out['vpe_edges']  = dists['rzVDist/abscissae/dim4'][:]
        out['vpe']        = np.linspace(0.5*(out['vpe_edges'][0]+out['vpe_edges'][1]), 0.5*(out['vpe_edges'][-2]+out['vpe_edges'][-1]), num=out['vpe_edges'].size-1)
        out['vpe_unit']   = 'm/s'
        out['nvpe']       = out['vpe'].size

        out['time_edges'] = dists['rzVDist/abscissae/dim5'][:]
        out['time']       = np.linspace(0.5*(out['time_edges'][0]+out['time_edges'][1]), 0.5*(out['time_edges'][-2]+out['time_edges'][-1]), num=out['time_edges'].size-1)
        out['time_unit']  = 's'
        out['ntime']      = out['time'].size

        # In HDF5 the ordinate shape is, for some reason, (1,1,vpe,vpa,z,R,1), one of those 1's is time but we ignore that
        out['ordinate']       = np.transpose(dists['rzVDist/ordinate'][0,0,:],(3,2,1,0,4))
        out['ordinate_names'] = ['R','z','vpa','vpe','time']

        temp["rzVDist"] = out
        out = temp

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

    path = "distributions/"

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

    f.close()
