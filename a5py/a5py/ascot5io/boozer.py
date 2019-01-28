"""
Boozer coordinate input IO.

File: boozer.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, gprof, qprof, Iprof, delta, nu, theta_bzr, theta_geo,
               desc=None):
    """
    Write boozer input to HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        desc : str, optional <br>
            Input's description.
    """

    parent = "boozer"
    group  = "Boozer"

    npsi       = qprof.shape[0]
    nR         = theta_geo.shape[0]
    nz         = theta_geo.shape[1]
    ntheta_bzr = delta.shape[1]
    ntheta_geo = theta_bzr.shape[1]

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("nR",         (1,), data=nR,         dtype="i8")
        g.create_dataset("nz",         (1,), data=nz,         dtype="i8")
        g.create_dataset("npsi",       (1,), data=npsi,       dtype="i8")
        g.create_dataset("ntheta_geo", (1,), data=ntheta_geo, dtype="i8")
        g.create_dataset("ntheta_bzr", (1,), data=ntheta_bzr, dtype="i8")

        g.create_dataset("g",         (npsi,),           data=gprof, dtype="f8")
        g.create_dataset("q",         (npsi,),           data=qprof, dtype="f8")
        g.create_dataset("I",         (npsi,),           data=Iprof, dtype="f8")
        g.create_dataset("delta",     (npsi,ntheta_bzr), data=delta, dtype="f8")
        g.create_dataset("nu",        (npsi,ntheta_bzr), data=nu,    dtype="f8")

        g.create_dataset("theta_bzr", (npsi,ntheta_geo), data=theta_bzr,
                         dtype="f8")
        g.create_dataset("theta_geo", (nR,nz),           data=theta_geo,
                         dtype="f8")


def read_hdf5(fn, qid):
    """
    Read Boozer input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            qid of the efield to be read.

    Returns:
        Dictionary containing Boozer data.
    """

    path = "boozer" + "/Boozer-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["desc"] = f[path].attrs["desc"]

        # Actual data.
        for k in f[path].keys():
            out[k] = f[path][k][:]

    return out

def write_hdf5_dummy(fn):
    """
    Write dummy boozer input.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """
    npsi       = 6
    nR         = 5
    nz         = 10
    ntheta_bzr = 8
    ntheta_geo = 12

    gprof     = np.ones((npsi,))
    qprof     = np.ones((npsi,))
    Iprof     = np.ones((npsi,))
    delta     = np.ones((npsi,ntheta_bzr))
    nu        = np.ones((npsi,ntheta_bzr))
    theta_bzr = np.ones((npsi,ntheta_geo))
    theta_geo = np.ones((nR,nz))
    write_hdf5(fn, gprof, qprof, Iprof, delta, nu, theta_bzr, theta_geo,
               desc="Dummy")
