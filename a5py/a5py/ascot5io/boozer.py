"""
Boozer coordinate input IO.

File: boozer.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

## ADD COMMENTS AND SIZE DEFINITIONS OF THE DATA FIELDS

def write_hdf5(fn, psimin, psimax, npsi, thetamin, thetamax, ntheta, rmin,
               rmax, nr, zmin, zmax, nz, r0, z0, psiin, psiout, psi_rz,
               theta_psithetageom, nu_psitheta, rs, zs, nrzs, desc=None):
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

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        # grid specifications
        g.create_dataset("psi_min",    (1,), data=psimin,     dtype="f8")
        g.create_dataset("psi_max",    (1,), data=psimax,     dtype="f8")
        g.create_dataset("npsi",       (1,), data=npsi,       dtype="i8")
        g.create_dataset("thetamin",   (1,), data=thetamin,   dtype="f8")
        g.create_dataset("thetamax",   (1,), data=thetamax,   dtype="f8")
        g.create_dataset("ntheta",     (1,), data=ntheta,     dtype="i8")
        g.create_dataset("r_min",      (1,), data=rmin,       dtype="f8")
        g.create_dataset("r_max",      (1,), data=rmax,       dtype="f8")
        g.create_dataset("nr",         (1,), data=nr,         dtype="i8")
        g.create_dataset("z_min",      (1,), data=zmin,       dtype="f8")
        g.create_dataset("z_max",      (1,), data=zmax,       dtype="f8")
        g.create_dataset("nz",         (1,), data=nz,         dtype="i8")
        g.create_dataset("r0",         (1,), data=r0,         dtype="f8")
        g.create_dataset("z0",         (1,), data=z0,         dtype="f8")
        g.create_dataset("nrzs",       (1,), data=nrzs,       dtype="i8")

        # the outermost poloidal psi-surface contour
        g.create_dataset("rs", (nrzs,), data=rs, dtype="f8")
        g.create_dataset("zs", (nrzs,), data=zs, dtype="f8")

        # psi data min and max values
        g.create_dataset("psi_inner", (1,), data=psiin,  dtype="f8")
        g.create_dataset("psi_outer", (1,), data=psiout, dtype="f8")

        # tabulated coordinates maps
        g.create_dataset("psi_rz", (nr,nz), data=psi_rz, dtype="f8")
        g.create_dataset("theta_psithetageom", (npsi,ntheta), data=theta_psithetageom, dtype="f8")
        g.create_dataset("nu_psitheta", (npsi,ntheta), data=nu_psitheta, dtype="f8")


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

def write_hdf5_dummy(fn, desc="Dummy"):
    """
    Write dummy boozer input.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """

    psimin     = 0
    psimax     = 1
    npsi       = 6
    thetamin   = 0.0
    thetamax   = 2*np.math.pi
    ntheta     = 10
    rmin       = 0.1
    rmax       = 10.0
    nr         = 5
    zmin       = -10.0
    zmax       = 10.0
    nz         = 10
    r0         = (rmax+rmin)/2.0
    z0         = (zmin+zmax)/2.0
    psiin      = 0
    psiout     = 1
    nrzs       = ntheta

    rs = np.cos(np.linspace(0,2*np.math.pi,ntheta))
    zs = np.sin(np.linspace(0,2*np.math.pi,ntheta))

    psi_rz    = np.ones((nr,nz))
    theta_psithetageom = np.ones((npsi,ntheta))
    nu_psitheta = np.ones((npsi,ntheta))

    write_hdf5(fn, psimin, psimax, npsi, thetamin, thetamax, ntheta, rmin,
               rmax, nr, zmin, zmax, nz, r0, z0, psiin, psiout, psi_rz,
               theta_psithetageom, nu_psitheta, rs, zs, nrzs, desc=desc)
