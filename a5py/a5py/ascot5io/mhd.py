"""
MHD input IO.

File: mhd.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, nmodes, mmodes, amplitude, omega, alpha, phi,
               desc=None):
    """
    Write MHD input to HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        desc : str, optional <br>
            Input's description.
    """

    parent = "mhd"
    group  = "MHD"

    n_modes = alpha.shape[0]
    npsi    = alpha.shape[1]

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)

        g.create_dataset("n_modes", (1,), data=n_modes, dtype="i8")
        g.create_dataset("npsi",    (1,), data=npsi,    dtype="i8")

        g.create_dataset("nmode", (n_modes,), data=nmodes, dtype="i8")
        g.create_dataset("mmode", (n_modes,), data=mmodes, dtype="i8")

        g.create_dataset("amplitude_nm", (n_modes,), data=amplitude, dtype="f8")
        g.create_dataset("omega_nm",     (n_modes,), data=omega, dtype="f8")

        g.create_dataset("alpha_nm", (n_modes,npsi), data=alpha, dtype="f8")
        g.create_dataset("phi_nm",   (n_modes,npsi), data=phi, dtype="f8")


def read_hdf5(fn, qid):
    """
    Read MHD input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            qid of the efield to be read.

    Returns:
        Dictionary containing MHD data.
    """

    path = "mhd" + "/MHD-" + qid

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
    Write dummy MHD input.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """
    n_modes = 2
    npsi    = 6

    nmodes    = np.array([1, 2])
    mmodes    = np.array([3, 4])
    amplitude = np.array([0.1, 2])
    omega     = np.array([1, 1.5])
    alpha     = np.ones((n_modes,npsi))
    phi       = np.ones((n_modes,npsi))
    write_hdf5(fn, nmodes, mmodes, amplitude, omega, alpha, phi,
               desc="Dummy")
