"""
MHD input IO.

File: mhd.py
"""
import h5py
import numpy as np

from . ascot5file import add_group
from . ascot5data import AscotData

def write_hdf5(fn, nmode, nmodes, mmodes, amplitude, omega, alpha, phi,
               npsi, psimin, psimax, ntime=4, tmin=-100, tmax=100,
               desc=None):
    """
    Write MHD input to HDF5 file.

    Giving time grid is optional, but if it is not given, default values are
    used as ASCOT5 requires time-dependent profiles.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        nmode : int <br>
            Number of modes.
        nmodes : array_like (nmode,) <br>
            Mode n (toroidal) numbers.
        mmodes : array_like (nmode,) <br>
            Mode m (poloidal) numbers.
        amplitude : array_like (nmode,) <br>
            Mode amplitudies.
        omega : array_like (nmode,) <br>
            Mode frequencies [rad/s].
        alpha : array_like (npsi, ntime, nmode) <br>
            Mode alpha. If no time grid, the shape should be (npsi, nmode).
        phi : array_like (npsi, ntime, nmode) <br>
            Mode phi. If no time grid, the shape should be (npsi, nmode).
        npsi : int <br>
            Number of psi grid points.
        psimin : float <br>
            Minimum value in psi grid.
        psimax : float <br>
            Maximum value in psi grid.
        ntime : int, optional <br>
            Number of time grid points.
        tmin : float, optional <br>
            Minimum value in time grid. Must be given for non-stationary input.
        tmax : float, optional <br>
            Maximum value in time grid. Must be given for non-stationary input.
        desc : str, optional <br>
            Input's description.

    Returns:
        Name of the new input that was written.
    """
    assert nmodes.size == nmode
    assert mmodes.size == nmode

    if alpha.ndim == 2:
        # Add quasi time-dependency.
        alpha = np.transpose( np.tile(alpha, (ntime,1,1)), (1,0,2) )
        phi   = np.transpose( np.tile(phi,   (ntime,1,1)), (1,0,2) )

    assert alpha.shape == (npsi,ntime,nmode)
    assert phi.shape   == (npsi,ntime,nmode)

    alpha = np.transpose(alpha, (2,1,0) )
    phi   = np.transpose(phi,   (2,1,0) )

    parent = "mhd"
    group  = "MHD"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("nmode",  (1,), data=nmode,  dtype="i4")
        g.create_dataset("npsi",   (1,), data=npsi,   dtype="i4")
        g.create_dataset("psimin", (1,), data=psimin, dtype="f8")
        g.create_dataset("psimax", (1,), data=psimax, dtype="f8")
        g.create_dataset("ntime",  (1,), data=ntime,  dtype="i4")
        g.create_dataset("tmin",   (1,), data=tmin,   dtype="f8")
        g.create_dataset("tmax",   (1,), data=tmax,   dtype="f8")

        g.create_dataset("nmodes", (nmode,), data=nmodes, dtype="i4")
        g.create_dataset("mmodes", (nmode,), data=mmodes, dtype="i4")

        g.create_dataset("amplitude", (nmode,), data=amplitude, dtype="f8")
        g.create_dataset("omega",     (nmode,), data=omega,     dtype="f8")

        g.create_dataset("alpha", (nmode,ntime,npsi), data=alpha, dtype="f8")
        g.create_dataset("phi",   (nmode,ntime,npsi), data=phi,   dtype="f8")

    return gname


def write_hdf5_dummy(fn, desc="Dummy"):
    """
    Write dummy MHD input.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """
    nmode = 2
    npsi  = 6

    nmodes    = np.array([1, 2])
    mmodes    = np.array([3, 4])
    amplitude = np.array([0.1, 2])
    omega     = np.array([1, 1.5])
    alpha     = np.ones((npsi,nmode))
    phi       = np.ones((npsi,nmode))
    psimin    = 0
    psimax    = 1
    write_hdf5(fn, nmode, nmodes, mmodes, amplitude, omega, alpha, phi,
               npsi, psimin, psimax,
               desc=desc)


def read_hdf5(fn, qid):
    """
    Read MHD input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "mhd/MHD_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    out["alpha"] = np.transpose(out["alpha"], (2,1,0) )
    out["phi"]   = np.transpose(out["phi"],   (2,1,0) )
    return out


class MHD(AscotData):
    """
    Object representing MHD data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())
