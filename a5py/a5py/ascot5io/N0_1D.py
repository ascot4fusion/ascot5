"""
Constant-on-flux-surfaces tokamak neutral HDF5 IO

File: N0_1D.py
"""
import numpy as np
import h5py

from . ascot5file import add_group
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, rhomin, rhomax, nrho,
               nspecies, anum, znum, density, temperature, maxwellian=1,
               desc=None):
    """
    Write 1D neutral input in HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        rhomin : float <br>
            Minimum value in rho grid [1].
        rhomax : float <br>
            Maximum value in rho grid [1].
        nrho : int <br>
            Number of rho grid points.
        nspecies : int <br>
            Number of neutral species.
        anum : array_like (nspecies,1) <br>
            Neutral species' atomic mass number.
        znum array_like (nspecies,1) <br>
            Neutral species' charge number.
        density array_like (nrho,nspecies) <br>
            Neutral species-wise density [m^-3].
        temperature array_like (nrho,nspecies) <br>
            Neutral species-wise temperature [eV].
        maxwellian array_like (nspecies,1) <br> :
            Whether species distribution is Maxwellian (1) of monoenergetic (0)
        desc : str, optional <br>
            Input description.

    Returns:
        Name of the new input that was written.
    """
    assert density.shape    == (nrho,nspecies)
    assert temperature.size == nrho
    assert anum.size == nspecies
    assert znum.size == nspecies
    assert (maxwellian == 1 or maxwellian.size == nspecies)

    parent = "neutral"
    group  = "N0_1D"
    gname  = ""

    # Transpose n0 from (rho, spec) to (spec, rho)
    density     = np.transpose(density)

    if maxwellian == 1:
        maxwellian = np.ones( (int(nspecies),1) )

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("rhomin",   (1,), data=rhomin,   dtype="f8")
        g.create_dataset("rhomax",   (1,), data=rhomax,   dtype="f8")
        g.create_dataset("nrho",     (1,), data=nrho,     dtype="i4")
        g.create_dataset("nspecies", (1,), data=nspecies, dtype="i4")

        g.create_dataset("anum",        (nspecies,),  data=anum,
                         dtype="i4")
        g.create_dataset("znum",        (nspecies,),  data=znum,
                         dtype="i4")
        g.create_dataset("maxwellian",  (nspecies,),  data=maxwellian,
                         dtype="i4")

        g.create_dataset("density",     (nspecies,nrho), data=density,
                         dtype="f8")
        g.create_dataset("temperature", (nspecies,nrho), data=temperature,
                         dtype="f8")

    return gname


def read_hdf5(fn, qid):
    """
    Read 1D neutral input from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

    Returns:
        Dictionary containing input data.
    """

    path = "neutral/N0_1D_" + qid

    out = {}
    with h5py.File(fn,"r") as f:
        for key in f[path]:
            out[key] = f[path][key][:]

    out["density"]     = np.transpose(out["density"],     (1,0))
    out["temperature"] = np.transpose(out["temperature"], (1,0))

    return out


def write_hdf5_dummy(fn, desc="Dummy"):
    N0rhomin = 0
    N0rhomax = 2
    N0nrho   = 100
    N0spec   = 1
    N0anum   = np.array([1])
    N0znum   = np.array([1])
    N0dens   = 5e+16*np.ones( (N0nrho,N0spec) )
    N0temp   = 1e+3 *np.ones( (N0nrho,N0spec) )
    write_hdf5(fn,
               N0rhomin, N0rhomax, N0nrho,
               N0spec, N0anum, N0znum,
               N0dens, N0temp,
               desc=desc)


class N0_1D(AscotData):
    """
    Object representing N0_1D data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)
