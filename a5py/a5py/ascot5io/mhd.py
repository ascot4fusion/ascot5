"""
MHD input IO.

This module reads and writes both MHD_STAT and MHD_NONSTAT data. The only
difference between these two is the time-dependency of the eigenfunctions,
that is, does the mode amplitude (alpha or Phi) depend only on psi or psi
and time. MHD_STAT assumes only psi dependency, making the interpolation
*much* faster. So don't use MHD_NONSTAT unless you really have to.

File: mhd.py
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

def write_hdf5(fn, nmode, nmodes, mmodes, amplitude, omega, phase, alpha, phi,
               nrho, rhomin, rhomax, ntime=None, tmin=None, tmax=None,
               desc=None):
    """
    Write MHD input to HDF5 file.

    This module writes both stationary and non-stationary MHD data. The latter
    is used when the time-grid is provided.

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
        omega : array_like (nmode,) <br>
            Mode phases [rad].
        alpha : array_like (nrho, ntime, nmode) <br>
            Magnetic perturbation eigenfunctions. If no time grid, the shape
            should be (nrho, nmode).
        phi : array_like (nrho, ntime, nmode) <br>
            Electric perturbation eigenfunctions. If no time grid, the shape
            should be (nrho, nmode).
        nrho : int <br>
            Number of rho grid points.
        rhomin : float <br>
            Minimum value in rho grid.
        rhomax : float <br>
            Maximum value in rho grid.
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

    assert (ntime is None and tmin is None and tmax is None) or \
    (ntime is not None and tmin is not None and tmax is not None)

    if ntime is None:
        assert alpha.shape == (nrho,nmode)
        assert phi.shape   == (nrho,nmode)
        alpha = np.transpose(alpha, (1,0) )
        phi   = np.transpose(phi,   (1,0) )
    else:
        assert alpha.shape == (nrho,ntime,nmode)
        assert phi.shape   == (nrho,ntime,nmode)
        alpha = np.transpose(alpha, (2,1,0) )
        phi   = np.transpose(phi,   (2,1,0) )

    parent = "mhd"
    group  = "MHD_STAT" if ntime is None else "MHD_NONSTAT"
    gname  = ""

    with h5py.File(fn, "a") as f:
        g = add_group(f, parent, group, desc=desc)
        gname = g.name.split("/")[-1]

        g.create_dataset("nmode",  (1,), data=nmode,  dtype="i4")
        g.create_dataset("nrho",   (1,), data=nrho,   dtype="i4")
        g.create_dataset("rhomin", (1,), data=rhomin, dtype="f8")
        g.create_dataset("rhomax", (1,), data=rhomax, dtype="f8")

        g.create_dataset("nmodes", (nmode,), data=nmodes, dtype="i4")
        g.create_dataset("mmodes", (nmode,), data=mmodes, dtype="i4")

        g.create_dataset("amplitude", (nmode,), data=amplitude, dtype="f8")
        g.create_dataset("omega",     (nmode,), data=omega,     dtype="f8")
        g.create_dataset("phase",     (nmode,), data=phase,     dtype="f8")

        if ntime is None:
            g.create_dataset("alpha", (nmode,nrho), data=alpha, dtype="f8")
            g.create_dataset("phi",   (nmode,nrho), data=phi,   dtype="f8")
        else:
            g.create_dataset("ntime", (1,), data=ntime, dtype="i4")
            g.create_dataset("tmin",  (1,), data=tmin,  dtype="f8")
            g.create_dataset("tmax",  (1,), data=tmax,  dtype="f8")
            g.create_dataset("alpha", (nmode,ntime,nrho), data=alpha,dtype="f8")
            g.create_dataset("phi",   (nmode,ntime,nrho), data=phi,  dtype="f8")

    return gname


def write_hdf5_dummy(fn, desc="Dummy"):
    """
    Write dummy MHD input.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
    """
    nmode = 2
    nrho  = 6

    nmodes    = np.array([1, 2])
    mmodes    = np.array([3, 4])
    amplitude = np.array([0.1, 2])
    omega     = np.array([1, 1.5])
    phase     = np.array([0, 3.141/4])
    alpha     = np.ones((nrho,nmode))
    phi       = np.ones((nrho,nmode))
    rhomin    = 0
    rhomax    = 1
    return write_hdf5(
        fn, nmode, nmodes, mmodes, amplitude, omega, phase, alpha, phi,
        nrho, rhomin, rhomax,
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

    path = "mhd/MHD_STAT_" + qid
    isstat = True

    out = {}
    with h5py.File(fn,"r") as f:
        if path not in f:
            path = "mhd/MHD_NONSTAT_" + qid
            isstat = False
        for key in f[path]:
            out[key] = f[path][key][:]

    if isstat:
        out["alpha"] = np.transpose(out["alpha"], (1,0) )
        out["phi"]   = np.transpose(out["phi"],   (1,0) )
    else:
        out["alpha"] = np.transpose(out["alpha"], (2,1,0) )
        out["phi"]   = np.transpose(out["phi"],   (2,1,0) )
    return out


class MHD(DataGroup):
    """
    Object representing MHD data.
    """

    def read(self):
        return read_hdf5(self._file, self.get_qid())

    def isstat(self):
        """
        Is the data STAT or NONSTAT.
        """
        path = "mhd/MHD_STAT_" + self.get_qid()
        with h5py.File(self._file, "r") as f:
            if path not in f:
                return False
        return True


    def plot_amplitude(self, amplitude="alpha", mode=None, ax=None):
        """
        Plot radial profile of all (or given mode) amplitudes.

        For NONSTAT the profiles are shown at earliest time in dataset.
        TODO: Animate the NONSTAT to show profiles as a function of time.

        Args:
            amplitude : str <br>
                Which amplitude is plotted: "alpha" (magnetic) or "phi" (electric)
            mode : tuple(int,int) <br>
                Mode (n,m) to be plotted or None to plot all modes.
            ax : axes <br>
                Axes where plotting is done or None to create and show new fig.
        """
        import matplotlib.pyplot as plt

        data = self.read()
        isstat = self.isstat()

        axes = ax
        if ax is None:
            fig  = plt.figure()
            axes = fig.add_subplot(1,1,1)

        axes.set_xlabel(r"$\rho$")
        axes.set_xlim(0,1)
        if amplitude == "alpha":
            axes.set_ylabel("Tm")
        else:
            axes.set_ylabel("V")

        amplitude = data[amplitude]
        rho = np.linspace(data["rhomin"], data["rhomax"], int(data["nrho"]))

        for n in range(int(data["nmode"])):
            if mode is None or \
               ( mode[0] == data["nmodes"][n] and mode[1] == data["mmodes"][n] ):
                if isstat:
                    axes.plot(rho, amplitude[:,n] * data["amplitude"][n])
                else:
                    #t = data["tmin"]
                    #tgrid = np.linspace(data["tmin"], data["tmax"], data["ntime"])
                    axes.plot(rho, amplitude[n,0,:] * data["amplitude"][n])

        if ax is None:
            plt.show()
