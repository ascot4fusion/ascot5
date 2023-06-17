"""Input containing MHD eigenmodes.

This input is used in simulations studying particle response to MHD.

Using this data in a simulation requires that boozer data, containing
mapping between real-space and Boozer coordinates, is present.

This module reads and writes both MHD_STAT and MHD_NONSTAT data. The only
difference between these two is the time-dependency of the eigenfunctions,
that is, does the mode amplitude (alpha or Phi) depend only on psi or psi
and time. MHD_STAT assumes only psi dependency, making the interpolation
*much* faster. So don't use MHD_NONSTAT unless you really have to.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

class MHD_STAT(DataGroup):
    """Stationary MHD eigenfunctions.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["alpha"] = np.transpose(out["alpha"], (1,0) )
        out["phi"]   = np.transpose(out["phi"],   (1,0) )
        return out

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

    @staticmethod
    def write_hdf5(fn, nmode, nmodes, mmodes, amplitude, omega, phase, alpha,
                   phi, nrho, rhomin, rhomax, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        nmode : int
            Number of modes.
        nmodes : array_like (nmode,)
            Mode n (toroidal) numbers.
        mmodes : array_like (nmode,)
            Mode m (poloidal) numbers.
        amplitude : array_like (nmode,)
            Mode amplitudies.
        omega : array_like (nmode,)
            Mode frequencies [rad/s].
        omega : array_like (nmode,)
            Mode phases [rad].
        alpha : array_like (nrho, nmode)
            Magnetic perturbation eigenfunctions.
        phi : array_like (nrho, nmode)
            Electric perturbation eigenfunctions.
        nrho : int
            Number of rho grid points.
        rhomin : float
            Minimum value in rho grid.
        rhomax : float
            Maximum value in rho grid.
        desc : str, optional
            Input's description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if nmodes.size != nmode:
            raise ValueError("Shape of nmodes is inconsistent with nmode")
        if mmodes.size != nmode:
            raise ValueError("Shape of mmodes is inconsistent with nmode")
        if alpha.shape != (nrho,nmode):
            raise ValueError("alpha has inconsistent shape")
        if phi.shape   != (nrho,nmode):
            raise ValueError("phi has inconsistent shape")
        alpha = np.transpose(alpha, (1,0) )
        phi   = np.transpose(phi,   (1,0) )

        parent = "mhd"
        group  = "MHD_STAT"
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

            g.create_dataset("alpha", (nmode,nrho), data=alpha, dtype="f8")
            g.create_dataset("phi",   (nmode,nrho), data=phi,   dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This method writes two modes with constant eigenfunctions.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return MHD_STAT.write_hdf5(
            fn=fn, nmode=2, nmodes=np.array([1, 2]), mmodes=np.array([3, 4]),
            amplitude=np.array([0.1, 2]), omega=np.array([1, 1.5]),
            phase=np.array([0, 3.141/4]), alpha=np.ones((6,2)),
            phi=np.ones((6,2)), nrho=6, rhomin=0, rhomax=1, desc="DUMMY")

class MHD_NONSTAT(DataGroup):
    """Time-dependent MHD eigenfunctions.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["alpha"] = np.transpose(out["alpha"], (2,1,0) )
        out["phi"]   = np.transpose(out["phi"],   (2,1,0) )
        return out

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

    @staticmethod
    def write_hdf5(fn, nmode, nmodes, mmodes, amplitude, omega, phase, alpha,
                   phi, nrho, rhomin, rhomax, ntime, tmin, tmax,
                   desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        nmode : int
            Number of modes.
        nmodes : array_like (nmode,)
            Mode n (toroidal) numbers.
        mmodes : array_like (nmode,)
            Mode m (poloidal) numbers.
        amplitude : array_like (nmode,)
            Mode amplitudies.
        omega : array_like (nmode,)
            Mode frequencies [rad/s].
        omega : array_like (nmode,)
            Mode phases [rad].
        alpha : array_like (nrho, ntime, nmode)
            Magnetic perturbation eigenfunctions.
        phi : array_like (nrho, ntime, nmode)
            Electric perturbation eigenfunctions.
        nrho : int
            Number of rho grid points.
        rhomin : float
            Minimum value in rho grid.
        rhomax : float
            Maximum value in rho grid.
        ntime : int, optional
            Number of time grid points.
        tmin : float
            Minimum value in time grid.
        tmax : float, optional
            Maximum value in time grid.
        desc : str, optional
            Input's description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if nmodes.size != nmode:
            raise ValueError("Shape of nmodes is inconsistent with nmode")
        if mmodes.size != nmode:
            raise ValueError("Shape of mmodes is inconsistent with nmode")
        if alpha.shape != (nrho,ntime,nmode):
            raise ValueError("alpha has inconsistent shape")
        if phi.shape   != (nrho,ntime,nmode):
            raise ValueError("phi has inconsistent shape")
        alpha = np.transpose(alpha, (2,1,0) )
        phi   = np.transpose(phi,   (2,1,0) )

        parent = "mhd"
        group  = "MHD_NONSTAT"
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

            g.create_dataset("ntime", (1,), data=ntime, dtype="i4")
            g.create_dataset("tmin",  (1,), data=tmin,  dtype="f8")
            g.create_dataset("tmax",  (1,), data=tmax,  dtype="f8")
            g.create_dataset("alpha", (nmode,ntime,nrho), data=alpha,dtype="f8")
            g.create_dataset("phi",   (nmode,ntime,nrho), data=phi,  dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This method writes two modes with constant eigenfunctions.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return MHD_NONSTAT.write_hdf5(
            fn=fn, nmode=2, nmodes=np.array([1, 2]), mmodes=np.array([3, 4]),
            amplitude=np.array([0.1, 2]), omega=np.array([1, 1.5]),
            phase=np.array([0, 3.141/4]), alpha=np.ones((6,3,2)),
            phi=np.ones((6,3,2)), nrho=6, rhomin=0, rhomax=1, ntime=3,
            tmin=0, tmax=1, desc="DUMMY")
