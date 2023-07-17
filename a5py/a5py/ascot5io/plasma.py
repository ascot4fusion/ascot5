"""Input data representing plasma background species.

Plasma input is required for simulations with collisions or collisional
diagnostics enabled. This includes BBNBI5 and ASCOT-BMC.
"""
import h5py
import numpy as np

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class plasma_1D(DataGroup):
    """Plasma profiles that have only radial dependency.
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
                if key in ["nion", "nrho"]:
                    out[key] = int(out[key])

        out["idensity"] = np.transpose(out["idensity"])
        return out

    def plot(self, pls=None):
        import matplotlib.pyplot as plt

        if pls is None:
            pls = self.read()

        plotyy=False # This didn't come out nice, but keep it just in case.

        if plotyy:
            fig, ax1 = plt.subplots()

            color = 'tab:red'
            ax1.set_xlabel('$\\rho $')
            ax1.set_ylabel('Density (10$^{19}$/m$^3$)', color=color)
            ax1.plot(pls['rho'],pls['edensity']*1e-19,'-' ,color=color,label='n$_e$')
            ax1.plot(pls['rho'],pls['idensity']*1e-19,'--',label='n$_i$')
            ax1.tick_params(axis='y', labelcolor=color)

            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            color = 'tab:blue'
            ax2.set_ylabel('Temperature(eV)', color=color)  # we already handled the x-label with ax1
            ax2.plot(pls['rho'],pls['itemperature'],'-' ,color=color,label='T$_i$')
            ax2.plot(pls['rho'],pls['etemperature'],'--',label='T$_e$')
            ax2.tick_params(axis='y', labelcolor=color)

            ax1.legend()
            ax2.legend()

            fig.tight_layout()  # otherwise the right y-label is slightly clipped

            plt.show()

        else:

            plt.subplot(2,1,1)
            plt.plot(pls['rho'],pls['edensity']*1e-19,'-',label='n$_e$')
            plt.plot(pls['rho'],pls['idensity']*1e-19,'--',label='n$_i$')
            plt.legend()
            plt.ylabel('Density (10$^{19}$/m$^3$)')

            plt.subplot(2,1,2)
            plt.plot(pls['rho'],pls['etemperature'],'-',label='T$_e$')
            plt.plot(pls['rho'],pls['itemperature'],'--',label='T$_i$')
            plt.legend()
            plt.ylabel('Temperature (eV)')

            plt.xlabel('$\\rho $')

            plt.show()

    @staticmethod
    def write_hdf5(fn, nrho, nion, anum, znum, mass, charge, rho,
                   edensity, etemperature, idensity, itemperature, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Path to hdf5 file.
        nrho : int
            Number of rho grid points.
        nion : int
            Number of ion species.
        anum : array_like (nion,1)
            Ion species atomic mass number
        znum : array_like (nion,1)
            Ion species charge number.
        mass : array_like (nion,1)
            Ion species mass [amu].
        charge : array_like (nion,1)
            Ion species charge [e].
        rho : array_like (nrho,1)
            rho grid, doesn't have to be uniform.
        edensity : array_like (nrho,1)
            Electron density [m^-3].
        etemperature : array_like (nrho,1)
            Electron temperature [eV].
        idensity : array_like (nrho,nion)
            Ion density [m^-3].
        itemperature : array_like (nrho,1)
            Ion temperature [ev].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if etemperature.size != nrho:
            raise ValueError("Invalid size for electron temperature.")
        if itemperature.size != nrho:
            raise ValueError("Invalid size for ion temperature.")
        if edensity.size != nrho:
            raise ValueError("Invalid size for electron density.")
        if idensity.shape != (nrho,nion):
            raise ValueError("Invalid size for ion density.")

        idensity = np.transpose(idensity)

        parent = "plasma"
        group  = "plasma_1D"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nion',   (1,1),    data=nion,   dtype='i4')
            g.create_dataset('nrho',   (1,1),    data=nrho,   dtype='i4')
            g.create_dataset('znum',   (nion,1), data=znum,   dtype='i4')
            g.create_dataset('anum',   (nion,1), data=anum,   dtype='i4')
            g.create_dataset('charge', (nion,1), data=charge, dtype='i4')
            g.create_dataset('mass',   (nion,1), data=mass,   dtype='f8')
            g.create_dataset('rho',    (nrho,1), data=rho,    dtype='f8')

            g.create_dataset('etemperature', (nrho,1),    data=etemperature,
                             dtype='f8')
            g.create_dataset('edensity',     (nrho,1),    data=edensity,
                             dtype='f8')
            g.create_dataset('itemperature', (nrho,1),    data=itemperature,
                             dtype='f8')
            g.create_dataset('idensity',     (nion,nrho), data=idensity,
                             dtype='f8')

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is an uniform hydrogen plasma.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return plasma_1D.write_hdf5(
            fn=fn, nrho=3, nion=1, znum=np.array([1]), anum=np.array([1]),
            mass=np.array([1]), charge=np.array([1]),
            rho=np.array([0, 0.5, 100]), edensity=1e20*np.ones((3,1)),
            etemperature=1e3*np.ones((3,1)), idensity=1e20*np.ones((3,1)),
            itemperature=1e20*np.ones((3,1)), desc="DUMMY")

class plasma_1DS(DataGroup):
    """Same input as :class:`plasma_1D` but interpolated with splines.

    Spline interpolation could yield negative value even though the inputs
    were all positive. To avoid this, ASCOT5 internally takes a logarithm
    of the input data before constructing the splines. The actual values
    are then computed after the spline-evaluation. The logarithmic also
    helps to interpolate the values at the edge more accurately.
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
                if key in ["nion", "nrho"]:
                    out[key] = int(out[key])

        out["idensity"] = np.transpose(out["idensity"])
        return out

    @staticmethod
    def write_hdf5(fn, nrho, nion, anum, znum, mass, charge, rho,
                   edensity, etemperature, idensity, itemperature, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Path to hdf5 file.
        nrho : int
            Number of rho grid points.
        nion : int
            Number of ion species.
        anum : array_like (nion,1)
            Ion species atomic mass number
        znum : array_like (nion,1)
            Ion species charge number.
        mass : array_like (nion,1)
            Ion species mass [amu].
        charge : array_like (nion,1)
            Ion species charge [e].
        rho : array_like (nrho,1)
            rho grid, doesn't have to be uniform.
        edensity : array_like (nrho,1)
            Electron density [m^-3].
        etemperature : array_like (nrho,1)
            Electron temperature [eV].
        idensity : array_like (nrho,nion)
            Ion density [m^-3].
        itemperature : array_like (nrho,1)
            Ion temperature [ev].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if etemperature.size != nrho:
            raise ValueError("Invalid size for electron temperature.")
        if itemperature.size != nrho:
            raise ValueError("Invalid size for ion temperature.")
        if edensity.size != nrho:
            raise ValueError("Invalid size for electron density.")
        if idensity.shape != (nrho,nion):
            raise ValueError("Invalid size for ion density.")

        idensity = np.transpose(idensity)

        parent = "plasma"
        group  = "plasma_1DS"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nion',   (1,1),    data=nion,   dtype='i4')
            g.create_dataset('nrho',   (1,1),    data=nrho,   dtype='i4')
            g.create_dataset('znum',   (nion,1), data=znum,   dtype='i4')
            g.create_dataset('anum',   (nion,1), data=anum,   dtype='i4')
            g.create_dataset('charge', (nion,1), data=charge, dtype='i4')
            g.create_dataset('mass',   (nion,1), data=mass,   dtype='f8')
            g.create_dataset('rho',    (nrho,1), data=rho,    dtype='f8')

            g.create_dataset('etemperature', (nrho,1),    data=etemperature,
                             dtype='f8')
            g.create_dataset('edensity',     (nrho,1),    data=edensity,
                             dtype='f8')
            g.create_dataset('itemperature', (nrho,1),    data=itemperature,
                             dtype='f8')
            g.create_dataset('idensity',     (nion,nrho), data=idensity,
                             dtype='f8')

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is an uniform hydrogen plasma.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return plasma_1DS.write_hdf5(
            fn=fn, nrho=3, nion=1, znum=np.array([1]), anum=np.array([1]),
            mass=np.array([1]), charge=np.array([1]),
            rho=np.array([0, 0.5, 100]), edensity=1e20*np.ones((3,1)),
            etemperature=1e3*np.ones((3,1)), idensity=1e20*np.ones((3,1)),
            itemperature=1e20*np.ones((3,1)), desc="DUMMY")
