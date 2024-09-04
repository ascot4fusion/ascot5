"""Input data representing plasma background species.

Plasma input is required for simulations with collisions or collisional
diagnostics enabled. This includes BBNBI5 and ASCOT-BMC.
"""
import h5py
import numpy as np

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

import a5py.routines.plotting as a5plt
from a5py.physlib.species import speciesdict

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

    def plot_radial(self, rholim=None, densitylim=None, temperaturelim=None,
                    axes=None):
        """Plot plasma profiles as a function of rho.

        Parameters
        ----------
        rholim : [float, float], optional
            Limits on x-axis.
        densitylim : [float, float, float], optional
            Limits on the first y axis where the middle value is when the scale
            changes from logarithmic to linear [m^-3].
        temperaturelim : [float, float, float], optional
            Limits on the second y axis where the middle value is when the scale
            changes from logarithmic to linear [eV].
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        """
        def formatspec(s, a, z, q):
            """Get species name in '^anum_znum A^charge' format"""
            return "$_{{{:d}}}^{{{:d}}}{}^{{{}}}$".format(a,z,s,"+" + str(q))

        pls = self.read()
        ndens = [pls['edensity']]
        y1legends = ["$n_e$"]
        for i in range(pls["nion"]):
            ndens.append(pls["idensity"][:,i])

            a = int(pls["anum"][i])
            z = int(pls["znum"][i])
            q = int(pls["charge"][i])
            y1legends.append(formatspec("A", a, z, q))
            for s, d in speciesdict.items():
                if d[0] == a and d[1] == z:
                    s = ''.join(c for c in s if not c.isnumeric())
                    y1legends[i+1] = formatspec(s, a, z, int(d[2].v))
                    break

        if rholim is None: rholim = [pls['rho'][0], pls['rho'][-1]]
        if densitylim is None:
            densitylim = [np.amin(pls["idensity"]), 1e19,
                          np.amax(pls["edensity"])]
        if temperaturelim is None:
            temperaturelim = [np.amin(pls["itemperature"]), 1e3,
                              np.amax(pls["itemperature"])]

        a5plt.radialprofile(
            pls['rho'], ndens, y2=[pls['etemperature'], pls['itemperature']],
            xlim=rholim, y1lim=densitylim, y2lim=temperaturelim,
            xlabel=r"$\rho$", y1label=r"Density [m$^{-3}$]",
            y1legends=y1legends, y2label=r"Temperature [eV]",
            y2legends=[r"$T_e$", r"$T_i$"], axes=axes)

    @staticmethod
    def write_hdf5(fn, nrho, nion, anum, znum, mass, charge, rho, vtor,
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
        vtor : array_like (nrho,1)
            Plasma rotation [rad/s].
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
        if vtor.size != nrho:
            raise ValueError("Invalid size for toroidal rotation.")
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
            g.create_dataset('vtor',   (nrho,1), data=vtor,   dtype='f8')

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
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is an uniform hydrogen plasma.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return {"nrho":3, "nion":1, "znum":np.array([1]), "anum":np.array([1]),
                "mass":np.array([1]), "charge":np.array([1]),
                "rho":np.array([0, 0.5, 100]), "vtor":np.zeros((3,1)),
                "edensity":1e20*np.ones((3,1)),
                "etemperature":1e3*np.ones((3,1)),
                "idensity":1e20*np.ones((3,1)),
                "itemperature":1e20*np.ones((3,1))}

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

    def plot_radial(self, rholim=None, densitylim=None, temperaturelim=None,
                    axes=None):
        """Plot plasma profiles as a function of rho.

        Parameters
        ----------
        rholim : [float, float], optional
            Limits on x-axis.
        densitylim : [float, float, float], optional
            Limits on the first y axis where the middle value is when the scale
            changes from logarithmic to linear [m^-3].
        temperaturelim : [float, float, float], optional
            Limits on the second y axis where the middle value is when the scale
            changes from logarithmic to linear [eV].
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        """
        def formatspec(s, a, z, q):
            """Get species name in '^anum_znum A^charge' format"""
            return "$_{{{:d}}}^{{{:d}}}{}^{{{}}}$".format(a,z,s,"+" + str(q))

        pls = self.read()
        ndens = [pls['edensity']]
        y1legends = ["$n_e$"]
        for i in range(pls["nion"]):
            ndens.append(pls["idensity"][:,i])

            a = int(pls["anum"][i])
            z = int(pls["znum"][i])
            q = int(pls["charge"][i])
            y1legends.append(formatspec("A", a, z, q))
            for s, d in speciesdict.items():
                if d[0] == a and d[1] == z:
                    s = ''.join(c for c in s if not c.isnumeric())
                    y1legends[i+1] = formatspec(s, a, z, int(d[2].v))
                    break

        if rholim is None: rholim = [pls['rhomin'][0], pls['rhomax'][0]]
        if densitylim is None:
            densitylim = [np.amin(pls["idensity"]), 1e19,
                          np.amax(pls["edensity"])]
        if temperaturelim is None:
            temperaturelim = [np.amin(pls["itemperature"]), 1e3,
                              np.amax(pls["itemperature"])]

        rho = np.linspace(pls['rhomin'], pls['rhomax'], pls['nrho'])
        a5plt.radialprofile(
            rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
            xlim=rholim, y1lim=densitylim, y2lim=temperaturelim,
            xlabel=r"$\rho$", y1label=r"Density [m$^{-3}$]",
            y1legends=y1legends, y2label=r"Temperature [eV]",
            y2legends=[r"$T_e$", r"$T_i$"], axes=axes)

    @staticmethod
    def write_hdf5(fn, nrho, nion, anum, znum, mass, charge, rhomin, rhomax,
                   vtor, edensity, etemperature, idensity, itemperature,
                   desc=None):
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
        rhomin : float
            Minimum rho grid value.
        rhomax : float
            Maximum rho grid value.
        vtor : array_like (nrho,1)
            Plasma rotation [rad/s].
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
        if vtor.size != nrho:
            raise ValueError("Invalid size for toroidal flow.")
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
            g.create_dataset('rhomin', (1,1),    data=rhomin, dtype='f8')
            g.create_dataset('rhomax', (1,1),    data=rhomax, dtype='f8')
            g.create_dataset('vtor',   (nrho,1), data=vtor,   dtype='f8')
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
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is an uniform hydrogen plasma.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return {"nrho":3, "nion":1, "znum":np.array([1]), "anum":np.array([1]),
                "mass":np.array([1]), "charge":np.array([1]),
                "rhomin":0, "rhomax":100, "vtor":np.zeros((3,1)),
                "edensity":1e20*np.ones((3,1)),
                "etemperature":1e3*np.ones((3,1)),
                "idensity":1e20*np.ones((3,1)),
                "itemperature":1e3*np.ones((3,1))}

class plasma_1Dt(DataGroup):
    """Time-dependent 1D plasma profiles.
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
                if key in ["nion", "nrho", "ntime"]:
                    out[key] = int(out[key])

        return out

    @staticmethod
    def write_hdf5(fn, nrho, ntime, nion, anum, znum, mass, charge, rho, time,
                   vtor, edensity, etemperature, idensity, itemperature,
                   desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Path to hdf5 file.
        nrho : int
            Number of rho grid points.
        ntime : int
            Number of time grid points.
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
        time : array_like (nrho,1)
            time grid, doesn't have to be uniform.
        vtor : array_like (ntime,nrho)
            Plasma rotation [rad/s].
        edensity : array_like (ntime,nrho)
            Electron density [m^-3].
        etemperature : array_like (ntime,nrho)
            Electron temperature [eV].
        idensity : array_like (ntime,nion,nrho)
            Ion density [m^-3].
        itemperature : array_like (ntime,nrho)
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
        if vtor.shape != (ntime,nrho):
            raise ValueError("Invalid shape for toroidal flow.")
        if etemperature.shape != (ntime,nrho):
            raise ValueError("Invalid shape for electron temperature.")
        if itemperature.shape != (ntime,nrho):
            raise ValueError("Invalid shape for ion temperature.")
        if edensity.shape != (ntime,nrho):
            raise ValueError("Invalid shape for electron density.")
        if idensity.shape != (ntime,nion,nrho):
            raise ValueError("Invalid shape for ion density.")

        parent = "plasma"
        group  = "plasma_1Dt"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nion',   (1,1),     data=nion,   dtype='i4')
            g.create_dataset('nrho',   (1,1),     data=nrho,   dtype='i4')
            g.create_dataset('ntime',  (1,1),     data=ntime,  dtype='i4')
            g.create_dataset('znum',   (nion,1),  data=znum,   dtype='i4')
            g.create_dataset('anum',   (nion,1),  data=anum,   dtype='i4')
            g.create_dataset('charge', (nion,1),  data=charge, dtype='i4')
            g.create_dataset('mass',   (nion,1),  data=mass,   dtype='f8')
            g.create_dataset('rho',    (nrho,1),  data=rho,    dtype='f8')
            g.create_dataset('time',   (ntime,1), data=time,   dtype='f8')
            g.create_dataset('vtor',   (ntime,nrho), data=vtor, dtype='f8')
            g.create_dataset('etemperature', (ntime,nrho),
                             data=etemperature, dtype='f8')
            g.create_dataset('edensity',     (ntime,nrho),
                             data=edensity, dtype='f8')
            g.create_dataset('itemperature', (ntime,nrho),
                             data=itemperature, dtype='f8')
            g.create_dataset('idensity',     (ntime,nion,nrho),
                             data=idensity, dtype='f8')

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is an uniform hydrogen plasma.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return {"nrho":3, "ntime":4, "nion":1, "znum":np.array([1]),
                "anum":np.array([1]), "mass":np.array([1]),
                "charge":np.array([1]), "rho":np.array([0, 0.5, 100]),
                "time":np.array([0, 0.2, 0.4, 0.6]), "vtor":np.zeros((4,3)),
                "edensity":1e20*np.ones((4,3)),
                "etemperature":1e3*np.ones((4,3)),
                "idensity":1e20*np.ones((4,1,3)),
                "itemperature":1e3*np.ones((4,3))}
