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

    def plot_radial(self,
                    rholim=None,
                    rholim_sol=None,
                    densitylim=None,
                    densitylim_sol=None,
                    temperaturelim=None,
                    temperaturelim_sol=None,
                    axes=None,
                    plot_sol_separately=False,
                    ):
        """Plot plasma profiles as a function of rho. If ``plot_sol_separately``
        is True, then the core and scrape-off layer are plotted separately.

        Parameters
        ----------
        rholim : [float, float], optional
            Limits on x-axis. If plot_sol_separately is True, then the limits
            are applied to the left side of x-axis.
        rholim_sol : [float, float], optional
            Limits on the right side of x-axis when ``plot_sol_separately`` is
            True.
        densitylim : [float, float, float] or [float, float], optional
            Limits on the first y axis where the middle value is when the scale
            changes from logarithmic to linear [m^-3]. If ``plot_sol_separately``
            is True, only the first and last values are used.
        densitylim_sol : [float, float], optional
            Limits on the sol y axis when ``plot_sol_separately`` is True.
        temperaturelim : [float, float, float] or [float, float], optional
            Limits on the second y axis where the middle value is when the scale
            changes from logarithmic to linear [eV]. If ``plot_sol_separately``
            is True, only the first and last values are used.
        temperaturelim_sol : [float, float], optional
            Limits on the sol y axis when ``plot_sol_separately`` is True.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
            If ``plot_sol_separately`` is True, there should be two axes:
            axes[0] for core and axes[1] for SOL.
        plot_sol_separately : bool, optional
            If True, plot core on left-hand side and SOL on the right-hand side.
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

        rho = pls["rho"]
        if rholim is None: rholim = [pls['rho'][0][0], pls['rho'][-1][0]]

        if plot_sol_separately:
            if rholim_sol is None:
                rholim_sol = [np.maximum(pls['rho'][0][0],0.95), pls['rho'][-1][0]]
            if axes is None:
                fig, axes = a5plt.plt.subplots(1, 2, figsize=(14, 5))
                fig.suptitle(self.get_desc())
            if densitylim is None:
                densitylim = [np.minimum(0, np.amin(pls["idensity"])),
                              1.1*np.amax(pls["edensity"])]
            else:
                # In case three elements are given accidentally
                densitylim = [densitylim[0], densitylim[-1]]
            if densitylim_sol is None:
                densitylim_sol = [np.amin(pls["idensity"]),
                                  1.1*np.amax(pls["edensity"])]
            else:
                densitylim_sol = [densitylim_sol[0], densitylim_sol[-1]]
            if temperaturelim is None:
                temperaturelim = [np.minimum(0, np.amin(pls["itemperature"])),
                                  1.1*np.amax(pls["itemperature"])]
            else:
                temperaturelim = [temperaturelim[0], temperaturelim[-1]]
            if temperaturelim_sol is None:
                temperaturelim_sol = [np.amin(pls["itemperature"]),
                                      1.1*np.amax(pls["itemperature"])]
            else:
                temperaturelim_sol = [temperaturelim_sol[0], temperaturelim_sol[-1]]

            a5plt.radialprofile_single_scale_axes(
                rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
                xlim=rholim,
                y1lim=densitylim, y2lim=temperaturelim,
                xlabel=r"$\rho$", y1label=r"Density [m$^{-3}$]",
                y1legends=None, y2label=None,
                y2legends=None, axes=axes[0], yscale="linear")

            a5plt.radialprofile_single_scale_axes(
                rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
                xlim=rholim_sol,
                y1lim=densitylim_sol, y2lim=temperaturelim_sol,
                xlabel=r"$\rho$", y1label=None,
                y1legends=y1legends, y2label=r"Temperature [eV]",
                y2legends=[r"$T_e$", r"$T_i$"], axes=axes[1], yscale="log")
            a5plt.plt.show()
        else:
            if densitylim is None:
                densitylim = [np.amin(pls["idensity"]),
                              1e19,
                              1.1*np.amax(pls["edensity"])]
            if temperaturelim is None:
                temperaturelim = [np.amin(pls["itemperature"]),
                                  1e3,
                                  1.1*np.amax(pls["itemperature"])]

            a5plt.radialprofile(
                rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
                xlim=rholim, y1lim=densitylim, y2lim=temperaturelim,
                xlabel=r"$\rho$", y1label=r"Density [m$^{-3}$]",
                y1legends=y1legends, y2label=r"Temperature [eV]",
                y2legends=[r"$T_e$", r"$T_i$"], axes=axes, title=self.get_desc())

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

class plasma_2D(DataGroup):
    """Plasma profiles with (R,z) dependency.
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

        for name in ["vtor", "edensity", "etemperature", "idensity",
                     "itemperature"]:
            out[name] = np.transpose(out[name])
        return out

    @staticmethod
    def write_hdf5(fn, nr, nz, nion, rmin, rmax, zmin, zmax, anum, znum, mass,
                   charge, vtor, edensity, etemperature, idensity, itemperature,
                   desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Path to hdf5 file.
        nr : int
            Number of R grid points.
        nz : int
            Number of z grid points.
        nion : int
            Number of ion species.
        rmin : float
            Minimum R grid value.
        rmax : float
            Maximum R grid value.
        zmin : float
            Minimum z grid value.
        zmax : float
            Maximum z grid value.
        anum : array_like (nion,1)
            Ion species atomic mass number
        znum : array_like (nion,1)
            Ion species charge number.
        mass : array_like (nion,1)
            Ion species mass [amu].
        charge : array_like (nion,1)
            Ion species charge [e].
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
        if vtor.shape != (nr,nz):
            raise ValueError("Invalid size for toroidal rotation.")
        if etemperature.shape != (nr,nz):
            raise ValueError("Invalid size for electron temperature.")
        if itemperature.shape != (nr,nz):
            raise ValueError("Invalid size for ion temperature.")
        if edensity.shape != (nr,nz):
            raise ValueError("Invalid size for electron density.")
        if idensity.shape != (nr,nz,nion):
            raise ValueError("Invalid size for ion density.")

        vtor = np.transpose(vtor)
        etemperature = np.transpose(etemperature)
        itemperature = np.transpose(itemperature)
        edensity = np.transpose(edensity)
        idensity = np.transpose(idensity)

        parent = "plasma"
        group  = "plasma_2D"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nion',   (1,1),    data=nion,   dtype='i4')
            g.create_dataset('nr',     (1,1),    data=nr,     dtype='i4')
            g.create_dataset('nz',     (1,1),    data=nz,     dtype='i4')
            g.create_dataset('rmin',   (1,1),    data=rmin,   dtype='f8')
            g.create_dataset('rmax',   (1,1),    data=rmax,   dtype='f8')
            g.create_dataset('zmin',   (1,1),    data=zmin,   dtype='f8')
            g.create_dataset('zmax',   (1,1),    data=zmax,   dtype='f8')
            g.create_dataset('znum',   (nion,1), data=znum,   dtype='i4')
            g.create_dataset('anum',   (nion,1), data=anum,   dtype='i4')
            g.create_dataset('charge', (nion,1), data=charge, dtype='i4')
            g.create_dataset('mass',   (nion,1), data=mass,   dtype='f8')
            g.create_dataset('vtor',   (nz,nr),  data=vtor,   dtype='f8')

            g.create_dataset('etemperature', (nz,nr),    data=etemperature,
                             dtype='f8')
            g.create_dataset('edensity',     (nz,nr),    data=edensity,
                             dtype='f8')
            g.create_dataset('itemperature', (nz,nr),    data=itemperature,
                             dtype='f8')
            g.create_dataset('idensity',     (nion,nz,nr), data=idensity,
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
        nr, nz = 3, 4
        return {
            "nr":nr, "nz":nz, "nion":1, "rmin":0, "rmax":1, "zmin":-1, "zmax":1,
            "znum":np.array([1]), "anum":np.array([1]),
            "mass":np.array([1]), "charge":np.array([1]),
            "vtor":np.zeros((nr,nz)),
            "edensity":np.full((nr,nz), 1e20),
            "etemperature":np.full((nr,nz), 1e3),
            "idensity":np.full((nr,nz,1), 1e20),
            "itemperature":np.full((nr,nz), 1e3)}

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

    def plot_radial(self,
                    rholim=None,
                    rholim_sol=None,
                    densitylim=None,
                    densitylim_sol=None,
                    temperaturelim=None,
                    temperaturelim_sol=None,
                    axes=None,
                    plot_sol_separately=False,
                    ):
        """Plot plasma profiles as a function of rho. If
        ``plot_sol_separately`` is True, then the core and scrape-off layer are
        plotted separately.

        Parameters
        ----------
        rholim : [float, float], optional
            Limits on x-axis. If plot_sol_separately is True, then the limits
            are applied to the left side of x-axis.
        rholim_sol : [float, float], optional
            Limits on the right side of x-axis when ``plot_sol_separately`` is True.
        densitylim : [float, float, float] or [float, float], optional
            Limits on the first y axis where the middle value is when the scale
            changes from logarithmic to linear [m^-3]. If ``plot_sol_separately``
            is True, only the first and last values are used.
        densitylim_sol : [float, float], optional
            Limits on the sol y axis when ``plot_sol_separately`` is True.
        temperaturelim : [float, float, float] or [float, float], optional
            Limits on the second y axis where the middle value is when the scale
            changes from logarithmic to linear [eV]. If ``plot_sol_separately``
            is True, only the first and last values are used.
        temperaturelim_sol : [float, float], optional
            Limits on the sol y axis when ``plot_sol_separately`` is True.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
            If ``plot_sol_separately`` is True, there should be two axes:
            axes[0] for core and axes[1] for SOL.
        plot_sol_separately : bool, optional
            If True, plot core on left-hand side and SOL on the right-hand side.
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

        rho = np.linspace(pls['rhomin'][0], pls['rhomax'][0], pls['nrho'])
        if rholim is None: rholim = [pls['rhomin'][0][0], pls['rhomax'][0][0]]
        if plot_sol_separately:
            if rholim_sol is None:
                rholim_sol = [np.maximum(pls['rhomin'][0][0],0.95), pls['rhomax'][0][0]]

            if axes is None:
                fig, axes = a5plt.plt.subplots(1, 2, figsize=(14, 5))
                fig.suptitle(self.get_desc())
            if densitylim is None:
                densitylim = [np.minimum(0, np.amin(pls["idensity"])),
                              1.1*np.amax(pls["edensity"])]
            else:
                # In case three elements are given accidentally
                densitylim = [densitylim[0], densitylim[-1]]
            if densitylim_sol is None:
                densitylim_sol = [np.amin(pls["idensity"]),
                                  1.1*np.amax(pls["edensity"])]
            else:
                densitylim_sol = [densitylim_sol[0], densitylim_sol[-1]]
            if temperaturelim is None:
                temperaturelim = [np.minimum(0, np.amin(pls["itemperature"])),
                                  1.1*np.amax(pls["itemperature"])]
            else:
                temperaturelim = [temperaturelim[0], temperaturelim[-1]]
            if temperaturelim_sol is None:
                temperaturelim_sol = [np.amin(pls["itemperature"]),
                                      1.1*np.amax(pls["itemperature"])]
            else:
                temperaturelim_sol = [temperaturelim_sol[0], temperaturelim_sol[-1]]

            a5plt.radialprofile_single_scale_axes(
                rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
                xlim=rholim,
                y1lim=densitylim, y2lim=temperaturelim,
                xlabel=r"$\rho$", y1label=r"Density [m$^{-3}$]",
                y1legends=None, y2label=None,
                y2legends=None, axes=axes[0], yscale="linear")

            a5plt.radialprofile_single_scale_axes(
                rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
                xlim=rholim_sol,
                y1lim=densitylim_sol, y2lim=temperaturelim_sol,
                xlabel=r"$\rho$", y1label=None,
                y1legends=y1legends, y2label=r"Temperature [eV]",
                y2legends=[r"$T_e$", r"$T_i$"], axes=axes[1], yscale="log")
            a5plt.plt.show()
        else:
            if densitylim is None:
                densitylim = [np.amin(pls["idensity"]),
                              1e19,
                              1.1*np.amax(pls["edensity"])]
            if temperaturelim is None:
                temperaturelim = [np.amin(pls["itemperature"]),
                                  1e3,
                                  1.1*np.amax(pls["itemperature"])]

            a5plt.radialprofile(
                rho, ndens, y2=[pls['etemperature'], pls['itemperature']],
                xlim=rholim, y1lim=densitylim, y2lim=temperaturelim,
                xlabel=r"$\rho$", y1label=r"Density [m$^{-3}$]",
                y1legends=y1legends, y2label=r"Temperature [eV]",
                y2legends=[r"$T_e$", r"$T_i$"], axes=axes, title=self.get_desc())

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
