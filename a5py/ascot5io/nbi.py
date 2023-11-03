"""Input data representing NBI.

Neutral beam injectors are used by BBNBI to generate NBI-ion source and
calculate shinethrough.
"""
import numpy as np
import h5py
import unyt

from a5py.routines.plotting import openfigureifnoaxes

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class Injector():
    """A single injector made up of beamlets

    Attributes
    ----------
    ids : int
        Numerical identifier for this injector which should be unique in
        an input bundle.
    anum : int
        Mass number of injected species.
    znum : int
        Nuclear charge number of injected species.
    mass : float
        Mass of the injected species [kg].
    energy : float
        Full injection energy [J].
    efrac : array_like (3,)
        Particle fractions for full, 1/2 and 1/3 energies.
    power : float
        Injected power [W].
    divh : float
        Horizontal divergence [rad].
    divv : float
        Vertical divergence [rad].
    divhalofrac : float
        Fraction of particles with halo divergence.
    divhaloh : float
        Horizontal divergence in halo [rad].
    divhalov : float
        Vertical divergence in halo [rad].
    nbeamlet : int
        Number of beamlets.
    beamletx : array_like, (nbeamlet,)
        x coordinates of beamlets [m].
    beamlety : array_like, (nbeamlet,)
        y coordinates of beamlets [m].
    beamletz : array_like, (nbeamlet,)
        z coordinates of beamlets [m].
    beamletdx : array_like, (nbeamlet,)
        x components of the unit direction vector of beamlets.
    beamletdy : array_like, (nbeamlet,)
        y components of the unit direction vector of beamlets.
    beamletdz : array_like, (nbeamlet,)
        z components of the unit direction vector of beamlets.
    """

    def __init__(self, ids, anum, znum, mass, energy, efrac, power, divh, divv,
                 divhalofrac, divhaloh, divhalov, nbeamlet, beamletx, beamlety,
                 beamletz, beamletdx, beamletdy, beamletdz):
        """Create a single injector.

        Parameters
        ----------
        ids : int
            Numerical identifier for this injector which should be unique in
            an input bundle.
        anum : int
            Mass number of injected species.
        znum : int
            Nuclear charge number of injected species.
        mass : float
            Mass of the injected species [kg].
        energy : float
            Full injection energy [J].
        efrac : array_like (3,)
            Particle fractions for full, 1/2 and 1/3 energies.
        power : float
            Injected power [W].
        divh : float
            Horizontal divergence [rad].
        divv : float
            Vertical divergence [rad].
        divhalofrac : float
            Fraction of particles with halo divergence.
        divhaloh : float
            Horizontal divergence in halo [rad].
        divhalov : float
            Vertical divergence in halo [rad].
        nbeamlet : int
            Number of beamlets.
        beamletx : array_like, (nbeamlet,)
            x coordinates of beamlets [m].
        beamlety : array_like, (nbeamlet,)
            y coordinates of beamlets [m].
        beamletz : array_like, (nbeamlet,)
            z coordinates of beamlets [m].
        beamletdx : array_like, (nbeamlet,)
            x components of the unit direction vector of beamlets.
        beamletdy : array_like, (nbeamlet,)
            y components of the unit direction vector of beamlets.
        beamletdz : array_like, (nbeamlet,)
            z components of the unit direction vector of beamlets.
        """
        if efrac.shape != (3,):
            raise ValueError("efrac must have shape (3,) but was "
                             + str(efrac.shape))
        for v in [beamletx, beamlety, beamletz, beamletdx, beamletdy,
                  beamletdz]:
            if v.shape != (nbeamlet,):
                raise ValueError(
                    "Beamlet data is inconsistent or has a wrong shape: "
                    + "expected (%d,)" % nbeamlet)

        self.ids         = ids
        self.anum        = anum
        self.znum        = znum
        self.mass        = mass
        self.energy      = energy
        self.efrac       = efrac
        self.power       = power
        self.divh        = divh
        self.divv        = divv
        self.divhalofrac = divhalofrac
        self.divhaloh    = divhaloh
        self.divhalov    = divhalov
        self.nbeamlet    = nbeamlet
        self.beamletx    = beamletx
        self.beamlety    = beamlety
        self.beamletz    = beamletz
        self.beamletdx   = beamletdx
        self.beamletdy   = beamletdy
        self.beamletdz   = beamletdz

    def plotbeamlets(self, view="3d", direction=True, axes=None):
        """Plot the beamlets of this injector.

        Parameters
        ----------
        view : {"rz", "xy", "3d"}, optional
            Control whether the plot is shown in (R,z) or (x,y) plane or in 3D.
        direction : bool, optional
            Flag for showing the beam direction with arrows originating from
            beamlets.

            If False, only the beamlet locations are shown.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        """
        if view == "3d":
            x = self.beamletx
            y = self.beamlety
            z = self.beamletz
            dx = self.beamletdx
            dy = self.beamletdy
            dz = self.beamletdz

            @openfigureifnoaxes("3d")
            def plot(axes=None):
                axes.scatter(x, y, z)
                axes.set_xlabel("x [m]")
                axes.set_ylabel("y [m]")
                axes.set_zlabel("z [m]")
                if direction:
                    axes.quiver3D(x, y, z, dx, dy, dz)
            plot(axes=axes)
        else:
            if view == "xy":
                x = self.beamletx
                y = self.beamlety
                dx = self.beamletdx
                dy = self.beamletdy
                xlabel = "x [m]"
                ylabel = "y [m]"
            elif view == "rz":
                x = np.sqrt( self.beamletx**2 + self.beamlety**2 )
                y = self.beamletz
                dx = np.sqrt( ( self.beamletx + self.beamletdx )**2
                            + ( self.beamlety + self.beamletdy )**2 ) - x
                dy = self.beamletdz
                xlabel = "R [m]"
                ylabel = "z [m]"
            else:
                raise ValueError(
                    "Argument 'view' must be one of the following:"\
                    + " 'rz', 'xy', '3d'")

            @openfigureifnoaxes()
            def plot(axes=None):
                axes.set_xlabel(xlabel)
                axes.set_ylabel(ylabel)
                axes.scatter(x, y)
                if direction:
                    axes.quiver(x, y, dx, dy)
            plot(axes=axes)

class NBI(DataGroup):
    """Object representing a bundle of injectors used as an input in BBNBI.

    An injector bundle means that the particles generated from those with BBNBI
    are not separable. If your work requires that the different NBI sources are
    kept separate, then use different bundles for those (a bundle can consists
    of a single injector). Otherwise it is more convenient to combine injectors
    in a single bundle since it then produces a single BBNBI output.
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
        with h5py.File(fn, "r") as f:
            ninj = int(f[path]["ninj"][0])
            out["injectors"] = [None] * ninj
            for i in range(ninj):
                kwargs = {}
                inj = f[path + "/inj" + str(i+1)]
                for key in inj:
                    kwargs[key] = inj[key][:]
                out["injectors"][i] = Injector(**kwargs)

        out["ninj"] = ninj
        return out

    @staticmethod
    def write_hdf5(fn, ninj, injectors, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        ninj : int
            Number of injectors.
        injectors : array of :class:`Injector`
            Array of injectors used in this bundle.
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        parent = "nbi"
        group  = "nbi"
        gname  = ""
        if ninj != len(injectors):
            raise ValueError("Number of injectors is not the same as len(nbi)")

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("ninj",  (1,), data=ninj,    dtype="i4")

            for i in range(ninj):
                inj      = injectors[i]
                nbeamlet = inj.nbeamlet
                ginj     = g.create_group("inj"+str(i+1))

                ginj.create_dataset("ids",  (1,), data=inj.ids,  dtype="i4")
                ginj.create_dataset("anum", (1,), data=inj.anum, dtype="i4")
                ginj.create_dataset("znum", (1,), data=inj.znum, dtype="i4")
                ginj.create_dataset("mass", (1,), data=inj.mass, dtype="f8")
                ginj.create_dataset("divh", (1,), data=inj.divh, dtype="f8")
                ginj.create_dataset("divv", (1,), data=inj.divv, dtype="f8")
                ginj.create_dataset(
                    "divhalofrac", (1,), data=inj.divhalofrac, dtype="f8")
                ginj.create_dataset(
                    "divhaloh",    (1,), data=inj.divhaloh,    dtype="f8")
                ginj.create_dataset(
                    "divhalov",    (1,), data=inj.divhalov,    dtype="f8")
                ginj.create_dataset(
                    "energy",      (1,), data=inj.energy,      dtype="f8")
                ginj.create_dataset(
                    "efrac",       (3,), data=inj.efrac,       dtype="f8")
                ginj.create_dataset(
                    "power",       (1,), data=inj.power,       dtype="f8")
                ginj.create_dataset(
                    "nbeamlet",  (1,),        data=inj.nbeamlet,  dtype="i4")
                ginj.create_dataset(
                    "beamletx",  (nbeamlet,), data=inj.beamletx,  dtype="f8")
                ginj.create_dataset(
                    "beamlety",  (nbeamlet,), data=inj.beamlety,  dtype="f8")
                ginj.create_dataset(
                    "beamletz",  (nbeamlet,), data=inj.beamletz,  dtype="f8")
                ginj.create_dataset(
                    "beamletdx", (nbeamlet,), data=inj.beamletdx, dtype="f8")
                ginj.create_dataset(
                    "beamletdy", (nbeamlet,), data=inj.beamletdy, dtype="f8")
                ginj.create_dataset(
                    "beamletdz", (nbeamlet,), data=inj.beamletdz, dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is a very large rectangular wall.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        nbi = {
            "id" : 1, "nbeamlet" : 1,
            "beamletx"  : np.array([1.0]), "beamlety"  : np.array([0.0]),
            "beamletz"  : np.array([0.0]),
            "beamletdx" : np.array([1.0]), "beamletdy" : np.array([0.0]),
            "beamletdz" : np.array([0.0]),
            "div_h" : 0.0, "div_v" : 0.0, "div_halo_frac" : 0.0,
            "div_halo_h" : 0.0, "div_halo_v" : 0.0,
            "anum" : 1.0,  "znum" : 1.0, "mass" : 1.0, "energy" : 1.0,
            "efrac" : [1,0,0], "power" : 1
        }
        return NBI.write_hdf5(fn, 1, [nbi], desc="DUMMY")

    @staticmethod
    def generate(r, phi, z, tanrad, focallen, dgrid, nbeamlet,
                 anum, znum, mass, energy, efrac, power, div, tilt=None,
                 tanz=None, divhalofrac=None, divhalo=None, ids=1, desc=None):
        """Generate a generic injector.

        The injector grid can be either circular or rectangular, and will be
        randomly populated with beamlets with a common focus point.

        The contents of this function can be used as a template for implementing
        machine-specific injectors with an actual geometry.

        Parameters
        ----------
        r : float
            Injector center point R-coordinate [m].
        phi : float
            Injector center point phi-coordinate [deg].
        z : float
            Injector center point z-coordinate [m].
        tanrad : float
            Tangency radius of the injector centerline [m].
        focallen : float
            Distance to the focus point from the center of the grid
        dgrid : float or array_like (2,)
            Diameter of circular grid or width and height of rectangular grid
        nbeamlet : int
            Number of beamlets in this injector
        anum : int
            Mass number of injected species
        znum : int
            Nuclear charge number of injected species
        mass : float
            Mass of the injected species [amu]
        energy : float
            Full injection energy [keV]
        efrac : array_like (3,)
            Particle fractions for full, 1/2 and 1/3 energies
        power : float
            Injected power [W]
        div : array_like (2,)
            Horizontal and vertical divergences [rad]
        tilt : float, optional
            Vertical tilt angle of the beam centerline [deg]
        tanz : float, optional
            Vertical shift of the tangency radius point [m]
        divhalofrac : float, optional
            Fraction of particles with halo divergence
        divhalo : float, optional
            Horizontal and vertical divergences [rad]
        ids : int, optional
            Unique numerical identifier for this injector.

        Returns
        -------
        inj : :class:`Injector`
            A single injector.
        """
        phi    = phi * np.pi/180.0
        phidir = np.pi + phi - np.arcsin(tanrad/r);
        center = np.array([r * np.cos(phi), r * np.sin(phi), z])

        unitd = np.array([np.cos(phidir), np.sin(phidir), 0.0])
        unitz = np.array([0.0, 0.0, 1.0])
        unith = np.cross(unitd, unitz)

        if tilt is not None:
            tilt = tilt * np.pi/180.0
        elif tanz is not None:
            tilt = tanz / r * np.cos(phidir)
        else:
            tilt = 0

        focus = center + focallen * unitd + focallen * np.tan(tilt) * unitz

        n = {}
        n["ids"]       = 1
        n["nbeamlet"]  = nbeamlet
        n["beamletx"]  = np.zeros(nbeamlet)
        n["beamlety"]  = np.zeros(nbeamlet)
        n["beamletz"]  = np.zeros(nbeamlet)
        n["beamletdx"] = np.zeros(nbeamlet)
        n["beamletdy"] = np.zeros(nbeamlet)
        n["beamletdz"] = np.zeros(nbeamlet)
        n["divh"]      = div[0]
        n["divv"]      = div[1]

        if divhalofrac is not None:
            n["divhalofrac"] = divhalofrac
            n["divhaloh"]    = divhalo[0]
            n["divhalov"]    = divhalo[1]
        else:
            n["divhalofrac"] = 0.0
            n["divhaloh"]    = 1e-10
            n["divhalov"]    = 1e-10

        # Injected species
        n["anum"]   = anum
        n["znum"]   = znum
        n["mass"]   = mass
        n["energy"] = energy
        n["efrac"]  = np.array(efrac)
        n["power"]  = power

        for i in range(nbeamlet):
            if np.array(dgrid).size == 1:
                # circular grid
                d = np.inf
                while d > 0.5 * dgrid:
                    dx = 0.5 * dgrid * ( 2 * np.random.rand(1) - 1 )
                    dy = 0.5 * dgrid * ( 2 * np.random.rand(1) - 1 )
                    d  = np.sqrt(dx**2 + dy**2)
            else:
                # rectangular grid
                dx = 0.5 * dgrid[0] * ( 2 * np.random.rand(1) - 1 )
                dy = 0.5 * dgrid[1] * ( 2 * np.random.rand(1) - 1 )

            beamletxyz = center + dx*unith + dy*unitz
            dirxyz = focus - beamletxyz
            dirxyz = dirxyz / np.linalg.norm(dirxyz)

            n["beamletx"][i]  = beamletxyz[0]
            n["beamlety"][i]  = beamletxyz[1]
            n["beamletz"][i]  = beamletxyz[2]
            n["beamletdx"][i] = dirxyz[0]
            n["beamletdy"][i] = dirxyz[1]
            n["beamletdz"][i] = dirxyz[2]

        return Injector(**n)
