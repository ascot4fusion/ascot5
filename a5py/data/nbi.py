"""Defines :class:`NbiBundle` neutral beam injection input class, the actual
injector class :class:`Nbi`, and the corresponding factory method.
"""
import ctypes
from typing import Tuple, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin

class NbiStruct(ctypes.Structure):

    _fields_ = [
        ('id', ctypes.c_int32),
        ('dummy', ctypes.c_int32),
        ('inj', ctypes.c_void_p),
        ]


class NbiBundle(InputVariant):
    """A bundle of neutral beam injectors."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in nbi.h."""

        _fields_ = [
            ('id', ctypes.c_int32),
            ('n_beamlet', ctypes.c_int32),
            ('beamlet_x', ctypes.POINTER(ctypes.c_double)),
            ('beamlet_y', ctypes.POINTER(ctypes.c_double)),
            ('beamlet_z', ctypes.POINTER(ctypes.c_double)),
            ('beamlet_dx', ctypes.POINTER(ctypes.c_double)),
            ('beamlet_dy', ctypes.POINTER(ctypes.c_double)),
            ('beamlet_dz', ctypes.POINTER(ctypes.c_double)),
            ('power', ctypes.c_double),
            ('energy', ctypes.c_double),
            ('efrac', ctypes.c_double * 3),
            ('div_h', ctypes.c_double),
            ('div_v', ctypes.c_double),
            ('div_halo_frac', ctypes.c_double),
            ('div_halo_h', ctypes.c_double),
            ('div_halo_v', ctypes.c_double),
            ('anum', ctypes.c_int32),
            ('znum', ctypes.c_int32),
            ('mass', ctypes.c_double),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="NbiBundle",
            struct=NbiBundle.Struct(),
            )
        self._bxyz: unyt.unyt_array
        self._axisrz: unyt.unyt_array
        self._psival: unyt.unyt_array
        self._rhoval: unyt.unyt_array
        self._jacobian: unyt.unyt_array

    @property
    def bxyz(self) -> unyt.unyt_array:
        """Magnetic field vector in Cartesian basis at origo."""
        if self._format == Format.HDF5:
            return self._read_hdf5("bxyz")
        if self._format == Format.CSTRUCT:
            return self._from_struct_("B", shape=(3,), units="T")
        return self._bxyz.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        """ginj.create_dataset("ids",  (1,), data=inj.ids,  dtype="i4")
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
            "beamletdz", (nbeamlet,), data=inj.beamletdz, dtype="f8")"""
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
        }
        return data

    def stage(self):
        init = LIBASCOT.NBI_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ndpointer(ctypes.c_double),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.exyz,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._exyz
            self._staged = True

    def unstage(self):
        free = LIBASCOT.NBI_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._exyz = self.exyz
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateNbiMixin(TreeMixin):
    """Mixin class used by :class:`Data` to create :class:`NbiBundle` input."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_nbibundle(
            self,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> NbiBundle:
        r"""Create one (or a bundle of) neutral beam injector(s).

        This input is used by BBNBI to launch neutral particles to the plasma
        and to generate NBI-ion source.

        Markers generated from this input are not separable, i.e. they
        are all treated as a single source. If you wish to have separate markers
        for different injectors, you should create separate NBI inputs for
        those.

        Parameters
        ----------
        injectors : array of :class:`Injector`
            Array of injectors used in this bundle.
        note : str, optional
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, optional
            Set this input as active on creation.
        dryrun : bool, optional
            Do not add this input to the `data` structure or store it on disk.

            Use this flag to modify the input manually before storing it.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        inputdata : NbiBundle
            Freshly minted input data object.

        Notes
        -----
        A single injector consists of a number of beamlets that have common
        weights (each beamlet produces the same number of particles) and
        common properties such as energy, divergence, etc. When a neutral marker
        is generated at the beamlet position, its velocity is diverted from the
        beamlet direction vector randomly, with the divergence angle given by

        .. math::
        """
        parameters = _variants.parse_parameters(
            bxyz, jacobian, axisrz, rhoval, psival,
        )
        _variants.validate_required_parameters(
            parameters,
            names=["bxyz", "jacobian", "axisrz", "rhoval"],
            units=["T", "T/m", "m", "1"],
            shape=[(3,), (3,3), (2,), ()],
            dtype="f8",
            default=[
                np.ones((3,)), np.zeros((3,3)), (1.0, 0.0), 0.5,
                ],
        )
        _variants.validate_optional_parameters(
            parameters,
            ["psival"], ["Wb/m"], [()], "f8", [parameters["rhoval"]],
        )
        meta = _variants.new_metadata("B_TC", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj


class Nbi():
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
        Mass of the injected species [amu].
    energy : float
        Full injection energy [eV].
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
        self._beamletx, self._beamlety, self._beamletz = (
            beamletx, beamlety, beamletz
            )
        self._beamletdx, self._beamletdy, self._beamletdz = (
            beamletdx, beamletdy, beamletdz
            )

    @property
    def beamletxyz(self) -> np.ndarray:
        """Beamlet locations in Cartesian coordinates.

        Returns
        -------
        beamletxyz : np.ndarray, shape (nbeamlet, 3)
            Beamlet locations in Cartesian coordinates.
        """
        return np.column_stack((self.beamletx, self.beamlety, self.beamletz))

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

    @staticmethod
    def write_hdf5(fn, ninj, injectors, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        
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

                

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is a very large rectangular wall.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        return {"id":1, "nbeamlet":1,
                "beamletx":np.array([1.0]), "beamlety":np.array([0.0]),
                "beamletz":np.array([0.0]),
                "beamletdx":np.array([1.0]), "beamletdy":np.array([0.0]),
                "beamletdz":np.array([0.0]),
                "div_h":0.0, "div_v":0.0, "div_halo_frac":0.0,
                "div_halo_h":0.0, "div_halo_v":0.0,
                "anum":1.0, "znum":1.0, "mass":1.0, "energy":1.0,
                "efrac":[1,0,0], "power":1}
