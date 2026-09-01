"""Marker initial and final phase-space positions and related quantities.
"""
import ctypes
import numpy as np
import unyt
import a5py.physlib as physlib
from a5py.libascot import DataStruct, LIBASCOT
from a5py.physlib.formulas import resolve_quantities, compute_quantities
from a5py.engine.interpolate import interpolate


class Structure(DataStruct):
    """Python wrapper for the particle state struct in particle.h."""

    _fields_ = [
        ("r", ctypes.c_double),
        ("z", ctypes.c_double),
        ("phi", ctypes.c_double),
        ("rprt", ctypes.c_double),
        ("zprt", ctypes.c_double),
        ("phiprt", ctypes.c_double),
        ("zeta", ctypes.c_double),
        ("ekin", ctypes.c_double),
        ("pitch", ctypes.c_double),
        ("pr", ctypes.c_double),
        ("pz", ctypes.c_double),
        ("pphi", ctypes.c_double),
        ("mass", ctypes.c_double),
        ("charge", ctypes.c_double),
        ("time", ctypes.c_double),
        ("theta", ctypes.c_double),
        ("weight", ctypes.c_double),
        ("mileage", ctypes.c_double),
        ("cputime", ctypes.c_double),
        ("id", ctypes.c_size_t),
        ("walltile", ctypes.c_size_t),
        ("endcond", ctypes.c_uint),
        ("err", ctypes.c_uint),
        ("anum", ctypes.c_int),
        ("znum", ctypes.c_int),
        ]


class MarkerState():
    """State of the markers at the fixed point in simulation workflow."""

    def __init__(self):
        self._file = None
        self._cdata = None

    @property
    def r(self) -> unyt.unyt_array:
        r"""Guiding-center :math:`R` coordinate."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].r
            return out
        return self._file.read("r")

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Guiding-center :math:`\phi` coordinate."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "rad")
            for i in range(out.size):
                out[i] = self._cdata[i].phi
            out.convert_to_units("deg")
            return out
        return self._file.read("phi")

    @property
    def z(self) -> unyt.unyt_array:
        r"""Guiding-center :math:`z` coordinate."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].z
            return out
        return self._file.read("z")

    @property
    def ekin(self) -> unyt.unyt_array:
        r"""Guiding-center kinetic energy."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "J")
            for i in range(out.size):
                out[i] = self._cdata[i].ekin
            return out.to("eV")
        return self._file.read("ekin")

    @property
    def pitch(self) -> unyt.unyt_array:
        r"""Guiding-center pitch."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "1")
            for i in range(out.size):
                out[i] = self._cdata[i].pitch
            return out
        return self._file.read("pitch")

    @property
    def zeta(self) -> unyt.unyt_array:
        r"""Guiding-center gyro-angle."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "rad")
            for i in range(out.size):
                out[i] = self._cdata[i].zeta
            return out
        return self._file.read("zeta")

    @property
    def rprt(self) -> unyt.unyt_array:
        r"""Particle :math:`R` coordinate."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].rprt
            return out
        return self._file.read("rprt")

    @property
    def phiprt(self) -> unyt.unyt_array:
        r"""Particle :math:`\phi` coordinate."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "rad")
            for i in range(out.size):
                out[i] = self._cdata[i].phiprt
            out.convert_to_units("deg")
            return out
        return self._file.read("phiprt")

    @property
    def zprt(self) -> unyt.unyt_array:
        r"""Particle :math:`z` coordinate."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].zprt
            return out
        return self._file.read("zprt")

    @property
    def pr(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`R` component."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pr
            return out
        return self._file.read("pr")

    @property
    def pr(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`R` component."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pr
            return out
        return self._file.read("pr")

    @property
    def pphi(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`\phi` component."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pphi
            return out
        return self._file.read("pphi")

    @property
    def pz(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`z` component."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pz
            return out
        return self._file.read("pz")

    @property
    def mass(self) -> unyt.unyt_array:
        r"""Marker mass."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "kg")
            for i in range(out.size):
                out[i] = self._cdata[i].mass
            return out
        return self._file.read("mass")

    @property
    def charge(self) -> unyt.unyt_array:
        r"""Marker charge."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "C")
            for i in range(out.size):
                out[i] = self._cdata[i].charge
            return out
        return self._file.read("charge")

    @property
    def anum(self) -> unyt.unyt_array:
        r"""Marker atomic mass number."""
        if self._file is None:
            out = np.array([0]*self.n)
            for i in range(out.size):
                out[i] = self._cdata[i].anum
            return out
        return self._file.read("anum")

    @property
    def znum(self) -> unyt.unyt_array:
        r"""Marker charge number."""
        if self._file is None:
            out = np.array([0]*self.n)
            for i in range(out.size):
                out[i] = self._cdata[i].znum
            return out
        return self._file.read("znum")

    @property
    def weight(self) -> unyt.unyt_array:
        r"""Marker weight."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "particles/s")
            for i in range(out.size):
                out[i] = self._cdata[i].weight
            return out
        return self._file.read("weight")

    @property
    def time(self) -> unyt.unyt_array:
        r"""Marker time."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "s")
            for i in range(out.size):
                out[i] = self._cdata[i].time
            return out
        return self._file.read("time")

    @property
    def mileage(self) -> unyt.unyt_array:
        r"""Marker mileage."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "s")
            for i in range(out.size):
                out[i] = self._cdata[i].mileage
            return out
        return self._file.read("mileage")

    @property
    def cputime(self) -> unyt.unyt_array:
        r"""Marker CPU time."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "s")
            for i in range(out.size):
                out[i] = self._cdata[i].cputime
            return out
        return self._file.read("cputime")

    @property
    def theta(self) -> unyt.unyt_array:
        r"""Marker theta."""
        if self._file is None:
            out = unyt.unyt_array([0.]*self.n, "deg")
            for i in range(out.size):
                out[i] = self._cdata[i].theta
            return out
        return self._file.read("theta")

    @property
    def ids(self) -> np.ndarray:
        """Unique identifier for each marker."""
        if self._cdata is not None:
            out = np.zeros((self.n,), dtype="i8")
            for i in range(out.size):
                out[i] = self._cdata[i].id
            return out
        if self._file is not None:
            return self._file.read("ids")

    @property
    def endcond(self) -> np.ndarray:
        """Marker end condition."""
        if self._cdata is not None:
            out = np.zeros((self.n,), dtype="i8")
            for i in range(out.size):
                out[i] = self._cdata[i].endcond
            return out
        if self._file is not None:
            return self._file.read("endcond")

    @property
    def walltile(self) -> np.ndarray:
        """ID of the wall tile if the marker hit one."""
        if self._cdata is not None:
            out = np.zeros((self.n,), dtype="i8")
            for i in range(out.size):
                out[i] = self._cdata[i].walltile
            return out
        if self._file is not None:
            return self._file.read("walltile")

    @property
    def err(self) -> np.ndarray:
        """Error code if marker."""
        if self._cdata is not None:
            out = np.zeros((self.n,), dtype="i8")
            for i in range(out.size):
                out[i] = self._cdata[i].walltile
            return out
        if self._file is not None:
            return self._file.read("walltile")

    @property
    def n(self) -> int:
        r"""Number of markers."""
        if self._file is None:
            return len(self._cdata)
        return self._file.read("ids").size

    def save(self, file):
        """"""
        for field in Structure._fields_:
            name = field[0]
            if name == "id":
                name = "ids"
            file.write(name, getattr(self, name))
        del self._cdata
        self._file = file

    def combine(self, file=None):
        """"""

    def evaluate(self, *qnt, filter=None, bfield=None, mode=None):
        marker_ids = self.ids
        if filter is not None:
            marker_ids = marker_ids[np.isin(marker_ids, filter)]
        sorter = np.argsort(marker_ids)

        marker_quantities = [
            "r", "z", "phi", "time", "charge", "weight", "vr", "vphi",
            "vz", "mass"
        ]
        marker_quantities += [
            "gyroangle", "pitch", "ekin",
        ]
        field_quantities = [
            "psidr", "psidphi", "psidz", "br", "bphi", "bz", "psi", "brdr",
            "brdphi", "brdz", "bphidr", "bphidphi", "bphidz", "bzdr", "bzdphi",
            "bzdz",
        ]

        available = marker_quantities + field_quantities
        fundamentals_needed, all_qnt = resolve_quantities(available, qnt)

        field_quantities = list(set(fundamentals_needed) & set(field_quantities))
        marker_quantities = list(set(fundamentals_needed) - set(field_quantities))

        if field_quantities:
            marker_quantities = list(set(marker_quantities + ["r", "z", "phi", "time"]))

        evaluated = {}
        for q in marker_quantities:
            evaluated[q] = getattr(self, q)[sorter]

        if field_quantities:
            evaluated.update(
                interpolate(
                    evaluated["r"], evaluated["phi"], evaluated["z"], evaluated["time"],
                    0, *field_quantities, bfield=bfield
                    )
            )
        evaluated.update(compute_quantities(evaluated, all_qnt))

        out = []
        for q in qnt:
            out.append(evaluated[q])
        return out

    @classmethod
    def from_params(cls, mrk):
        """"""
        obj = cls()
        cdata = (Structure * len(mrk))()
        for i in range(len(mrk)):
            for field in Structure._fields_:
                val = getattr(mrk[i], field[0])
                setattr(cdata[i], field[0], val)
        obj._cdata = cdata
        return obj

    @classmethod
    def from_markerinput(cls, mrk, bjacb, idx=None):
        """Initialize a new state that matches the marker input.

        Parameters
        ----------
        mrk
            The marker input.
        bjacb
            Magnetic field vector and gradient at each marker position.
        idx
            Index of the markers to include.
        """
        if idx is None:
            idx = np.arange(mrk.n)

        obj = cls()
        obj._cdata = (Structure * idx.size)()

        isfieldline = mrk._cdata[0].mass == 0
        isparticle = not isfieldline and mrk._cdata[0].ekin == 0
        isguidingcenter = not isfieldline and not isparticle
        if isfieldline:
            zeros = [
                "pr", "pphi", "pz", "mass", "charge", "anum", "znum", "weight",
                "theta", "endcond", "walltile", "cputime", "mileage", "ekin",
                "zeta",
                ]
            for i, j in enumerate(idx):
                ctypes.memmove(
                    ctypes.addressof(obj._cdata[i]),
                    ctypes.addressof(mrk._cdata[j]),
                    ctypes.sizeof(Structure),
                    )
                obj._cdata[i].rprt = obj._cdata[i].r
                obj._cdata[i].phiprt = obj._cdata[i].phi
                obj._cdata[i].zprt = obj._cdata[i].z
                for zero in zeros:
                    setattr(obj._cdata[i], zero, 0)

        elif isparticle:
            zeros = ["theta", "endcond", "walltile", "cputime", "mileage"]
            r, phi, z, ppar, mu, zeta = (
                np.array([0.0]), np.array([0.0]), np.array([0.0]),
                np.array([0.0]), np.array([0.0]), np.array([0.0]),
            )
            for i, j in enumerate(idx):
                ctypes.memmove(
                    ctypes.addressof(obj._cdata[i]),
                    ctypes.addressof(mrk._cdata[j]),
                    ctypes.sizeof(Structure),
                    )
                LIBASCOT.gctransform_particle2guidingcenter(
                    obj._cdata[i].mass,
                    obj._cdata[i].charge,
                    bjacb,
                    obj._cdata[i].rprt,
                    obj._cdata[i].phiprt,
                    obj._cdata[i].zprt,
                    obj._cdata[i].pr,
                    obj._cdata[i].pphi,
                    obj._cdata[i].pz,
                    r,
                    phi,
                    z,
                    ppar,
                    mu,
                    zeta
                )
                obj._cdata[i].ekin = 0.0
                obj._cdata[i].pitch = 0.0
                for zero in zeros:
                    setattr(obj._cdata[i], zero, 0)

        elif isguidingcenter:
            from a5py.physlib import Quantity
            r, phi, z, pr, pphi, pz, pparprt, muprt, zetaprt = (
                np.array([0.0]), np.array([0.0]), np.array([0.0]),
                np.array([0.0]), np.array([0.0]), np.array([0.0]),
                np.array([0.0]), np.array([0.0]), np.array([0.0]),
            )
            zeros = ["theta", "endcond", "walltile", "cputime", "mileage"]
            for i, j in enumerate(idx):
                ctypes.memmove(
                    ctypes.addressof(obj._cdata[i]),
                    ctypes.addressof(mrk._cdata[j]),
                    ctypes.sizeof(Structure),
                    )
                pnorm = Quantity.registry["pnorm"].compute(
                    mass=obj._cdata[i].mass*unyt.kg,
                    ekin=obj._cdata[i].ekin*unyt.J,
                    ).to("kg*m/s")
                ppar = Quantity.registry["ppar"].compute(
                    pnorm=pnorm,
                    pitch=obj._cdata[i].pitch,
                    )
                mu = Quantity.registry["mu"].compute(
                    mass=obj._cdata[i].mass*unyt.kg,
                    pitch=obj._cdata[i].pitch,
                    pnorm=pnorm,
                    bnorm=np.sqrt(bjacb[i,0]**2+bjacb[i,1]**2+bjacb[i,2]**2)*unyt.T,
                    ).to("J/T")
                bjacb0 = np.zeros((12,))
                bjacb0[:] = bjacb[i,:]
                LIBASCOT.gctransform_guidingcenter2particle(
                    obj._cdata[i].mass,
                    obj._cdata[i].charge,
                    bjacb0,
                    obj._cdata[i].r,
                    obj._cdata[i].phi,
                    obj._cdata[i].z,
                    ppar,
                    mu,
                    obj._cdata[i].zeta,
                    r,
                    phi,
                    z,
                    pparprt,
                    muprt,
                    zetaprt,
                )
                LIBASCOT.gctransform_pparmuzeta2prpphipz(
                    obj._cdata[i].mass,
                    obj._cdata[i].charge,
                    bjacb0,
                    phi[0],
                    pparprt[0],
                    muprt[0],
                    zetaprt[0],
                    pr,
                    pphi,
                    pz,
                )
                obj._cdata[i].rprt = r[0]
                obj._cdata[i].phiprt = phi[0]
                obj._cdata[i].zprt = z[0]
                obj._cdata[i].pr = pr[0]
                obj._cdata[i].pphi = pphi[0]
                obj._cdata[i].pz = pz[0]
                for zero in zeros:
                    setattr(obj._cdata[i], zero, 0)

        return obj
