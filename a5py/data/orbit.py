"""Orbit diagnostics.
"""
import numpy as np
import unyt

import ctypes
from typing import Optional

from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.data.options import OrbitParams
from a5py.physlib import formulas
from a5py.libascot import DataStruct
from a5py.engine.interpolate import evaluate


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    _fields_ = [
        ("r", ctypes.POINTER(ctypes.c_double)),
        ("z", ctypes.POINTER(ctypes.c_double)),
        ("phi", ctypes.POINTER(ctypes.c_double)),
        ("p1", ctypes.POINTER(ctypes.c_double)),
        ("p2", ctypes.POINTER(ctypes.c_double)),
        ("p3", ctypes.POINTER(ctypes.c_double)),
        ("mileage", ctypes.POINTER(ctypes.c_double)),
        ("stamp", ctypes.POINTER(ctypes.c_double)),
        ("id", ctypes.POINTER(ctypes.c_size_t)),
        ("idx", ctypes.POINTER(ctypes.c_size_t)),
        ("charge", ctypes.POINTER(ctypes.c_int)),
        ("poincare", ctypes.POINTER(ctypes.c_int)),
        ("simmode", ctypes.POINTER(ctypes.c_int)),
        ("npoint", ctypes.c_size_t),
        ("ntoroidal", ctypes.c_size_t),
        ("npoloidal", ctypes.c_size_t),
        ("nradial", ctypes.c_size_t),
        ("interval", ctypes.c_double),
        ("toroidal", ctypes.POINTER(ctypes.c_double)),
        ("poloidal", ctypes.POINTER(ctypes.c_double)),
        ("radial", ctypes.POINTER(ctypes.c_double)),
    ]


class Orbit():
    """Orbit diagnostics that collect marker phase-space coordinates and related
    quantities at certain points along the marker orbit.
    """

    GYROORBIT     = 1
    GUIDINGCENTER = 2
    FIELDLINE     = 3

    def __init__(self):
        self._file = None
        self._cdata = None

    @property
    def ntotal(self) -> int:
        r"""Total size of the orbit data."""
        if self._file is None:
            return self._cdata.id_ref.size
        return self._file.read("id").size

    @property
    def id(self) -> np.ndarray:
        r"""Unique identifier for each marker."""
        if self._file is None:
            return self._cdata.id_ref
        return self._file.read("id")

    @property
    def r(self) -> unyt.unyt_array:
        r"""Marker :math:`R` coordinate."""
        if self._file is None:
            return self._cdata.r_ref * unyt.m
        return self._file.read("r")

    @property
    def z(self) -> unyt.unyt_array:
        r"""Marker :math:`z` coordinate."""
        if self._file is None:
            return self._cdata.z_ref * unyt.m
        return self._file.read("z")

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Marker :math:`\phi` coordinate."""
        if self._file is None:
            return self._cdata.phi_ref * unyt.rad
        return self._file.read("phi")

    @property
    def p1(self) -> unyt.unyt_array:
        r"""Marker first momentum coordinate."""
        if self._file is None:
            return self._cdata.p1_ref
        return self._file.read("p1")

    @property
    def p2(self) -> unyt.unyt_array:
        r"""Marker second momentum coordinate."""
        if self._file is None:
            return self._cdata.p2_ref
        return self._file.read("p2")

    @property
    def p3(self) -> unyt.unyt_array:
        r"""Marker third momentum coordinate."""
        if self._file is None:
            return self._cdata.p3_ref
        return self._file.read("p3")

    @property
    def charge(self) -> np.ndarray:
        r"""Marker charge."""
        if self._file is None:
            return self._cdata.charge_ref * unyt.e
        return self._file.read("charge")

    @property
    def mileage(self) -> np.ndarray:
        r"""Marker mileage."""
        if self._file is None:
            return self._cdata.mileage_ref
        return self._file.read("mileage")

    @property
    def poincare(self) -> np.ndarray:
        r"""Poincare classification."""
        if self._file is None:
            return self._cdata.poincare_ref
        return self._file.read("poincare")

    @property
    def simmode(self) -> np.ndarray:
        r"""Simulation mode."""
        if self._file is None:
            return self._cdata.simmode_ref
        return self._file.read("simmode")

    @classmethod
    def from_params(cls, nmrk: int, params: OrbitParams):
        """Initialize orbit diagnostics from simulation parameters."""
        if params.collect == "no":
            return None

        cdata = Struct()
        for field in ["r", "z", "phi", "p1", "p2", "p3", "mileage", "stamp",
                      "id", "idx", "charge", "poincare", "simmode"]:
            dtype, ctype = np.float64, ctypes.c_double
            if field in ["id", "idx"]:
                dtype, ctype = np.int64, ctypes.c_size_t
            elif field in ["charge", "poincare", "simmode"]:
                dtype, ctype = np.int32, ctypes.c_int
            size = nmrk * params.buffer_size
            if field in ["idx", "stamp"]:
                size = nmrk
            arr = np.zeros(size, dtype=dtype)
            setattr(cdata, field, arr.ctypes.data_as(ctypes.POINTER(ctype)))
            setattr(cdata, field + "_ref", arr)

        cdata.npoint = params.buffer_size
        cdata.interval = params.interval

        poincaredata = {
            "toroidal": params.toroidal_angles,
            "poloidal": params.poloidal_angles,
            "radial": params.radial_distances,
            }
        for poincare, data in poincaredata.items():
            if data is None:
                setattr(cdata, "n" + poincare, 0)
                continue

            setattr(cdata, "n" + poincare, len(data))
            arr = np.zeros(len(data), dtype=np.float64)
            np.copyto(arr, data)
            setattr(cdata, poincare, arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
            setattr(cdata, poincare + "_ref", arr)

        orbit = cls()
        orbit._cdata = cdata
        return orbit


    def save(self, file):
        """"""
        for field in ["r", "z", "phi", "p1", "p2", "p3", "mileage",
                      "id", "charge", "poincare", "simmode"]:
            file.write(field, getattr(self, field))
        for field in ["stamp", "idx"]:
            file.write(field, getattr(self._cdata, field + "_ref"))

        del self._cdata
        self._file = file


    def extract_fields(self, *quantities, filter=None, poincare=None):
        """Return recorded orbit data.

        This function returns the actual data that was recorded during the run.
        The data maybe stored in the struct or in the file. The data is returned
        as an array ordered by marker ID (major) and mileage (minor).

        The evaluated quantity at a given position corresponds to that mode
        which was active at the simulation:

        - GO simulations return particle phase-space.
        - GC guiding-center simulations return guiding-center phase-space.
        - In hybrid simulations, the quantity corresponds to the simulation mode
          at the moment the orbit position was recorded.

        Parameters
        ----------
        *quantities : str
            Names of the queried quantities.
        filter : array_like, optional
            Return values for specific marker IDs.
        poincare : int, optional
            Return only points that are on the given Poincaré plane.

        Returns
        -------
        **values : dict[str, array_like]
            Dictionary values of the quantities as an array ordered by marker ID
            (major) and mileage (minor).
        """
        ids = self.id
        mileage = self.mileage

        order = np.lexsort((mileage, ids))
        if filter is not None:
            filter = np.asarray(filter, dtype=ids.dtype)
            mask = np.isin(ids[order], filter)
        else:
            mask = ids[order] > 0
        if poincare is not None:
            mask &= self.poincare == poincare
        order = order[mask]

        arrays = {}
        for name in quantities:
            arrays[name] = getattr(self, name)[order]

        return arrays

    def evaluate_quantity(self, mode, state, bfield, *quantities, filter=None, poincare=None):
        """Calculate orbit quantities using the helper functions and data.

        Parameters
        ----------
        mode : str
            Simulation mode.
        state
            State data used in evaluating derived quantities.
        bfield
            Magnetic field data used in evaluating derived quantities.
        *quantities : str
            Names of the quantities.
        filter : array_like, optional
            Return values for specific marker IDs.
        poincare : int, optional
            Return only points that are on the given Poincaré plane.

        Returns
        -------
        *value : array_like
            The quantities as arrays ordered by marker ID (major) and mileage
            (minor).
        """
        # Check what fields are available based on the simulation mode
        fields = ["r", "z", "phi", "mileage", "id",]
        if mode == "particle":
            fields += ("pr", "pphi", "pz", "charge", "mass")
        elif mode == "guidingcenter":
            fields += ("ppar", "mu", "zeta", "charge", "mass")
        elif mode == "hybrid":
            pass

        # Add magnetic field to available quantities
        bfield_fundamentals = (
            "psi", "rho", "br", "bphi", "bz",
        )

        # Filter non-physical
        extract = []
        physical_quantities = []
        for qnt in quantities:
            if qnt in ["id", "poincare", "simmode", ]:
                extract.append(qnt)
            elif qnt in ["time", "connectionlength",]:
                extract.append("mileage")
                extract.append("id")
            else:
                physical_quantities.append(qnt)

        # Resolve requested quantities
        fundamentals = fields + list(bfield_fundamentals)
        needed, computethese = formulas.resolve_quantities(fundamentals, physical_quantities)

        map_quantity_to_field = {
            "particle": {
                "p1": "pr",
                "p2": "pphi",
                "p3": "pz",
            },
            "guidingcenter": {
                "p1": "ppar",
                "p2": "mu",
                "p3": "zeta",
            },
            "fieldline": {
                "p1": "pitch",
            },
        }

        evaluate_bfield = []
        for qnt in needed:
            if qnt in fields and qnt != "mass":
                if qnt in map_quantity_to_field[mode].values():
                    qnt = [q for q, v in map_quantity_to_field[mode].items() if v == qnt][0]
                extract.append(qnt)
            elif qnt in bfield_fundamentals:
                evaluate_bfield.append(qnt)

        if len(evaluate_bfield) > 0:
            for qnt in ["r", "phi", "z", "mileage", "id"]:
                if qnt not in extract:
                    extract.append(qnt)

        out = self.extract_fields(*extract, filter=filter, poincare=poincare)
        for qnt in list(out.keys()):
            if qnt in map_quantity_to_field[mode]:
                name = map_quantity_to_field[mode][qnt]
                if name in ["pr", "pphi", "pz", "ppar",]:
                    out[qnt] *= unyt.kg * unyt.m / unyt.s
                elif name in ["mu",]:
                    out[qnt] *= unyt.eV / unyt.T
                elif name in ["zeta",]:
                    out[qnt] *= unyt.rad
                out[name] = out[qnt]
                del out[qnt]

        # Set time
        if "time" in quantities or len(evaluate_bfield) > 0:
            idx = np.searchsorted(state.ids, out["id"])
            out["time"] = state.time[idx] + (mode != "fieldline") * out["mileage"] * unyt.s
        if "connectionlength" in quantities:
            out["connectionlength"] = 0.0

        if len(evaluate_bfield) > 0:
            bout = evaluate(out["r"], out["phi"], out["z"], out["time"], *evaluate_bfield, bfield=bfield)
            for qnt, val in zip(evaluate_bfield, bout):
                out[qnt] = val

        # Get mass and time from state
        if mode != "fieldline" and "mass" in needed:
            out["mass"] = state.mass[0]

        # Calculate requested quantities
        computed = formulas.compute_quantities(out, computethese)

        # Organize and return
        out = out | computed
        arr = []
        for qnt in quantities:
            if qnt == "mileage":
                if mode == "fieldline":
                    out[qnt] *= unyt.m
                else:
                    out[qnt] *= unyt.s

            arr.append(out[qnt])

        if len(arr) == 1:
            return arr[0]
        return tuple(arr)
