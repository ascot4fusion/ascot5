"""Python wrapper for libascot.so and its helper functions and classes.

Importing this module will load libascot.so.
"""

from __future__ import annotations

import ctypes
from typing import Optional, overload
from pathlib import Path

import unyt
import numpy as np

LIBASCOT: ctypes.CDLL
"""The ctypes.CDLL object for libascot.so."""


def init_fun(name, *argtypes, restype=None):
    """Initialize a function in libascot.so."""
    fun = getattr(LIBASCOT, name)
    fun.argtypes = argtypes
    fun.restype = restype


def input_category(cls):

    def use(self, variant):
        for enumval, field in enumerate(cls._fields_):
            if field[0] in variant.variant.lower():
                self.type = ctypes.c_int32(1 + enumval) # Type enum starts at 1
                setattr(self, field[0], ctypes.pointer(variant._cdata))
                break

    cls.use = use
    return cls


def _find_libascot():
    """Find and load libascot.so.

    If the code was installed with a package manager, libascot.so should be
    found in ``LD_LIBRARY_PATH``. If this is a local install, then libascot.so
    can be found in the relative path "../build/libascot.so" assuming that
    ``a5py`` was installed with ``pip install -e .`` option.

    This function is called when this module is initialized.
    """
    global LIBASCOT
    err = 0
    libpath = str(Path(__file__).absolute().parent.parent) + "/build/libascot.so"
    try:
        LIBASCOT = ctypes.CDLL(libpath)
    except OSError as error:
        # Either libascot.so was not in "../build/" or it could not be loaded
        # due to missing dependencies.
        err = error
        # Missing dependencies error would contain the name of the missing
        # library, not libascot.so.
        unresolved_dependencies = not "libascot.so" in str(err)
        # Undefined symbol error should only concern developers, but it can be
        # a major headache.
        undefined_symbol = "undefined symbol" in str(err)

    if err:
        if unresolved_dependencies:
            raise ImportError(str(err))
        elif undefined_symbol:
            raise ImportError("Error in the source: " + str(err))
        else:
            try:
                # This looks for libascot.so in LD_LIBRARY_PATH.
                LIBASCOT = ctypes.CDLL("libascot.so")
            except OSError as error:
                final_err = error
            if final_err:
                raise ImportError(
                    "Failed to load libascot.so: make sure it is compiled "
                    "(if installed from source) or that it is included in "
                    "LD_LIBRARY_PATH."
                )


_find_libascot()


# pylint: disable=too-few-public-methods
class DataStruct(ctypes.Structure):
    """Expands the ctypes.Structure class to better support data storage."""

    @overload
    def readonly_carray(self, name: str, shape: tuple[int, ...]) -> np.ndarray: ...

    @overload
    def readonly_carray(
        self, name: str, shape: tuple[int, ...], units: str
    ) -> unyt.unyt_array: ...

    def readonly_carray(
        self,
        name: str,
        shape: tuple[int, ...],
        units: Optional[str] = None,
    ) -> np.ndarray | unyt.unyt_array:
        """Return a read-only view on the array within the struct.

        Parameters
        ----------
        name : str
            Name of the array within the struct.
        shape : tuple
            Shape of the array.
        units : str, *optional*
            Units of the array.

        Returns
        -------
        arr : np.ndarray or unyt.array
            Read-only view on the array.
        """
        # as_array only reads data based on shape but doesn't actually reshape
        # the data, which is why we need to reshape explicitly
        arr = np.ctypeslib.as_array(getattr(self, name), shape=shape).reshape(shape)
        arr.flags.writeable = False
        if units is not None:
            return unyt.unyt_array(arr, units)
        return arr

    def readonly_grid(
        self,
        coordinate: str,
        units: str,
        member: Optional[str] = None,
    ) -> unyt.unyt_array:
        """Return a read-only view on the grid within the struct.

        The uniform grids in libascot.so are stored as ``(xlim, nx)``.
        This function reads the values and converts them to (read-only) grid
        array.

        Parameters
        ----------
        coordinate : str
            Name of the grid coordinate.
        units : str
            Units of the grid coordinate.
        member : str, *optional*
            A member struct (e.g. spline struct) that contains the grid if it is
            not in the main struct.

        Returns
        -------
        grid : unyt.array
            Read-only view on the grid.
        """
        struct = self if member is None else getattr(self, member)
        grid = np.linspace(
            getattr(struct, f"{coordinate}lim")[0],
            getattr(struct, f"{coordinate}lim")[1],
            getattr(struct, f"n{coordinate}"),
        )
        grid.flags.writeable = False
        return unyt.unyt_array(grid, units)

    def readonly_interp(self, name: str, units: str, idx: int = None):
        """Read the tabulated values from interpolation struct.

        This does not read all values but only those that correspond to the
        values of the function being interpolated (i.e. other coefficients are
        ignored). The returned array is read-only.

        Parameters
        ----------
        name : str
            Name of the interpolation struct attribute.
        units : str
            Units of the array.
        idx : int, optional
            In case the attribute is an array of structs, this is the index on
            the array which is read.

        Returns
        -------
        arr : np.ndarray
            Read-only view on the array.
        """
        if idx is None:
            attribute = getattr(self, name)
        else:
            attribute = getattr(self, name)[idx]
        match attribute:
            case Linear1D():
                shape, skip = (attribute.nx,), 1
            case Linear3D():
                shape, skip = (attribute.nx, attribute.ny, attribute.nz), 1
            case Spline1D():
                shape, skip = (attribute.nx,), 2
            case Spline2D():
                shape, skip = (attribute.nx, attribute.ny), 4
            case Spline3D():
                shape, skip = (attribute.nx, attribute.ny, attribute.nz), 8
            case _:
                raise ValueError(
                    f"Unknown interpolation struct type: {type(attribute)} "
                    "Please submit an issue to the ASCOT5 GitHub repository."
                )
        arr = np.ctypeslib.as_array(
            attribute.c[: np.multiply.reduce(shape) * skip : skip],
        )
        arr.flags.writeable = False
        return unyt.unyt_array(arr.reshape(shape), units)


class Linear1D(ctypes.Structure):
    """Python wrapper for the 1D linear-interpolator defined in interp.h."""

    _fields_ = [
        ("nx", ctypes.c_size_t),
        ("xbc", ctypes.c_int32),
        ("dx", ctypes.c_double),
        ("xlim", ctypes.c_double * 2),
        ("c", ctypes.POINTER(ctypes.c_double)),
    ]


class Linear2D(ctypes.Structure):
    """Python wrapper for the 2D linear-interpolator defined in interp.h."""

    _fields_ = [
        ("nx", ctypes.c_size_t),
        ("ny", ctypes.c_size_t),
        ("xbc", ctypes.c_int32),
        ("ybc", ctypes.c_int32),
        ("dx", ctypes.c_double),
        ("dy", ctypes.c_double),
        ("xlim", ctypes.c_double * 2),
        ("ylim", ctypes.c_double * 2),
        ("c", ctypes.POINTER(ctypes.c_double)),
    ]


class Linear3D(ctypes.Structure):
    """Python wrapper for the 3D linear-interpolator defined in interp.h."""

    _fields_ = [
        ("nx", ctypes.c_size_t),
        ("ny", ctypes.c_size_t),
        ("nz", ctypes.c_size_t),
        ("xbc", ctypes.c_int32),
        ("ybc", ctypes.c_int32),
        ("zbc", ctypes.c_int32),
        ("dx", ctypes.c_double),
        ("dy", ctypes.c_double),
        ("dz", ctypes.c_double),
        ("xlim", ctypes.c_double * 2),
        ("ylim", ctypes.c_double * 2),
        ("zlim", ctypes.c_double * 2),
        ("c", ctypes.POINTER(ctypes.c_double)),
    ]


class Spline1D(ctypes.Structure):
    """Python wrapper for the 1D spline-interpolator defined in interp.h."""

    _fields_ = [
        ("nx", ctypes.c_size_t),
        ("xbc", ctypes.c_int32),
        ("dx", ctypes.c_double),
        ("xlim", ctypes.c_double * 2),
        ("c", ctypes.POINTER(ctypes.c_double)),
    ]


class Spline2D(ctypes.Structure):
    """Python wrapper for the 2D spline-interpolator defined in interp.h."""

    _fields_ = [
        ("nx", ctypes.c_size_t),
        ("ny", ctypes.c_size_t),
        ("xbc", ctypes.c_int32),
        ("ybc", ctypes.c_int32),
        ("dx", ctypes.c_double),
        ("dy", ctypes.c_double),
        ("xlim", ctypes.c_double * 2),
        ("ylim", ctypes.c_double * 2),
        ("c", ctypes.POINTER(ctypes.c_double)),
    ]


class Spline3D(ctypes.Structure):
    """Python wrapper for the 3D spline-interpolator defined in interp.h."""

    _fields_ = [
        ("nx", ctypes.c_size_t),
        ("ny", ctypes.c_size_t),
        ("nz", ctypes.c_size_t),
        ("xbc", ctypes.c_int32),
        ("ybc", ctypes.c_int32),
        ("zbc", ctypes.c_int32),
        ("dx", ctypes.c_double),
        ("dy", ctypes.c_double),
        ("dz", ctypes.c_double),
        ("xlim", ctypes.c_double * 2),
        ("ylim", ctypes.c_double * 2),
        ("zlim", ctypes.c_double * 2),
        ("c", ctypes.POINTER(ctypes.c_double)),
    ]
