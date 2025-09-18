"""Defines Wall2D first wall model input class and the corresponding factory
method.
"""
import ctypes
from typing import Dict, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from ..access import _variants, InputVariant, Format, TreeCreateClassMixin
from ...libascot import LIBASCOT
from ...exceptions import AscotIOException


class Wall2D(InputVariant):
    """Simple wall model where the wall is assumed to be axisymmetric."""

    # pylint: disable=too-few-public-methods
    class Struct(ctypes.Structure):
        """Python wrapper for the struct in wall_2d.h."""
        _pack_ = 1
        _fields_ = [
            ('n', ctypes.c_int32),
            ('PADDING_0', ctypes.c_ubyte * 4),
            ('wall_r', ctypes.POINTER(ctypes.c_double)),
            ('wall_z', ctypes.POINTER(ctypes.c_double)),
            ('flag', ctypes.POINTER(ctypes.c_int32)),
            ]


    def __init__(self, qid, date, note) -> None:
        super().__init__(
            qid=qid, date=date, note=note, variant="Wall2D",
            struct=Wall2D.Struct(),
            )
        self._r: unyt.unyt_array
        self._z: unyt.unyt_array
        self._flag: unyt.unyt_array
        self._labels: Dict[str,int]

    @property
    def r(self) -> unyt.unyt_array:
        r""":math:`R` coordinates of the wall vertices."""
        if self._staged:
            n = self._from_struct_("n", shape=())
            return self._from_struct_("wall_r", shape=(n,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("r")
        return self._r.copy()

    @property
    def z(self) -> unyt.unyt_array:
        r""":math:`z` coordinates of the wall vertices."""
        if self._staged:
            n = self._from_struct_("n", shape=())
            return self._from_struct_("wall_z", shape=(n,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("z")
        return self._z.copy()

    @property
    def flag(self) -> unyt.unyt_array:
        r"""Integer label for grouping wall elements together."""
        if self._staged:
            n = self._from_struct_("n", shape=())
            return self._from_struct_("flag", shape=(n,), units="")
        if self._format == Format.HDF5:
            return self._read_hdf5("flag")
        return self._flag.copy()

    @property
    def labels(self) -> Dict[str,int]:
        r"""Integer label for grouping wall elements together."""
        if self._staged:
            # This information is not currently stored/used within libascot.so
            return self._labels
        if self._format == Format.HDF5:
            return {key.decode("utf-8"): int(value) for key, value in
                    zip(self._read_hdf5("labelkeys"),
                        self._read_hdf5("labelvalues"))
                    }
        return self._labels.copy()

    def _export_hdf5(self):
        """Export data to HDF5 file."""
        if self._format == Format.HDF5:
            raise AscotIOException("Data is already stored in the file.")
        data = self.export()
        data["labelkeys"] = list(data["labels"].keys())
        data["labelvalues"] = list(data["labels"].values())
        del data["labels"]
        self._treemanager.hdf5manager.write_datasets(
            self.qid, self.variant, data,
            )
        self._format = Format.HDF5

    def export(self):
        data = {
            "r": self._r,
            "z": self._z,
            "flag": self._flag,
            "labels": self._labels,
        }
        return data

    def stage(self):
        init = LIBASCOT.wall_2d_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_int32),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                self.r.size,
                self.r.v,
                self.z.v,
                self.flag,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._r
                del self._z
                del self._flag
            self._staged = True

    def unstage(self):
        free = LIBASCOT.wall_2d_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._r = self.r
                self._z = self.z
                self._flag = self.flag
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateWall2DMixin(TreeCreateClassMixin):
    """Mixin class used by `Data` to create Wall2D input."""

    #pylint: disable=protected-access, too-many-arguments
    def create_wall2d(
            self,
            r: np.ndarray | None = None,
            z: np.ndarray | None = None,
            flag: Optional[np.ndarray] = None,
            labels: Optional[Dict[str,int]] = None,
            note: Optional[str] = None,
            activate: bool = False,
            dryrun: bool = False,
            store_hdf5: Optional[bool] = None,
            ) -> Wall2D:
        r"""Create closed 2D polygon that forms the (axisymmetric) first wall.

        This is a simple wall model which is easy to create when 3D wall is not
        available or if the exact number of losses or wall loads are not of
        interest. Performance-wise, there is little difference between the 2D
        and 3D wall models except that the latter may consume significantly more
        memory if the mesh is fine.

        Parameters
        ----------
        r : array_like (n,)
            :math:`R` coordinates of the wall vertices.
        z : array_like (n,)
            :math:`z` coordinates of the wall vertices.
        flag : array_like (n,), optional
            Integer label for grouping wall elements together.

            By grouping elements together (e.g. all elements that belong to
            the divertor form one group), investigating losses and wall loads
            component-wise becomes easier in post-processing.
        labels : dict[str,int], optional
            Human readable labels for the flag values.

            For example ``labels`` = ``{"divertor":1, "firstwall":2}``.
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
        inputdata : ~a5py.data.wall.Wall2D
            Freshly minted input data object.

        Notes
        -----
        During the simulation when using this input, wall impacts are checked
        with the point-in-polygon method known as *modified axis crossing
        method* [1]_. In a nutshell, this method works by moving the origo to
        the marker's :math:`(R,z)`-position. Next, the number of wall segments
        crossing the positive :math:`x`-axis (i.e. :math:`R`-axis on our case)
        are calculated. If the result is odd, the point is inside the polygon.

        This method is fast but assumes that the wall is closed. Once wall
        impact is detected, a more slower algorithm is run to determine which
        wall segment was crossed.

        [1] D.G. Alciatore, R. Miranda. A Winding Number and Point-in-Polygon
            Algorithm. Technical report, Colorado State University, 1995.
            https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf
        """
        # This function can't deal with dictionaries, so use None to omit it.
        parameters = _variants.parse_parameters(
            r, z, flag, None,
        )
        parameters["labels"] = labels

        default_r = np.array([0.0, 1.0, 1.0, 0.0])
        default_z = np.array([0.0, 0.0, 1.0, 1.0])
        n = default_r.size if parameters["r"] is None else parameters["r"].size
        _variants.validate_required_parameters(
            parameters,
            names=["r", "z"],
            units=["m", "m"],
            shape=[(n,), (n,)],
            dtype="f8",
            default=[default_r, default_z],
        )
        _variants.validate_optional_parameters(
            parameters, ["flag"], [""], [(n,)], "i4",
            [np.zeros(parameters["r"].shape, dtype="i4")],
        )
        if parameters["labels"] is None:
            parameters["labels"] = {"wall": 0}
        for label, id in parameters["labels"].items():
            if not isinstance(label, str):
                raise ValueError("Labels must be strings.")
            if not isinstance(id, int):
                raise ValueError("Flags in `labels` must be integers.")
        for id in np.unique(parameters["flag"]):
            if id not in parameters["labels"].values():
                raise ValueError(f"Flag {id} does not have a label.")
        meta = _variants.new_metadata("Wall2D", note=note)
        obj = self._treemanager.enter_input(
            meta, activate=activate, dryrun=dryrun, store_hdf5=store_hdf5,
            )
        for parameter, value in parameters.items():
            setattr(obj, f"_{parameter}", value)
            # No need to make dictionary immutable as copy is always returned
            if parameter != "labels":
                getattr(obj, f"_{parameter}").flags.writeable = False

        if store_hdf5:
            obj._export_hdf5()
        return obj
