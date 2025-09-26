"""Defines Wall3D triangular wall mesh input class and the corresponding factory
method.
"""
import ctypes
from typing import Dict, Optional

import unyt
import numpy as np

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in wall_3d.h."""

    _fields_ = [
        ('n', ctypes.c_int32),
        ('PADDING_0', ctypes.c_ubyte * 4),
        ('xmin', ctypes.c_double),
        ('xmax', ctypes.c_double),
        ('xgrid', ctypes.c_double),
        ('ymin', ctypes.c_double),
        ('ymax', ctypes.c_double),
        ('ygrid', ctypes.c_double),
        ('zmin', ctypes.c_double),
        ('zmax', ctypes.c_double),
        ('zgrid', ctypes.c_double),
        ('depth', ctypes.c_int32),
        ('ngrid', ctypes.c_int32),
        ('wall_tris', ctypes.POINTER(ctypes.c_double)),
        ('flag', ctypes.POINTER(ctypes.c_int32)),
        ('tree_array', ctypes.POINTER(ctypes.c_int32)),
        ('tree_array_size', ctypes.c_int32),
        ('PADDING_1', ctypes.c_ubyte * 4),
        ]


@Leaf.register
class Wall3D(InputVariant):
    """Wall mesh made of triangular elements."""

    @property
    def vertices(self) -> unyt.unyt_array:
        r""":math:`z` coordinates of the wall vertices."""
        if self._staged:
            n = self._from_struct_("n", shape=())
            return self._from_struct_("wall_tris", shape=(9*n,), units="m")
        if self._format == Format.HDF5:
            return self._read_hdf5("vertices")

    @property
    def flag(self) -> unyt.unyt_array:
        r"""Integer label for grouping wall elements together."""
        if self._staged:
            n = self._from_struct_("n", shape=())
            return self._from_struct_("flag", shape=(n,), units="")
        if self._format == Format.HDF5:
            return self._read_hdf5("flag")

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
            "vertices": self._vertices,
            "flag": self._flag,
            "labels": self._labels,
        }
        return data

    def stage(self):
        init = LIBASCOT.wall_3d_init
        init.restype = ctypes.c_int32
        init.argtypes = [
            ctypes.POINTER(__class__.Struct),
            ctypes.c_int32,
            ndpointer(ctypes.c_double),
            ndpointer(ctypes.c_int32),
            ]
        if not self._staged:
            if init(
                ctypes.byref(self._struct_),
                int(self.vertices.size / 9),
                self.vertices.v,
                self.flag,
            ):
                raise AscotIOException("Failed to stage data.")
            if self._format is Format.MEMORY:
                del self._vertices
                del self._flag
            self._staged = True

    def unstage(self):
        free = LIBASCOT.wall_3d_free
        free.restype = None
        free.argtypes = [ctypes.POINTER(__class__.Struct)]

        if self._staged:
            if self._format is Format.MEMORY:
                self._vertices = self.vertices
                self._flag = self.flag
            free(ctypes.byref(self._struct_))
            self._staged = False


# pylint: disable=too-few-public-methods
class CreateWall3DMixin(TreeMixin):
    """Mixin class used by `Data` to create Wall3D input."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_wall3d(
            self,
            vertices: np.ndarray=None,
            flag: Optional[np.ndarray]=None,
            labels: Optional[dict[str,int]]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> Wall3D:
        r"""Create triangular wall mesh represting the 3D wall.

        This input is a 3D wall model where the wall consists of small triangles
        that form a surface mesh. A wall collision happens when a marker
        intersects a wall triangle's plane.

        Parameters
        ----------
        vertices : array_like (9*n,)
            Each triangle's vertices' :math:`(x,y,z)` coordinates.

            The format of the array is :math:`[i*9 + j*3 + k]`, where :math:`i`
            is the index of the triangle, :math:`j` index of the vertex, and
            :math:`k \in [1,2,3]` marks the coordinate :math:`(x,y,z)`.
        flag : array_like (n,1), optional
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
        preview : bool, *optional*
            If True, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : ~a5py.data.wall.Wall3D
            Freshly minted input data object.

        Notes
        -----
        During initialization, the volume that the wall occupies is divided into
        smaller volumes using *octree*. Collision checks are only made with
        respect to triangles that are in the cells where the marker has been
        present in the ongoing time-step.
        """
        # This function can't deal with dictionaries, so use None to omit it.
        parameters = _variants.parse_parameters(
            vertices, flag, None,
        )
        parameters["labels"] = labels

        default_xyz = np.array([1., 0., 0., 2., 0., -1., 2., 0., 1.])
        if parameters["vertices"] is not None:
            n = int(parameters["vertices"].size / 9)
            if n != parameters["vertices"].size / 9:
                raise ValueError("´vertices´ have an invalid format.")
        else:
            n = int(default_xyz.size / 9)
        _variants.validate_required_parameters(
            parameters, ["vertices"], ["m"], [(9*n,)], "f8", [default_xyz],
        )
        _variants.validate_optional_parameters(
            parameters, ["flag"], [""], [(n,)], "i4",
            [np.zeros(parameters["vertices"].shape, dtype="i4")],
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
        meta = _variants.new_metadata("Wall3D", note=note)
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
