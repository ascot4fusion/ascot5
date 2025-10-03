"""Defines triangular 3D wall mesh input and the corresponding factory method.
"""
import copy
import ctypes
from typing import Sized, Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct, init_fun
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import Leaf, TreeMixin
from .common import WallVariant


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in wall_3d.h."""

    _fields_ = [
        ("n", ctypes.c_int32),
        ("depth", ctypes.c_int32),
        ("ngrid", ctypes.c_int32),
        ("tree_array_size", ctypes.c_int32),
        ("flag", ctypes.POINTER(ctypes.c_int32)),
        ("xmin", ctypes.c_double),
        ("xmax", ctypes.c_double),
        ("xgrid", ctypes.c_double),
        ("ymin", ctypes.c_double),
        ("ymax", ctypes.c_double),
        ("ygrid", ctypes.c_double),
        ("zmin", ctypes.c_double),
        ("zmax", ctypes.c_double),
        ("zgrid", ctypes.c_double),
        ("vertices", ctypes.POINTER(ctypes.c_double)),
        ("tree_array", ctypes.POINTER(ctypes.c_int32)),
        ]


init_fun(
    "WallTriangular3D_init",
    ctypes.POINTER(Struct),
    ctypes.c_int32,
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_int32),
    restype=ctypes.c_int32,
    )

init_fun("WallTriangular3D_free")


@Leaf.register
class WallTriangular3D(WallVariant):
    """Wall mesh made of triangular elements."""

    @property
    def n(self) -> int:
        if self._cdata is not None:
            return self._cdata.readonly_carray("n", ()) // 9
        assert self._file is not None
        return self._file.read("vertices").size // 9

    @property
    def vertices(self) -> unyt.unyt_array:
        r""":math:`(x, y, z)` coordinates of the wall vertices."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("vertices", (9*self.n,), "m")
        assert self._file is not None
        return self._file.read("vertices")

    # pylint: disable=too-many-arguments
    def _stage(
            self, vertices: unyt.unyt_array,
            flag: np.ndarray,
            labels: dict[str, int],
            ) -> None:
        self._cdata = Struct()
        self._labels = labels
        if LIBASCOT.WallTriangular3D_init(
            ctypes.byref(self._cdata), vertices.size, vertices.v, flag,
            ):
            self._cdata = None
            self._labels = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("vertices", self.vertices)
        self._save_flags_and_labels()

    def export(self) -> dict[str, unyt.unyt_array | np.ndarray]:
        fields = ["vertices", "flag", "labels",]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.WallTriangular3D_free(ctypes.byref(self._cdata))
        self._cdata = None
        self._labels = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_walltriangular3d(
            self,
            vertices: utils.ArrayLike,
            flag: Optional[np.ndarray]=None,
            labels: Optional[dict[str,int]]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> WallTriangular3D:
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
        flag : array_like (n,), *optional*
            Integer label for grouping wall elements.

            By grouping elements together (e.g. all elements that belong to
            the divertor form one group), investigating losses and wall loads
            component-wise becomes easier in post-processing. If ``flag`` is
            given, corresponding ``labels`` must be provided as well. By default
            all elements are simply labeled as "{wall: 0}".
        labels : dict[str,int], *optional*
            Human readable labels for the flag values.

            For example ``labels = {"divertor": 1, "firstwall": 2}``.
        note : str, *optional*
            A short note to document this data.

            The first word of the note is converted to a tag which you can use
            to reference the data.
        activate : bool, *optional*
            Set this input as active on creation.
        preview : bool, *optional*
            If ``True``, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : bool, *optional*
            Store this input to disk.

        Returns
        -------
        obj : :class:`.WallTriangular3D`
            Input variant created from the given parameters.

        Notes
        -----
        During initialization, the volume that the wall occupies is divided into
        smaller volumes using *octree*. Collision checks are only made with
        respect to triangles that are in the cells where the marker has been
        present in the ongoing time-step.
        """
        if not isinstance(vertices, Sized):
            raise ValueError(
                "'vertices' must be an array with at least nine elements."
                )
        nvertices, ntriangles = utils.size(vertices), vertices.size // 9
        if nvertices % 9 != 0:
            raise ValueError(
                "Number of elements in 'vertices' must be a multiple of 9.")

        with utils.validate_variables() as v:
            vertices = v.validate("vertices", vertices, (nvertices,), "m")
            flag = v.validate(
                "flag", flag, (ntriangles,), dtype="i4", default=np.full(ntriangles, 0),
                )

        labels = WallVariant.validate_labels(labels, flag)
        leaf = WallTriangular3D(note=note)
        leaf._stage(vertices=vertices, flag=flag, labels=copy.deepcopy(labels))
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="wall",
            )
        return leaf
