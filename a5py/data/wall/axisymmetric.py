"""Defines the 2D wall input and the corresponding factory method.
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
from a5py.data.access import InputVariant, Leaf, TreeMixin


# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    """Python wrapper for the struct in wall_2d.h."""

    _fields_ = [
        ("n", ctypes.c_int32),
        ("wall_r", ctypes.POINTER(ctypes.c_double)),
        ("wall_z", ctypes.POINTER(ctypes.c_double)),
        ("flag", ctypes.POINTER(ctypes.c_int32)),
        ]


init_fun(
    "wall_2d_init",
    ctypes.POINTER(Struct),
    ctypes.c_int32,
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_double),
    ndpointer(ctypes.c_int32),
    restype=ctypes.c_int32,
    )

init_fun("wall_2d_free")

@Leaf.register
class Wall2D(InputVariant):
    """Simple wall model where the wall is assumed to be axisymmetric."""

    _labels: dict[str, int] | None
    """Human readable labels for the flag values.

    These are not needed and, hence, not stored in libascot.so. Therefore we
    store them here when this instance is staged.
    """

    @property
    def n(self) -> int:
        """Number of wall vertices."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("n", ())
        assert self._file is not None
        return self._file.read("r").size

    @property
    def r(self) -> unyt.unyt_array:
        r""":math:`R` coordinates of the wall vertices."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("wall_r", (self.n,), "m")
        assert self._file is not None
        return self._file.read("r")

    @property
    def z(self) -> unyt.unyt_array:
        r""":math:`z` coordinates of the wall vertices."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("wall_z", (self.n,), "m")
        assert self._file is not None
        return self._file.read("z")

    @property
    def flag(self) -> np.ndarray:
        r"""Integer label for grouping wall elements."""
        if self._cdata is not None:
            return self._cdata.readonly_carray("flag", shape=(self.n,))
        assert self._file is not None
        return self._file.read("flag")

    @property
    def labels(self) -> dict[str,int]:
        r"""Human readable labels for the flag values."""
        if self._cdata is not None:
            assert self._labels is not None
            return copy.deepcopy(self._labels)
        assert self._file is not None
        return {label.decode("utf-8"): int(flag) for label, flag in
                zip( self._file.read("labelkeys"), self._file.read("labelvalues") )
                }

    # pylint: disable=too-many-arguments
    def _stage(
            self, r: unyt.unyt_array,
            z: unyt.unyt_array,
            flag: np.ndarray,
            labels: dict[str, int],
            ) -> None:
        self._cdata = Struct()
        self._labels = labels
        if LIBASCOT.wall_2d_init(
            ctypes.byref(self._cdata), r.size, r.v, z.v, flag,
            ):
            self._cdata = None
            self._labels = None
            raise AscotMeltdownError("Could not initialize struct.")

    def _save_data(self) -> None:
        assert self._file is not None
        self._file.write("r", self.r)
        self._file.write("z", self.z)
        self._file.write("flag", self.flag)

        labels = np.char.encode(list( self.labels.keys() ), "utf-8")
        flags = np.array( list(self.labels.values()) )
        self._file.write("labelkeys", labels)
        self._file.write("labelvalues", flags)

    def stage(self) -> None:
        super().stage()
        self._stage(**self.export())

    def unstage(self) -> None:
        super().unstage()
        assert self._cdata is not None
        LIBASCOT.wall_2d_free(ctypes.byref(self._cdata))
        self._cdata = None
        self._labels = None


# pylint: disable=too-few-public-methods
class CreateWall2DMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_wall2d(
            self,
            r: utils.ArrayLike,
            z: utils.ArrayLike,
            flag: Optional[utils.ArrayLike]=None,
            labels: Optional[dict[str,int]]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> Wall2D:
        r"""Create closed 2D polygon that forms the (axisymmetric) first wall.

        This is a simple wall model which is easy to create when 3D wall is not
        available or if the exact number of losses or wall loads are not of
        interest. Performance-wise, there is little difference between the 2D
        and 3D wall models except that the latter may consume significantly more
        memory if the mesh has fine resolution.

        Parameters
        ----------
        r : array_like (n,)
            :math:`R` coordinates of the wall vertices.
        z : array_like (n,)
            :math:`z` coordinates of the wall vertices.
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
        obj : :class:`.Wall2D`
            Input variant created from the given parameters.

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

        .. [1] D.G. Alciatore, R. Miranda. A Winding Number and Point-in-Polygon
           Algorithm. Technical report, Colorado State University, 1995.
           https://www.engr.colostate.edu/~dga/documents/papers/point_in_polygon.pdf
        """
        if not isinstance(r, Sized) or not isinstance(z, Sized):
            raise ValueError(
                "'r' and 'z' must be arrays with at least three elements."
                )
        nr, nz = utils.size(r), utils.size(z)
        if nr != nz:
            raise ValueError("'r' and 'z' must have same size.")
        n = nr
        if n < 3:
            raise ValueError(
                "'r' and 'z' must be arrays with at least three elements."
                )

        with utils.validate_variables() as v:
            r = v.validate("r", r, (n,), "m")
            z = v.validate("z", z, (n,), "m")
            flag = v.validate(
                "flag", flag, (n,), dtype="i4", default=np.full(n, 0),
                )

        if labels is None:
            labels = {"wall": 0}

        for label, label_flag in labels.items():
            if not isinstance(label, str):
                raise ValueError("Labels must be strings.")
            if not isinstance(label_flag, int):
                raise ValueError("Flags in 'labels' must be integers.")

        assert isinstance(flag, np.ndarray)
        unlabeled = [f for f in flag if f not in labels.values()]
        if unlabeled:
            raise ValueError(
                f"Flags `{unlabeled}` do not have a corresponding label."
                )

        leaf = Wall2D(note=note)
        leaf._stage(r=r, z=z, flag=flag, labels=copy.deepcopy(labels))
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="wall",
            )
        return leaf
