"""Defines magnetic field line marker input class and the corresponding factory
method.
"""
from typing import Optional, Any, Sized, cast

import unyt
import numpy as np

from a5py import utils
from a5py.data.access import InputVariant, Leaf, TreeMixin

from .state import Structure


@Leaf.register
class FieldlineMarker(InputVariant):
    """Marker input in magnetic field line (3D) phase-space."""

    @property
    def ids(self) -> np.ndarray:
        """Unique identifier for each marker."""
        if self._cdata is not None:
            out = np.zeros((self.n,), dtype="i8")
            for i in range(out.size):
                out[i] = self._cdata[i].id
            return out
        assert self._file is not None
        return self._file.read("ids")

    @property
    def r(self) -> unyt.unyt_array:
        r"""Field-line :math:`R` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].r
            return out
        assert self._file is not None
        return self._file.read("r")

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Field-line :math:`\phi` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "rad")
            for i in range(out.size):
                out[i] = self._cdata[i].phi

            out.convert_to_units("deg")
            return out
        assert self._file is not None
        return self._file.read("phi")

    @property
    def z(self) -> unyt.unyt_array:
        r"""Field-line :math:`z` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].z
            return out
        assert self._file is not None
        return self._file.read("z")

    @property
    def direction(self) -> unyt.unyt_array:
        r"""The direction this marker travels along the field line.

        A negative value means opposite to the magnetic field vector.
        """
        if self._cdata is not None:
            out = np.zeros((self.n,), dtype="f8")
            for i in range(out.size):
                out[i] = self._cdata[i].pitch
            return out
        assert self._file is not None
        return self._file.read("direction")

    @property
    def time(self) -> unyt.unyt_array:
        """Time instant (e.g. with respect to the plasma pulse) when the marker
        simulation is launched.
        """
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "s")
            for i in range(out.size):
                out[i] = self._cdata[i].time
            return out
        assert self._file is not None
        return self._file.read("time")

    @property
    def n(self) -> int:
        """Number of markers."""
        if self._cdata is not None:
            return len(self._cdata)
        assert self._file is not None
        return self._file.read("ids").size

    def _fill_state(self, bfield):
        pass

    #pylint: disable=too-many-arguments
    def _stage(
            self, r: unyt.unyt_array,
            z: unyt.unyt_array,
            phi: unyt.unyt_array,
            direction: unyt.unyt_array,
            time: unyt.unyt_array,
            ids: np.ndarray,
            ) -> None:
        self._cdata = (Structure * ids.size)()
        for i in range(ids.size):
            setattr(self._cdata[i], "r", r[i])
            setattr(self._cdata[i], "z", z[i])
            setattr(self._cdata[i], "id", ids[i])
            setattr(self._cdata[i], "phi", phi[i].to("rad"))
            setattr(self._cdata[i], "time", time[i])
            setattr(self._cdata[i], "pitch", direction[i])

    def _save_data(self) -> None:
        assert self._file is not None
        for field in ["r", "z", "phi", "direction", "time", "ids"]:
            self._file.write(field, getattr(self, field))

    def export(self) -> dict[str, np.ndarray | unyt.unyt_array]:
        fields = ["r", "z", "phi", "direction", "time", "ids"]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            r=self.r, z=self.z, phi=self.phi, direction=self.direction,
            time=self.time, ids=self.ids,
        )

    def unstage(self) -> None:
        super().unstage()
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_fieldlinemarker(
            self,
            r: utils.ArrayLike,
            z: utils.ArrayLike,
            phi: Optional[utils.ArrayLike]=None,
            direction: Optional[utils.ArrayLike]=None,
            time: Optional[utils.ArrayLike]=None,
            ids: Optional[utils.ArrayLike]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> FieldlineMarker:
        r"""Create markers representing magnetic field lines.

        Except for ``ids``, all input arrays can also be scalars in which case
        the same value is used for all markers. The number of markers is deduced
        from the length of the input arrays, which must be consistent: if one
        has length `n`, the others must be scalars or also have length
        `n`.

        Parameters
        ----------
        r : array_like (n,) *or* scalar
            Field line :math:`R` coordinate.
        z : array_like (n,) *or* scalar
            Field line :math:`z` coordinate.
        phi : array_like (n,) *or* scalar, *optional*
            Field line :math:`\phi` coordinate.

            If not given, the values are randomly sampled from a uniform
            distribution between 0 and 360 degrees.
        direction : array_like (n,) *or* scalar, *optional*
            The direction this marker travels along the field line.

            The default value is 1. A negative value means opposite to the
            magnetic field vector.
        time : array_like (n,) *or* scalar, *optional*
            Time instant (e.g. with respect to the plasma pulse) when the marker
            simulation is launched.

            Relevant mainly for time-dependent simulations. Default value is
            zero.
        ids : array_like (n,) *or* scalar, *optional*
            Unique identifier for each marker (must be a positive integer).

            If not given, the IDs are 1, 2, 3, ..., `n`.
        note : *optional*
                A short note to document this data.

                The first word of the note is converted to a tag which you can use
                to reference the data.
        activate : *optional*
            Set this input as active on creation.
        preview : *optional*
            If ``True``, the input is created but it is not included in the data
            structure nor saved to disk.

            The input cannot be used in a simulation but it can be previewed.
        save : *optional*
            Store this input to disk.

        Returns
        -------
        inputdata : FieldlineMarker
            Input variant created from the given parameters.

        Notes
        -----
        This input can only be used in magnetic field-line tracing since the
        markers don't have physical quantities such as mass or charge.
        """
        n = 1
        for param in [r, z, phi, direction, time, ids]:
            if isinstance(param, Sized):
                if n != 1 and n != utils.size(param):
                    raise ValueError(
                        "Input arrays have inconsistent sizes."
                        )
                n = utils.size(param)

        def to_array(x: Any) -> unyt.unyt_array | np.ndarray:
            """Convert x to an array of length n."""
            if x is None:
                return None
            if not isinstance(x, Sized):
                return np.full(n, x)

            if utils.size(x) == 1:
                if hasattr(x, "units"):
                    return np.full(n, x) * x.units
                return np.full(n, x)
            return x

        r, z, phi, direction, time = (
            to_array(r), to_array(z), to_array(phi), to_array(direction),
            to_array(time),
            )

        with utils.validate_variables() as v:
            r = v.validate("r", r, (n,), "m")
            z = v.validate("z", z, (n,), "m")
        with utils.validate_variables() as v:
            ids = v.validate(
                "ids", ids, (n,), dtype="i8", default=np.arange(1, n+1),
                )
            phi = v.validate(
                "phi", phi, (n,), "deg", default=np.random.rand(n) * 360
                )
            direction = v.validate(
                "direction", direction, (n,), default=np.ones((n,))
                )
            time = v.validate("time", time, (n,), "s", default=np.zeros((n,)))

        leaf = FieldlineMarker(note=note)
        leaf._stage(
            r=r, z=z, phi=phi, direction=direction, time=time,
            ids=cast(np.ndarray, ids),
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="marker",
            )
        return leaf
