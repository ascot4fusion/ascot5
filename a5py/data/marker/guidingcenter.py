"""Defines guiding center marker input class and the corresponding factory
method.
"""
from typing import Optional, Any, Sized, cast

import unyt
import numpy as np
from numpy.typing import ArrayLike

from a5py import utils
from a5py.physlib import Species
from a5py.data.access import InputVariant, Leaf, TreeMixin

from .state import Structure


@Leaf.register
class GuidingcenterMarker(InputVariant):
    """Marker input in guiding-center (5D) phase-space."""

    @property
    def species(self) -> Species:
        """Marker species."""
        if self._cdata is not None:
            return Species.from_znumanum(
                self._cdata[0].znum, self._cdata[0].anum
                )
        assert self._file is not None
        return Species.from_znumanum(
            cast(int, self._file.read("znum")),
            cast(int, self._file.read("anum")),
            )

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
    def charge(self) -> unyt.unyt_array:
        r"""Marker charge."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "e")
            for i in range(out.size):
                out[i] = self._cdata[i].charge
            return out
        assert self._file is not None
        return self._file.read("charge")

    @property
    def r(self) -> unyt.unyt_array:
        r"""Guiding center :math:`R` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].r
            return out
        assert self._file is not None
        return self._file.read("r")

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Guiding center :math:`\phi` coordinate."""
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
        r"""Guiding center :math:`z` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].z
            return out
        assert self._file is not None
        return self._file.read("z")

    @property
    def ekin(self) -> unyt.unyt_array:
        r"""Guiding center kinetic energy."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "J")
            for i in range(out.size):
                out[i] = self._cdata[i].ekin
            out.convert_to_units("eV")
            return out
        assert self._file is not None
        return self._file.read("ekin")

    @property
    def pitch(self) -> unyt.unyt_array:
        r"""Guiding center pitch."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "1")
            for i in range(out.size):
                out[i] = self._cdata[i].pitch
            return out
        assert self._file is not None
        return self._file.read("pitch")

    @property
    def gyroangle(self) -> unyt.unyt_array:
        r"""Guiding center gyro-angle."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "rad")
            for i in range(out.size):
                out[i] = self._cdata[i].zeta
            return out
        assert self._file is not None
        return self._file.read("gyroangle")

    @property
    def time(self) -> unyt.unyt_array:
        """The time-instant when the marker simulation starts."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "s")
            for i in range(out.size):
                out[i] = self._cdata[i].time
            return out
        assert self._file is not None
        return self._file.read("time")

    @property
    def weight(self) -> unyt.unyt_array:
        """How many physical particles each marker represents."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "particles/s")
            for i in range(out.size):
                out[i] = self._cdata[i].weight
            return out
        assert self._file is not None
        return self._file.read("weight")

    @property
    def n(self) -> int:
        """Number of markers."""
        if self._cdata is not None:
            return len(self._cdata)
        assert self._file is not None
        return self._file.read("ids").size

    #pylint: disable=too-many-arguments
    def _stage(
            self, r: unyt.unyt_array,
            z: unyt.unyt_array,
            phi: unyt.unyt_array,
            ekin: unyt.unyt_array,
            pitch: unyt.unyt_array,
            charge: unyt.unyt_array,
            time: unyt.unyt_array,
            ids: np.ndarray,
            gyroangle: unyt.unyt_array,
            weight: unyt.unyt_array,
            species: Species,
            ) -> None:
        self._cdata = (Structure * ids.size)()
        for i in range(ids.size):
            setattr(self._cdata[i], "r", r[i])
            setattr(self._cdata[i], "z", z[i])
            setattr(self._cdata[i], "id", ids[i])
            setattr(self._cdata[i], "phi", phi[i])
            setattr(self._cdata[i], "time", time[i])
            setattr(self._cdata[i], "ekin", ekin[i].to("J"))
            setattr(self._cdata[i], "pitch", pitch[i])
            setattr(self._cdata[i], "zeta", gyroangle[i])
            setattr(self._cdata[i], "charge", charge[i])
            setattr(self._cdata[i], "weight", weight[i])
            setattr(self._cdata[i], "anum", species.anum)
            setattr(self._cdata[i], "znum", species.znum)
            setattr(self._cdata[i], "mass", species.mass)

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "r", "z", "phi", "charge", "time", "ids", "weight", "ekin", "pitch",
            "gyroangle",
            ]:
            self._file.write(field, getattr(self, field))
        self._file.write("anum", np.array(self.species.anum))
        self._file.write("znum", np.array(self.species.znum))

    def export(self) -> dict[str, Species | np.ndarray | unyt.unyt_array]:
        fields = [
            "r", "z", "phi", "charge", "time", "ids", "weight", "ekin", "pitch",
            "gyroangle", "species",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            species=self.species, r=self.r, z=self.z, phi=self.phi,
            ekin=self.ekin, pitch=self.pitch, charge=self.charge,
            time=self.time, ids=self.ids, gyroangle=self.gyroangle,
            weight=self.weight,
        )

    def unstage(self) -> None:
        super().unstage()
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_guidingcentermarker(
            self,
            species: str | Species,
            charge: utils.Numerical,
            r: utils.Numerical,
            z: utils.Numerical,
            ekin: utils.Numerical,
            pitch: utils.Numerical,
            phi: Optional[utils.Numerical]=None,
            gyroangle: Optional[utils.Numerical]=None,
            weight: Optional[utils.Numerical]=None,
            time: Optional[utils.Numerical]=None,
            ids: Optional[utils.Numerical]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> GuidingcenterMarker:
        r"""Create marker input in guiding-center (5D) phase-space.

        This input can be used in all simulation modes as long as the marker
        is not neutral, but for gyro-orbit simulations the
        :class:`ParticleMarker` input is recommended.

        Parameters
        ----------
        species : str *or* :class:`.Species`
            Marker species.
        charge : array_like (n,) *or* scalar
            Marker charge.
        r : array_like (n,) *or* scalar
            Guiding center :math:`R` coordinate.
        z : array_like (n,) *or* scalar
            Guiding center :math:`z` coordinate.
        ekin : array_like (n,) *or* scalar
            Guiding center kinetic energy.
        pitch : array_like (n,) *or* scalar
            Guiding center pitch.
        phi : array_like (n,) *or* scalar, *optional*
            Guiding center :math:`\phi` coordinate.

            If not given, the values are randomly sampled from a uniform
            distribution between 0 and 360 degrees.
        gyroangle : array_like (n,), *optional*
            Guiding center gyro angle.

            Usually this doesn't matter as the gyro angle is not relevant for
            guiding center dynamics. This is only used if the simulation
            includes a guiding-center to particle transformation. If not given,
            the values are randomly sampled from an uniform distribution between
            0 and 2:math:`\pi` rad.
        weight : array_like (n,), *optional*
            How many physical particles each marker represents.
        time : array_like (n,), *optional*
            Time instant (e.g. with respect to the plasma pulse) when the marker
            simulation is launched.

            Relevant mainly for time-dependent simulations. Default value is
            zero.
        ids : array_like (n,), *optional*
            Unique identifier for each marker (must be a positive integer).

            If not given, the IDs are 1, 2, 3, ..., `n`.
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
        inputdata : GuidingcenterMarker
            Input variant created from the given parameters.
        """
        if isinstance(species, str):
            species = Species.from_string(species)

        n = 1
        for param in [
            r, z, phi, ekin, pitch, gyroangle, weight, time, ids, charge,
            ]:
            if isinstance(param, Sized):
                if utils.size(param) == 1:
                    continue
                if n != 1 and n != utils.size(param):
                    raise ValueError(
                        "Input arrays have inconsistent sizes."
                        )
                n = utils.size(param)

        to_array = lambda x: utils.to_array(x, n)

        r, z, phi, ekin, pitch, charge, weight, time, gyroangle = (
            to_array(r), to_array(z), to_array(phi), to_array(ekin),
            to_array(pitch), to_array(charge), to_array(weight), to_array(time),
            to_array(gyroangle),
            )

        with utils.validate_variables() as v:
            r = v.validate("r", r, (n,), "m")
            z = v.validate("z", z, (n,), "m")
            ekin = v.validate("ekin", ekin, (n,), "eV")
            pitch = v.validate("pitch", pitch, (n,), "1")
            charge = v.validate("charge", charge, (n,), "e")
        with utils.validate_variables() as v:
            ids = v.validate(
                "ids", ids, (n,), dtype="i8", default=np.arange(1, n+1),
                )
            phi = v.validate(
                "phi", phi, (n,), "deg", default=np.random.rand(n) * 360
                )
            weight = v.validate(
                "weight", weight, (n,), "particles/s", default=np.ones((n,))
                )
            time = v.validate("time", time, (n,), "s", default=np.zeros((n,)))
            gyroangle = v.validate(
                "gyroangle", gyroangle, (n,), "rad",
                default=2*np.pi*np.random.rand(n)
                )
        leaf = GuidingcenterMarker(note=note)
        leaf._stage(
            r=r, z=z, phi=phi, ekin=ekin, pitch=pitch, charge=charge, time=time,
            ids=cast(np.ndarray, ids), gyroangle=gyroangle, weight=weight,
            species=species,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="marker",
            )
        return leaf
