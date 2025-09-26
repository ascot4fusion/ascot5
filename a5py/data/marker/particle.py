"""Defines particle (gyro-orbit) input class and the corresponding factory
method.
"""
from typing import Optional, Any, Sized, cast

import unyt
import numpy as np

from a5py import utils
from a5py.physlib import Species
from a5py.data.access import InputVariant, Leaf, TreeMixin

from .state import Structure

@Leaf.register
class ParticleMarker(InputVariant):
    """Marker input in particle (6D) phase-space."""

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
        r"""Particle :math:`R` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].r
            return out
        assert self._file is not None
        return self._file.read("r")

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Particle :math:`\phi` coordinate."""
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
        r"""Particle :math:`z` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m")
            for i in range(out.size):
                out[i] = self._cdata[i].z
            return out
        assert self._file is not None
        return self._file.read("z")

    @property
    def vr(self) -> unyt.unyt_array:
        r"""Particle velocity :math:`R` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pr
            return out
        assert self._file is not None
        return self._file.read("vr")

    @property
    def vphi(self) -> unyt.unyt_array:
        r"""Particle velocity :math:`\phi` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pphi
            return out
        assert self._file is not None
        return self._file.read("vphi")

    @property
    def vz(self) -> unyt.unyt_array:
        r"""Particle velocity :math:`z` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pz
            return out
        assert self._file is not None
        return self._file.read("vz")

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
            charge: unyt.unyt_array,
            time: unyt.unyt_array,
            ids: np.ndarray,
            vr: unyt.unyt_array,
            vphi: unyt.unyt_array,
            vz: unyt.unyt_array,
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
            setattr(self._cdata[i], "pr", vr[i])
            setattr(self._cdata[i], "pphi", vphi[i])
            setattr(self._cdata[i], "pz", vz[i])
            setattr(self._cdata[i], "charge", charge[i])
            setattr(self._cdata[i], "weight", weight[i])
            setattr(self._cdata[i], "anum", species.anum)
            setattr(self._cdata[i], "znum", species.znum)
            setattr(self._cdata[i], "mass", species.mass)

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "r", "z", "phi", "charge", "time", "ids", "weight", "vr", "vz",
            "vphi",
            ]:
            self._file.write(field, getattr(self, field))
        self._file.write("anum", np.array(self.species.anum))
        self._file.write("znum", np.array(self.species.znum))

    def export(self) -> dict[str, Species | np.ndarray | unyt.unyt_array]:
        fields = [
            "r", "z", "phi", "charge", "time", "ids", "weight", "vr", "vz",
            "vphi", "species",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            r=self.r, z=self.z, phi=self.phi, charge=self.charge,
            time=self.time, ids=self.ids, vr=self.vr, vphi=self.vphi,
            vz=self.vz, weight=self.weight, species=self.species,
        )

    def unstage(self) -> None:
        super().unstage()
        self._cdata = None


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_particlemarker(
            self,
            species: str | Species,
            charge: utils.Numerical,
            r: utils.Numerical,
            z: utils.Numerical,
            vr: utils.Numerical,
            vphi: utils.Numerical,
            vz: utils.Numerical,
            phi: Optional[utils.Numerical]=None,
            weight: Optional[utils.Numerical]=None,
            time: Optional[utils.Numerical]=None,
            ids: Optional[utils.Numerical]=None,
            note: Optional[str]=None,
            activate: bool=False,
            preview: bool=False,
            save: Optional[bool]=None,
            ) -> ParticleMarker:
        r"""Create marker input in particle (6D) phase-space.

        This is the recommended marker input, especially for simulations that
        trace the whole particle orbit to avoid information loss due to
        guiding-center to particle transformations.

        This input can be used in all simulation modes.

        Parameters
        ----------
        species : str *or* :class:`.Species`
            Marker species.
        charge : array_like (n,) *or* scalar
            Marker charge.
        r : array_like (n,) *or* scalar
            Particle :math:`R` coordinate.
        z : array_like (n,) *or* scalar
            Particle :math:`z` coordinate.
        vr : array_like (n,) *or* scalar
            Velocity :math:`R` component.
        vphi : array_like (n,) *or* scalar
            Velocity :math:`\phi` component.
        vz : array_like (n,) *or* scalar
            Velocity :math:`z` component.
        phi : array_like (n,) *or* scalar, *optional*
            Particle :math:`\phi` coordinate.

            If not given, the values are randomly sampled from an uniform
            distribution between 0 and 360 degrees.
        weight : array_like (n,) *or* scalar, *optional*
            How many physical particles each marker represents.

            Default value is one.
        time : array_like (n,) *or* scalar, *optional*
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
        inputdata : ParticleMarker
            Input variant created from the given parameters.
        """
        if isinstance(species, str):
            species = Species.from_string(species)

        n = 1
        for param in [
            r, z, phi, vr, vphi, vz, weight, time, ids, charge,
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

        r, z, phi, vr, vphi, charge, weight, time, vz = (
            to_array(r), to_array(z), to_array(phi), to_array(vr),
            to_array(vphi), to_array(charge), to_array(weight), to_array(time),
            to_array(vz),
            )

        with utils.validate_variables() as v:
            r = v.validate("r", r, (n,), "m")
            z = v.validate("z", z, (n,), "m")
            vr = v.validate("vr", vr, (n,), "m/s")
            vz = v.validate("vz", vz, (n,), "m/s")
            vphi = v.validate("vphi", vphi, (n,), "m/s")
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

        leaf = ParticleMarker(note=note)
        leaf._stage(
            r=r, z=z, phi=phi, vr=vr, vphi=vphi, vz=vz, charge=charge,
            time=time, ids=cast(np.ndarray, ids), weight=weight,
            species=species,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="marker",
            )
        return leaf
