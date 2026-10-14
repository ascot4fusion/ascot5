"""Defines particle (gyro-orbit) input class and the corresponding factory
method.
"""
from typing import Optional, Any, Sized, cast

import unyt
import numpy as np

from a5py import utils
from a5py.physlib import Species
from a5py.data.access import InputVariant, Leaf, TreeMixin
from a5py.physlib.formulas import resolve_quantities, compute_quantities
from a5py.engine.interpolate import interpolate

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
    def mass(self) -> unyt.unyt_array:
        """Marker species mass."""
        return self.species.mass * np.ones((self.n,))

    @property
    def anum(self) -> unyt.unyt_array:
        """Marker species atomic mass number."""
        return self.species.anum * np.ones((self.n,))

    @property
    def znum(self) -> unyt.unyt_array:
        """Marker species charge number."""
        return self.species.znum * np.ones((self.n,))

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
                out[i] = self._cdata[i].rprt
            return out
        assert self._file is not None
        return self._file.read("r")

    @property
    def phi(self) -> unyt.unyt_array:
        r"""Particle :math:`\phi` coordinate."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "rad")
            for i in range(out.size):
                out[i] = self._cdata[i].phiprt

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
                out[i] = self._cdata[i].zprt
            return out
        assert self._file is not None
        return self._file.read("z")

    @property
    def pr(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`R` component."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pr
            return out
        assert self._file is not None
        return self._file.read("pr")

    @property
    def pphi(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`\phi` component."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pphi
            return out
        assert self._file is not None
        return self._file.read("pphi")

    @property
    def pz(self) -> unyt.unyt_array:
        r"""Particle momentum :math:`z` component."""
        if self._cdata is not None:
            out = unyt.unyt_array([0.]*self.n, "kg*m/s")
            for i in range(out.size):
                out[i] = self._cdata[i].pz
            return out
        assert self._file is not None
        return self._file.read("pz")

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
            pr: unyt.unyt_array,
            pphi: unyt.unyt_array,
            pz: unyt.unyt_array,
            weight: unyt.unyt_array,
            species: Species,
            ) -> None:
        self._cdata = (Structure * ids.size)()
        for i in range(ids.size):
            setattr(self._cdata[i], "rprt", r[i])
            setattr(self._cdata[i], "zprt", z[i])
            setattr(self._cdata[i], "id", ids[i])
            setattr(self._cdata[i], "phiprt", phi[i].to("rad"))
            setattr(self._cdata[i], "time", time[i])
            setattr(self._cdata[i], "pr", pr[i])
            setattr(self._cdata[i], "pphi", pphi[i])
            setattr(self._cdata[i], "pz", pz[i])
            setattr(self._cdata[i], "charge", int(charge[i].to("e")))
            setattr(self._cdata[i], "weight", weight[i])
            setattr(self._cdata[i], "anum", species.anum)
            setattr(self._cdata[i], "znum", species.znum)
            setattr(self._cdata[i], "mass", species.mass.to("kg"))

    def _save_data(self) -> None:
        assert self._file is not None
        for field in [
            "r", "z", "phi", "charge", "time", "ids", "weight", "pr", "pz",
            "pphi",
            ]:
            self._file.write(field, getattr(self, field))
        self._file.write("anum", np.array(self.species.anum))
        self._file.write("znum", np.array(self.species.znum))

    def export(self) -> dict[str, Species | np.ndarray | unyt.unyt_array]:
        fields = [
            "rprt", "zprt", "phiprt", "charge", "time", "ids", "weight", "pr", "pz",
            "pphi", "species",
            ]
        return {field: getattr(self, field) for field in fields}

    def stage(self) -> None:
        super().stage()
        self._stage(
            r=self.r, z=self.z, phi=self.phi, charge=self.charge,
            time=self.time, ids=self.ids, pr=self.pr, pphi=self.pphi,
            pz=self.pz, weight=self.weight, species=self.species,
        )

    def unstage(self) -> None:
        super().unstage()
        self._cdata = None

    def evaluate(self, *qnt, filter=None, bfield=None, go2gc=False):
        marker_ids = self.ids
        if filter is not None:
            marker_ids = marker_ids[np.isin(marker_ids, filter)]
        sorter = np.argsort(marker_ids)

        marker_quantities = [
            "r", "z", "phi", "time", "charge", "weight", "pr", "pphi",
            "pz", "mass"
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


# pylint: disable=too-few-public-methods
class CreateMixin(TreeMixin):
    """Provides the factory method."""

    #pylint: disable=protected-access, too-many-arguments, too-many-locals
    def create_particlemarker(
            self,
            species: str | Species,
            r: utils.ArrayLike,
            z: utils.ArrayLike,
            pr: utils.ArrayLike,
            pphi: utils.ArrayLike,
            pz: utils.ArrayLike,
            phi: Optional[utils.ArrayLike]=None,
            charge: Optional[utils.ArrayLike]=None,
            weight: Optional[utils.ArrayLike]=None,
            time: Optional[utils.ArrayLike]=None,
            ids: Optional[utils.ArrayLike]=None,
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
        pr : array_like (n,) *or* scalar
            Momentum :math:`R` component.
        pphi : array_like (n,) *or* scalar
            Momentum :math:`\phi` component.
        pz : array_like (n,) *or* scalar
            Momentum :math:`z` component.
        phi : array_like (n,) *or* scalar, *optional*
            Particle :math:`\phi` coordinate.

            If not given, the values are randomly sampled from an uniform
            distribution between 0 and 360 degrees.
        charge : array_like (n,) *or* scalar, *optional*
            Marker charge (state).

            Default is to assume that the marker is fully ionized.
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
            r, z, phi, pr, pphi, pz, weight, time, ids, charge,
            ]:
            if isinstance(param, Sized):
                if utils.size(param) == 1:
                    continue
                if n != 1 and n != utils.size(param):
                    raise ValueError(
                        "Input arrays have inconsistent sizes."
                        )
                n = utils.size(param)

        to_array = lambda x: utils.scalar2array(x, n)

        r, z, phi, pr, pphi, charge, weight, time, pz = (
            to_array(r), to_array(z), to_array(phi), to_array(pr),
            to_array(pphi), to_array(charge), to_array(weight), to_array(time),
            to_array(pz),
            )

        with utils.validate_variables() as v:
            r = v.validate("r", r, (n,), "m")
            z = v.validate("z", z, (n,), "m")
            pr = v.validate("pr", pr, (n,), "kg*m/s")
            pz = v.validate("pz", pz, (n,), "kg*m/s")
            pphi = v.validate("pphi", pphi, (n,), "kg*m/s")

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
            default_charge = species.znum if species.znum != 0 else -1
            charge = v.validate(
                "charge", charge, (n,), "e",
                default=np.full((n,), default_charge),
                )

        leaf = ParticleMarker(note=note)
        leaf._stage(
            r=r, z=z, phi=phi, pr=pr, pphi=pphi, pz=pz, charge=charge,
            time=time, ids=cast(np.ndarray, ids), weight=weight,
            species=species,
            )
        if preview:
            return leaf
        self._treemanager.enter_leaf(
            leaf, activate=activate, save=save, category="marker",
            )
        return leaf
