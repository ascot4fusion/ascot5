"""Interface between IMAS IDS and ASCOT5.

The purpose of this module is to provide a set of tools to read from IDS any
data that is needed to run ASCOT5, and write any output to IDS that ASCOT5 can
supply.

Note that reading from and writing to actual IDS file is separated from reading
and filling the IDS object.

This interface should work with access layers versions from AL3 to AL5.

Notes
-----
- `Python Access Layer documentation \
<https://sharepoint.iter.org/departments/POP/CM/IMDesign/Code%20Documentation/\
ACCESS-LAYER-doc/python/5.4/index.html#/L>`_
- `IMAS Data Dictionary documentation \
<https://sharepoint.iter.org/departments/POP/CM/IMDesign/Data%20Model/CI/\
imas-3.42.0/html_documentation.html>`_
"""
import warnings

import unyt
import numpy as np

import a5py.physlib as physlib
from a5py.ascot5io.nbi import Injector

try:
    import imas
except ImportError:
    imas = None

_TIMEINDEX = 0
"""Index of the time slice to be read from the IDS.

Slicing is done at the moment when the IDS is read, so the time slice is always
the first one.
"""

_STRUCTURED_SPACES = 10
"""Grid identifier in equilibrium IDS corresponding to structured spaces.

This is needed to read stellarator data.
"""

_SPECIES_TYPE = {
    0:"unspecified",
    1:"electron",
    2:"ion",
    3:"ion_state",
    4:"neutral",
    5:"neutral_state",
    6:"neutron",
    7:"photon",
}
"""Type of the species in the distribution sources IDS."""

_GYRO_TYPE = {
    1:"actual",
    2:"gyrocentre",
}
"""Flag indicating whether markers are particles or guiding centers in
   the distribution sources IDS."""

_COORDINATE_IDENTIFIERS = {
    "unspecified":(0, "1"),
    "x":(1, "m"),
    "y":(2, "m"),
    "z":(3, "m"),
    "r":(4, "m"),
    "phi":(5, "rad"),
    "psi":(10, "T*m^2"),
    "rho_tor":(11, "m"), # sqrt( (Phi-Phi_axis) / (pi*B0) )
    "rho_tor_norm":(12, "1"), # sqrt( (Phi-Phi_axis) / (Phi_lcf-Phi_axis) )
    "rho_pol":(13, "m"),
    "rho_pol_norm":(14, "m"),
    "theta":(20, "rad"),
    "theta_straight":(21, "rad"),
    "theta_equal_arc":(22, "rad"),
    "velocity":(100, "m/s"),
    "vx":(101, "m/s"),
    "vy":(102, "m/s"),
    "vz":(103, "m/s"),
    "velocity_phi":(104, "m/s"),
    "velocity_parallel":(105, "m/s"),
    "velocity_perpendicular":(106, "m/s"),
    "velocity_thermal":(107, "m/s"), # normalized to local thermal velocity
    "momentum":(200, "N*s"),
    "momentum_parallel":(201, "N*s"),
    "momentum_perpendicular":(202, "N*s"),
    "canonical_momemtum_phi":(203, "N*m*s"),
    "energy_hamiltonian":(300, "eV"),
    "energy_kinetic":(301, "eV"),
    "magnetic_moment":(302, "J/T"),
    "lambda":(400,"1/T"), # lambda = magnetic_moment / energy_hamiltonian
    "pitch_angle":(401, "1"), # angle between velocity and magnetic field
    "pitch":(403, "1"),
    "pitch_at_min_b":(404, "1"),
}
"""Coordinates and corresponding identifiers and units in coordinate arrays.

These are specified in utilities/coordinate_identifier.xml in the Data
Dictionary repository.

Notes
-----

- IMAS specifications state that `rho_tor_norm` and `rho_pol_norm` are in units
  of 'm' but this can't be correct.
- Likewise, `velocity_thermal` is in units of 'm/s' but this also a normalized
  quantity.
- `pitch_angle` has no units in the specification.
"""

_MOMENTS = {
    "density":"density_fast",
    "toroidalcurrent":"current_fast_phi",
    "pressure":"pressure_fast",
    "electronpowerdep":"collisions.electrons.powerthermal",
}
"""Mapping of ASCOT5 moment names to IMAS distribution names."""

def read_ids(
        ids_name, query, backend="hdf5", occurrence=0, time_requested=0,
        host=None,
        ):
    """Read contents of an IDS.

    Parameters
    ----------
    ids_name : str
        Name of the IDS database to open, e.g. "equilibrium", "core_profiles",
        etc.
    query : str or dict[str,str]
        Either path or legacy keys to specify where the IDS data is located.

        Path (only available for AL5) can either be absolute or relative path to
        the IDS file, e.g.:

        query="/absolute/path/to/data"

        The legacy keys are: "user", "tokamak", "version", "shot"
        ("pulse" in AL5), and "run", e.g.:

        query={
            "user":"public",
            "shot":131024,
            "run":41,
            "database":"ITER",
            "version":3,
            }
    backend : {"hdf5", "mdsplus", "memory", "ascii"}, optional
        The backend in which the IDS is stored.
    occurrence : int, optional
        If there are multiple instances of the same dataset, which one to use.
    time_requested: float, optional
        The time slice at which the data is read.
    host : str, optional
        The host where the IDS database is located if it is not local.

    Returns
    -------
    ids: imas.ids_base.IDSBase
        The IDS database object which was read.
    """
    path = "imas:"
    if host:
        path += f"[//{host}/]"
    path += f"{backend}?"

    if isinstance(query, str):
        path += f"path={query}"
    else:
        for key, value in query.items():
            if key == "shot":
                path += f"pulse={value};"
            elif key in ["user", "tokamak", "version", "pulse", "run"]:
                path += f"{key}={value};"
            else:
                raise ValueError(f"Unknown query key: {key}")
    try: # AL5
        with imas.DBEntry(path, mode="r") as imas_entry:
            try: # IMAS-Py
                interp = imas.ids_defs.CLOSEST_INTERP
            except: # Python access layer
                interp = imas.imasdef.CLOSEST_SAMPLE
            ids = imas_entry.get_slice(
                ids_name, time_requested, interp, autoconvert=False,
                )
    except AttributeError: # Legacy AL3/4
        if isinstance(query, str):
            raise ValueError("`query` must use the legacy format for AL < 5")
        ids = imas.ids(query["shot"], query["run"])
        ids.open_env(query["user"], query["tokamak"], query["version"])

        idsdata = ids.__dict__[ids_name]
        if "get" not in dir(idsdata):
            idsdata = ids.__dict__[ids_name + 'Array']
            idsdata.get(occurrence)
        else:
            idsdata.get(occurrence)
    return ids

def write_ids(ids, query, backend="hdf5", occurrence=0, host=None):
    """Write (pre-filled) IDS object to a file.

    Parameters
    ----------
    ids : imas.ids_base.IDSBase
        The IDS object to be written to the file.
    query : str or dict[str,str]
        Either path or legacy keys to specify where the IDS data is located.

        Path (only available for AL5) can either be absolute or relative path to
        the IDS file, e.g.:

        query="/absolute/path/to/data"

        The legacy keys are: "user", "tokamak", "version", "shot"
        ("pulse" in AL5), and "run", e.g.:

        query={
            "user":"public",
            "shot":131024,
            "run":41,
            "database":"ITER",
            "version":3,
            }
    backend : {"hdf5", "mdsplus", "memory", "ascii"}, optional
        The backend in which the IDS is stored.
    occurrence : int, optional
        The occurrence of the created data set.
    host : str, optional
        The host where the IDS database is located if it is not local.
    """
    path = "imas:"
    if host:
        path += f"[//{host}/]"
    path += f"{backend}?"

    if isinstance(query, str):
        path += f"path={query}"
    else:
        for key, value in query.items():
            if key == "shot":
                path += f"pulse={value};"
            elif key in ["user", "tokamak", "version", "pulse", "run"]:
                path += f"{key}={value};"
            else:
                raise ValueError(f"Unknown query key: {key}")
    try: # AL5
        with imas.DBEntry(path, mode="w") as imas_entry:
            imas_entry.put(ids, occurrence=occurrence)
    except AttributeError: # Legacy AL3/4
        if isinstance(query, str):
            raise ValueError("`query` must use the legacy format for AL < 5")
        ids = imas.ids(query["shot"], query["run"])
        ids.open_env(query["user"], query["tokamak"], query["version"])
        ids.put(ids, occurrence=occurrence)

class ImportImas():
    """Class for importing data from IMAS IDS."""

    def imas_b2ds(self, equilibrium_ids=None):
        """Import axisymmetric tokamak magnetic field from IMAS IDS.

        Parameters
        ----------
        equilibrium_ids : imas.ids_base.IDSBase
            The IDS database object which contains the equilibrium data.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        def fix_cocos(*psivals):
            """Convert psi to correct COCOS."""
            psivals_new = [None] * len(psivals)
            for i, psi in enumerate(psivals):
                psivals_new[i] =  psi / ( 2*np.pi )
            return psivals_new

        def identify_rzgrid(profiles_2d):
            """Identify the rectangular (R,z) grid in the profiles_2d."""
            RZGRID_INDEX = 1
            for index, profile in enumerate(profiles_2d):
                if  profile.grid_type.index == RZGRID_INDEX:
                    psi = profile.psi
                    rzgrid = profile.grid
                    try:
                        bphi = profile.b_field_tor
                    except ValueError:
                        bphi = profile.b_field_phi
                    return psi, bphi, rzgrid
            raise IOError(
                "No rectangular (R,z) grid found in `equilibrium` IDS "
                "(profiles_2d)"
            )

        profiles_2d, global_quantities = (
            equilibrium_ids.time_slice[_TIMEINDEX].profiles_2d,
            equilibrium_ids.time_slice[_TIMEINDEX].global_quantities,
        )

        psi, bphi, rzgrid = identify_rzgrid(profiles_2d)
        nr, rmin, rmax = rzgrid.dim1.size, rzgrid.dim1[0], rzgrid.dim1[-1]
        nz, zmin, zmax = rzgrid.dim2.size, rzgrid.dim2[0], rzgrid.dim2[-1]

        br = np.zeros_like(bphi)
        bz = np.zeros_like(bphi)
        psi0, psi1, axisr, axisz = (
            global_quantities.psi_axis,
            global_quantities.psi_boundary,
            global_quantities.magnetic_axis.r,
            global_quantities.magnetic_axis.z,
            )
        psi1, psi0, psi = fix_cocos(psi1, psi0, psi)

        b2d = {
            "nr":nr,
            "nz":nz,
            "br":br,
            "bz":bz,
            "psi":psi,
            "psi0":psi0,
            "psi1":psi1,
            "bphi":bphi,
            "rmin":rmin,
            "rmax":rmax,
            "zmin":zmin,
            "zmax":zmax,
            "axisr":axisr,
            "axisz":axisz,
        }
        return "B_2DS", b2d

    def imas_bsts(self, equilibrium_ids=None, ggd_index=0):
        """Import stellarator magnetic field from IMAS IDS.

        Parameters
        ----------
        equilibrium_ids : imas.ids_base.IDSBase
            The IDS database object which contains the equilibrium data.
        ggd_index : int, optional
            The index in the `description_ggd` array in the IDS from which the
            data is read.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """

        psi0, psi1 = 0., 1.
        warnings.warn("Making a bold assumption that psi0 = 0.0 and psi1 = 1.0")
        if len(equilibrium_ids.grids_ggd[_TIMEINDEX].grid) < 2:
            raise ValueError(
                "The IDS must contain both (R,phi,z) grid and "
                "the magnetic axis grid."
            )

        # By convention, index 0 is the (R,phi,z) grid.
        grid = equilibrium_ids.grids_ggd[_TIMEINDEX].grid[0]
        ggd = equilibrium_ids.time_slice[_TIMEINDEX].ggd[ggd_index]
        space = grid.space

        if grid.identifier.index != _STRUCTURED_SPACES:
            raise ValueError(
                "Expected grid identifier '{_STRUCTURED_SPACES}' "
                "(structured spaces) but got '{grid.identifier.index}' instead."
                )
        if (len(space) == 3 and
            space[0].coordinates_type[0] == _COORDINATE_IDENTIFIERS["r"][0] and
            space[1].coordinates_type[0] == _COORDINATE_IDENTIFIERS["phi"][0] and
            space[2].coordinates_type[0] == _COORDINATE_IDENTIFIERS["z"][0]
        ):
            raise ValueError("The grid is not (R,phi,z).")
        nr, nphi, nz = (
            len(grid.space[0].objects_per_dimension[0].object),
            len(grid.space[1].objects_per_dimension[0].object),
            len(grid.space[2].objects_per_dimension[0].object),
        )
        if not all((nr, nphi, nz)):
            raise ValueError("The grid has null dimension(s).")

        r = np.array([
            space[0].objects_per_dimension[0].object[i].geometry[0]
            for i in range(nr)
            ]) * unyt.unyt_quantity.from_string(_COORDINATE_IDENTIFIERS["r"][1])
        z = np.array([
            space[2].objects_per_dimension[0].object[i].geometry[0]
            for i in range(nz)
            ]) * unyt.unyt_quantity.from_string(_COORDINATE_IDENTIFIERS["z"][1])
        phi = np.array([
            space[1].objects_per_dimension[0].object[i].geometry[0]
            for i in range(nphi)
            ]) * unyt.unyt_quantity.from_string(_COORDINATE_IDENTIFIERS["phi"][1])

        order, shape = "F", (nr, nphi, nz)
        br, bz, bphi, psi = (
            np.reshape(ggd.b_field_r[0].values, newshape=shape, order=order),
            np.reshape(ggd.b_field_tor[0].values, newshape=shape, order=order),
            np.reshape(ggd.b_field_z[0].values, newshape=shape, order=order),
            np.reshape(ggd.psi[0].values, newshape=shape, order=order),
        )

        # By convention, index 1 is the magnetic axis.
        axisgrid = equilibrium_ids.grids_ggd[_TIMEINDEX].grid[1]
        space = axisgrid.space
        axis_nphi = len(space[0].objects_per_dimension[0].object)
        if not axis_nphi:
            raise ValueError("The axis grid has null dimension.")
        if (len(space) == 3 and
            space[0].coordinates_type[0] == _COORDINATE_IDENTIFIERS["r"][0] and
            space[1].coordinates_type[0] == _COORDINATE_IDENTIFIERS["phi"][0] and
            space[2].coordinates_type[0] == _COORDINATE_IDENTIFIERS["z"][0]
        ):
            raise ValueError("The axis grid is not (R,phi,z).")

        axis_r = np.array([
            space[0].objects_per_dimension[0].object[i].geometry[0]
            for i in range(nr)
            ]) * unyt.unyt_quantity.from_string(_COORDINATE_IDENTIFIERS["r"][1])
        axis_z = np.array([
            space[2].objects_per_dimension[0].object[i].geometry[0]
            for i in range(nz)
            ]) * unyt.unyt_quantity.from_string(_COORDINATE_IDENTIFIERS["z"][1])
        axis_phi = np.array([
            space[1].objects_per_dimension[0].object[i].geometry[0]
            for i in range(nphi)
            ]) * unyt.unyt_quantity.from_string(_COORDINATE_IDENTIFIERS["phi"][1])

        b = {
            "b_rmin":r[0],
            "b_rmax":r[-1],
            "b_nr":nr,
            "b_phimin":phi[0],
            "b_phimax":phi[-1],
            "b_nphi":nphi,
            "b_zmin":z[0],
            "b_zmax":z[-1],
            "b_nz":nz,
            "psi_rmin":r[0],
            "psi_rmax":r[-1],
            "psi_nr":nr,
            "psi_phimin":phi[0],
            "psi_phimax":phi[-1],
            "psi_nphi":nphi,
            "psi_zmin":z[0],
            "psi_zmax":z[-1],
            "psi_nz":nz,
            "axis_phimin":axis_phi[0],
            "axis_phimax":axis_phi[-1],
            "axis_nphi":axis_nphi,
            "axisr":axis_r,
            "axisz":axis_z,
            "br":br,
            "bphi":bphi,
            "bz":bz,
            "psi":psi,
            "psi0":psi0,
            "psi1":psi1,
        }
        return "B_STS", b

    def imas_wall2d(self, wall_ids=None, unit=0):
        """Import wall contour from IMAS IDS.

        Parameters
        ----------
        wall_ids : imas.ids_base.IDSBase
            The IDS database object which contains the wall data.
        unit : int, optional
            The index of the wall unit to be read.

            Note that the limiter may be made up of several wall units. This is
            not currently supported in ASCOT5.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        outline = wall_ids.description_2d[_TIMEINDEX].limiter.unit[unit].outline
        w = {"r":outline.r, "z":outline.z, "nelements":outline.r.size}
        return "wall_2D", w

    def imas_wall3d(self, wall_ids=None, ggd_index=0, space_index=0):
        """Import wall triangles from IMAS IDS.

        Parameters
        ----------
        wall_ids : imas.ids_base.IDSBase
            The IDS database object which contains the wall data.
        ggd_index : int, optional
            The index in the `description_ggd` array in the IDS from which the
            data is read.
        space_index : int, optional
            The index in the `grid_ggd.space` array in the IDS from which the
            data is read.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        desc_ggd = wall_ids.description_ggd[ggd_index]
        grid_ggd = desc_ggd.grid_ggd[_TIMEINDEX]
        space = grid_ggd.space[space_index]

        valid_coordinates = [
            _COORDINATE_IDENTIFIERS["x"][0],
            _COORDINATE_IDENTIFIERS["y"][0],
            _COORDINATE_IDENTIFIERS["z"][0],
            ]
        if not np.isin(space.coordinates_type, valid_coordinates).all():
            raise ValueError(
                "The wall elements must be given in (R, z, phi) coordinates."
            )

        nnodes = len(space.objects_per_dimension[0].object)
        node_coordinates = np.zeros((3, nnodes), dtype=float)
        for i in range(nnodes):
            node_coordinates[:,i] = space.objects_per_dimension[0].object[i].geometry

        nelements = len(space.objects_per_dimension[2].object)
        xyz = np.zeros((nelements, 3, 3), dtype=float)
        for itri in range( nelements ):
            triangle_indexes = space.objects_per_dimension[2].object[itri].nodes
            for i in range(3):
                # The element indexes in IMAS start from 1, not 0.
                xyz[itri,i,:] = node_coordinates[:,triangle_indexes[i]-1]
        label = (
            f"Wall IDS. ggd index: {ggd_index}, "
            "time index: {_TIMEINDEX}, "
            "space index: {space_index}."
        )
        w = {
            "flag":np.ones((nelements, 1)),
            "labels":{label:1},
            "x1x2x3":xyz[:,:,0],
            "y1y2y3":xyz[:,:,1],
            "z1z2z3":xyz[:,:,2],
            "nelements":[[nelements]],
        }
        return "wall_3D", w

    def imas_plasma(
            self, equilibrium_ids=None, core_profiles_ids=None, psi0=None,
            psi1=None
            ):
        """Import plasma profiles from IMAS IDS.

        The data in IDS is tabulated with respect to rho_tor, which is why
        the equilibrium IDS is needed to convert rho_tor to rho_pol.

        Parameters
        ----------
        equilibrium_ids : imas.ids_base.IDSBase
            The IDS database object which contains the equilibrium data.
        core_profiles_ids : imas.ids_base.IDSBase
            The IDS database object which contains the core profiles.
        psi0 : float, optional
            The value of psi at the magnetic axis.

            Might be needed to normalize psi grid.
        psi1 : float, optional
            The value of psi at the plasma boundary.

            Might be needed to normalize psi grid.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        iElement  = 0
        profiles_1d = core_profiles_ids.profiles_1d[_TIMEINDEX]
        grid = profiles_1d.grid
        eqprof = equilibrium_ids.time_slice[_TIMEINDEX].profiles_1d
        if (grid.rho_pol_norm):
            rho = grid.rho_pol_norm
        elif len(grid.psi) and grid.psi_magnetic_axis != grid.psi_boundary:
            rho = np.sqrt(
                (grid.psi - grid.psi_magnetic_axis) /
                (grid.psi_boundary - grid.psi_magnetic_axis)
                )
        elif len(grid.rho_tor) and len(eqprof.rho_tor) and len(eqprof.psi_norm):
            rho = np.interp(grid.rho_tor, eqprof.rho_tor, eqprof.psi_norm)
        elif len(grid.rho_tor) and len(eqprof.rho_tor) and len(eqprof.psi):
            psi = np.interp(grid.rho_tor, eqprof.rho_tor, eqprof.psi) / ( 2*np.pi )
            rho = np.sqrt( (psi - psi0) / (psi1 - psi0) )
        else:
            raise ValueError(
                "No sufficient data to initialize radial grid."
            )

        nrho, nion = rho.size, len(profiles_1d.ion)
        anum = np.array([
            int(profiles_1d.ion[i].element[iElement].a) for i in range(nion)
        ])
        znum = np.array([
            profiles_1d.ion[i].element[iElement].z_n for i in range(nion)
        ])
        mass = np.array([
            physlib.species.autodetect(a, z)["mass"] for a, z in zip(anum, znum)
        ]) / unyt.amu
        warnings.warn(
            "Assuming fully ionized ions in IMAS import. "
            "This may not be correct for all cases."
            )
        charge = znum * unyt.e

        edensity = profiles_1d.electrons.density #/ unyt.m**3
        idensity = np.zeros(shape=(nrho,nion)) #/ unyt.m**3
        for i in range(nion):
            idensity[:,i] = profiles_1d.ion[i].density_thermal

        etemperature = profiles_1d.electrons.temperature
        itemperature = profiles_1d.ion[0].temperature
        species_averaged_temperature_exists = len(profiles_1d.t_i_average) > 0
        if species_averaged_temperature_exists:
            itemperature = profiles_1d.t_i_average

        ivtor = np.zeros_like(idensity) * unyt.rad/unyt.s
        if len(profiles_1d.ion[i].rotation_frequency_tor) > 0:
            ivtor[:,i] = (
                profiles_1d.ion[i].rotation_frequency_tor
            )
        vtor = ivtor.mean(axis=1)
        valid_ions = idensity[0, :] > 10
        idensity = idensity[:, valid_ions]
        nion = idensity.shape[1]
        pls = {
            "rho":rho,
            "nrho":nrho,
            "nion":nion,
            "vtor":vtor,
            "anum":anum[valid_ions],
            "znum":znum[valid_ions],
            "mass":mass[valid_ions],
            "charge":charge[valid_ions],
            "idensity":idensity,
            "edensity":edensity,
            "itemperature":itemperature,
            "etemperature":etemperature,
        }

        pls = self._ascot.data.create_input(
            "import plasma profiles", dryrun=True, pls=pls, extrapolate=3.0
        )
        return "plasma_1D", pls

    def imas_marker(self, distribution_sources_ids=None):
        """Import markers from IMAS IDS.

        Parameters
        ----------
        distribution_sources_ids : imas.ids_base.IDSBase
            The IDS database object which contains the kinetic particle sources.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        mrk_all = {}
        time = distribution_sources_ids.time[_TIMEINDEX]
        for source in distribution_sources_ids.source:
            markers = source.markers[_TIMEINDEX]
            if len(markers.weights) == 0:
                continue
            if(not _SPECIES_TYPE[int(source.species.type.index)]
               in ["ion", "ion_state"]):
                warnings.warn(
                    "Skipping a source with species type: "
                    f"{_SPECIES_TYPE[source.species.type.index]} "
                    "(only ions are supported)."
                )
                continue

            n = len(markers.weights)
            mrk = {}
            indices = np.array([x.index for x in markers.coordinate_identifier ])
            if _GYRO_TYPE[int(source.gyro_type)] == "actual":
                for field in ["r", "z", "phi", "vx", "vy", "vz"]:
                    units = unyt.unyt_quantity.from_string(
                        _COORDINATE_IDENTIFIERS[field][1]
                        )
                    idx = np.argwhere(
                        indices == _COORDINATE_IDENTIFIERS[field][0]
                        )[0][0]
                    mrk[field] = markers.positions[:,idx] * units
                x, y = physlib.pol2cart(mrk["r"], mrk["phi"])
                mrk["vr"], mrk["vphi"] = (
                    physlib.cart2pol_vec(mrk["vx"], x, mrk["vy"], y)
                )
                del mrk["vx"], mrk["vy"]
            elif _GYRO_TYPE[int(source.gyro_type)] == "gyrocentre":
                for field in ["r", "z", "phi", "energy_kinetic", "pitch"]:
                    units = unyt.unyt_quantity.from_string(
                        _COORDINATE_IDENTIFIERS[field][1]
                        )
                    idx = np.argwhere(
                        indices == _COORDINATE_IDENTIFIERS[field][0]
                        )[0][0]
                    if field == "energy_kinetic":
                        mrk["energy"] = markers.positions[:,idx] * units
                    else:
                        mrk[field] = markers.positions[:,idx] * units
            else:
                raise ValueError(f"Unknown gyrotype = '{source.gyro_type}'")

            anum = np.rint(source.species.ion.element[0].a)
            znum = np.rint(source.species.ion.element[0].z_n)
            mass = physlib.species.autodetect(anum, znum)["mass"]
            charge = physlib.species.autodetect(anum, znum)["charge"]

            mrk.update({
                "n":n,
                "ids":np.arange(n) + 1,
                "time":time * np.ones((n,), dtype="f8") * unyt.s,
                "anum":anum * np.ones((n,), dtype="i8"),
                "znum":znum * np.ones((n,), dtype="i8"),
                "mass":mass * np.ones((n,), dtype="f8"),
                "charge":charge * np.ones((n,), dtype="f8"),
                "weight":np.array(markers.weights) * unyt.particles / unyt.s,
            })
            for key, value in mrk.items():
                if not key in mrk_all:
                    mrk_all[key] = value
                elif key == "n":
                    mrk_all[key] += value
                elif key == "ids":
                    ids = mrk_all[key][-1]
                    mrk_all[key] = np.concatenate((mrk_all[key], ids + value))
                else:
                    mrk_all[key] = np.concatenate((mrk_all[key], value))
            del mrk

        mrk_all["phi"] = mrk_all["phi"].to("deg")
        mrk_type = "Prt" if _GYRO_TYPE[int(source.gyro_type)] == "actual" else "GC"
        return mrk_type, mrk_all

    def imas_nbi(self, nbi_ids=None):
        """Import neutral beam injectors from IMAS IDS.

        Parameters
        ----------
        nbi_ids : imas.ids_base.IDSBase
            The IDS database object which contains the NBI specifications.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        nunit = len(nbi_ids.unit)
        injs = []
        for unit, injector in enumerate(nbi_ids.unit):
            inj = {
                "anum":injector.species.a,
                "znum":injector.species.z_n,
                "power":injector.power_launched.data,
                "energy":injector.energy.data,
                "efrac":injector.beam_current_fraction.data.ravel(),
                "powerfrac":injector.beam_power_fraction.data.ravel(),
            }
            inj["mass"] = physlib.species.autodetect(
                inj["anum"], inj["znum"])["mass"]

            divergence = injector.beamlets_group[0].divergence_component
            f1, f2 = (
                divergence[0].particles_fraction,
                divergence[1].particles_fraction
                )
            # Ascot defines divergence as 1/e folding in power whereas in IDS it
            # is 1/e folding in intensity. The difference is sqrt(2).
            fix_divergence = np.sqrt(2)
            inj.update({
                "divh":divergence[0].horizontal*fix_divergence,
                "divv":divergence[0].vertical*fix_divergence,
                "halofraction":f2 / ( f1 + f2 ),
                "halodivh":divergence[1].horizontal*fix_divergence,
                "halodivv":divergence[1].vertical*fix_divergence,
            })

            r, phi, z, tangencyradii, angles = (
                np.array([]), np.array([]), np.array([]), np.array([]),
                np.array([])
            )
            direction = injector.beamlets_group[0].direction
            for beamletgroup in injector.beamlets_group:
                r, z, phi, tangencyradii, angles = (
                    np.append(r, beamletgroup.beamlets.positions.r),
                    np.append(z, beamletgroup.beamlets.positions.z),
                    np.append(phi, beamletgroup.beamlets.positions.phi),
                    np.append(tangencyradii, beamletgroup.beamlets.tangency_radii),
                    np.append(angles, beamletgroup.beamlets.angles),
                )

            x, y = physlib.pol2cart(r, phi)
            beta = np.arccos(tangencyradii / r)
            alpha = np.pi + phi - direction * ( 0.5 * np.pi - beta )

            # Fixing origo to the beamlet position, we can calculate the beam
            # center line vector from spherical coordinates
            dxyz = physlib.sph2cart(
                r=1.0,
                phi=alpha,
                theta=np.pi/2 - angles,
            )

            injs.append(
                Injector(
                    unit+1,
                    inj["anum"],
                    inj["znum"],
                    inj["mass"],
                    inj["energy"],
                    inj["efrac"],
                    inj["power"],
                    inj["divh"],
                    inj["divv"],
                    inj["halofraction"],
                    inj["halodivh"],
                    inj["halodivv"],
                    x.size,
                    x,
                    y,
                    z,
                    dxyz[0],
                    dxyz[1],
                    dxyz[2]
                )
            )
        return ("nbi", {"ninj":nunit, "injectors":injs})

class ExportIMAS():
    """Class for exporting data to IMAS IDS.

    This class will be inherited by the run object, which has access to all data
    associated with the simulation.
    """

    def imas_distributions(self, distributions_ids=None):
        """Export distributions to IMAS distributions IDS.

        Parameters
        ----------
        distributions_ids : imas.ids_base.IDSBase, optional
            The distributions IDS which will be filled here.

            If not given, a new one is created.

        Returns
        -------
        distributions_ids : imas.ids_base.IDSBase
            The filled distributions IDS object.

            Same object as the input parameter `distributions_ids` if it was
            given.
        """
        SPECIES_TYPE_INV = {v: k for k, v in _SPECIES_TYPE.items()}
        if distributions_ids is None:
            distributions_ids = imas.distributions_ids()
            distributions_ids.ids_properties.homogeneous_time = (
                imas.imasdef.IDS_TIME_MODE_HOMOGENEOUS
            )
        codeversion = self.getcodeversion()
        for field in ["name", "description", "commit", "version", "repository"]:
            setattr(distributions_ids.code, field, codeversion[field])

        distributions_ids.code.parameters = self.options.toxml()
        distributions_ids.code.output_flag = np.array([self.get_run_success()])

        species = self.getspecies()
        anum, znum, mass = (
            float(species["anum"]),
            float(species["znum"]),
            float(species["mass"]),
        )

        gyro_type = self.getsimmode()
        if not (gyro_type == 1 or gyro_type == 2) :
            raise ValueError(
                "Only pure gyro or guiding center simulations can be exported "
                "to IMAS."
                )

        if "5d" in self.getdist_list()[0]:
            d5d = self.getdist('5d')
            charges = d5d.abscissa('charge')
        else:
            charges = []

        dist_index = 0 # Distributions are appended after each other
        distributions_ids.distribution.resize(len(charges))
        for charge in charges:
            d = distributions_ids.distribution[dist_index]
            d.is_delta_f = int(True)
            d.gyro_type = gyro_type

            if anum == 0:
                species.type.index, species.type.name, species.type.description = (
                    SPECIES_TYPE_INV["electron"], "electron", "Electron"
                )
            elif round(charge) == 0:
                species.type.index, species.type.name, species.type.description = (
                    SPECIES_TYPE_INV["neutral"], "neutral",
                    "Neutral species in a single/average state; refer to "
                    "neutral-structure"
                )
            else:
                species.type.index, species.type.name, species.type.description = (
                    SPECIES_TYPE_INV["ion"], "ion",
                    "Ion species in a single/average state; refer to "
                    "ion-structure"
                )

            if _SPECIES_TYPE[species.type.index] in ["electron", "ion"]:
                element = species.ion.element
            else:
                element = species.neutral.element
            element.resize(1) # Note that we assume single nucleus species
            element[0].a, element[0].z_n, element[0].atoms_n = (
                float(mass), float(znum), 1
            )
            species.z_ion = float(charge)

            prof2d = d.profiles_2d.resize(_TIMEINDEX+1)
            prof2d = d.profiles_2d[_TIMEINDEX]
            _type = prof2d.grid.type
            _type.index, _type.name, _type.description = (
                0,
                "rectangular RZ",
                "Rectangular grid in the (R,Z) coordinates;",
                )
            prof2d.grid.r, prof2d.grid.z = (
                d5d.abscissa("r").v, d5d.abscissa("z").v
            )

            moment = self.getdist_moments(self.getdist("5d"), *_MOMENTS.keys())
            for ascotname, imasname in _MOMENTS.items():
                ordinate = moment.ordinate(ascotname, toravg=True)
                setattr(prof2d, imasname, np.array(ordinate))

            warningtext = (
                "5D distribution output still WIP; fill in e.g. profiles_2d:\n"
                "density --> density_fast\n"
                "toroidalcurrent --> current_fast_phi\n"
                "pressure --> pressure_fast\n"
                "electronpowerdep --> collisions.electrons.powerthermal\n"
                "\n"
                "The following may not be useful:\n"
                "parallelcurrent : Parallel current\n"
                "chargedensity : Charge density\n"
                "energydensity : Energy density\n"
                "powerdep : Total deposited power\n"
                "ionpowerdep : Power deposited to ions (should be per species)\n"
            )
            warnings.warn(warningtext)

        warnings.warn("Rho distribution output not implemented")
        #for iDistribution in range(len(charges)):
        #    d = self.ids.distribution[iDistribution + distoffset]
        #distoffset += len(charges)
        return distributions_ids
