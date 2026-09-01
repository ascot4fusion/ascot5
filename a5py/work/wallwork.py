import numpy as np

def evaluate_0d(flags=None):
    """Get peak power load and other 0D quantities related to wall loads.

    Parameters
    ----------
    flags : array_like, optional
        Filter output to include only the elements whose flag is in this
        list (the values can either be integers or strings).

    Returns
    -------
    warea : float
        Total wetted area.
    pload : float
        Peak power load.
    """
    _, area, edepo, _, _ = self.getwall_loads(flags=flags)
    wetted_total = np.sum(area)
    energy_peak  = np.amax(edepo/area)
    return wetted_total, energy_peak

def evaluate_loads(run, weights=True, p_ids=None, flags=None):
    """Get wall loads and associated quantities.

    This method does not return loads on all wall elements (as usually most
    of them receive no loads) but only those that are affected and their
    IDs.

    For 2D walls, iangle, angle of incidence is not calculated.

    Parameters
    ----------
    weights : bool, optional
        Include marker weights to get physical results (otherwise particle
        deposition would be just the number of markers that hit the tile).

        Dropping weights is useful to check how many markers hit a tile
        which tells us how good the statistics are.
    p_ids : array_like, optional
        Calculate wall loads only for the particles with the given IDs.
    flags : str, int, or array_like, optional
        Filter output to include only the elements whose flag is in this
        list (the values can either be integers or strings).

    Returns
    -------
    ids : array_like
        Indices of loaded wall tiles.
    area : array_like
        Area of each tile.
    edepo : array_like
        Energy/power deposition per tile.
    pdepo : array_like
        Particle/flux deposition per tile.
    iangle : array_like
        Angle of incidence i.e. the angle between particle velocity and
        the surface normal.

        When weights are included, the angle of incidence is the average of
        all markers that hit the tile weighted by the marker energy and
        weight. Otherwise it is the mean value of all markers without any
        weights.
    """
    self._require("_endstate")
    ids, energy, weight, pr, pphi, pz, pnorm, phi = self.getstate(
        "walltile", "ekin", "weight", "pr", "pphi", "pz", "pnorm", "phi",
        state="end", mode="prt", endcond="wall", ids=p_ids)
    if flags is not None:
        w = self.wall.read()
        flag = w["flag"]
        labels = w["labels"]
        if not isinstance(flags, list):
            flags = [flags]
        mask = np.array([
            labels[f] if f in labels else f for f in flags
        ])
        idx = np.where(np.isin(flag[ids-1], mask))[0]
        ids, energy, weight, pr, pphi, pz, pnorm, phi = (
            ids[idx], energy[idx], weight[idx], pr[idx], pphi[idx],
            pz[idx], pnorm[idx], phi[idx]
        )

    energy.convert_to_units("J")
    eunit = (energy.units * weight.units)
    try:
        eunit = (1 * eunit).to("W")
    except unyt.exceptions.UnitConversionError:
        pass
    wetted = np.unique(ids)
    edepo  = np.zeros(wetted.shape)
    pdepo  = np.zeros(wetted.shape)
    iangle = np.zeros(wetted.shape)
    if weights:
        edepo = edepo * eunit
        pdepo = pdepo * weight.units
    else:
        edepo = edepo * energy.units

    area, nvec = self.wall.area(normal=True)
    area = area[wetted-1]
    punit = np.array([pr, pphi, pz]) / pnorm

    # Sort markers by wall tile ID to make processing faster
    idx    = np.argsort(ids)
    ids    = ids[idx]
    energy = energy[idx]
    weight = weight[idx]
    punit  = punit[:,idx].v
    idx = np.append(np.argwhere(ids[1:] - ids[:-1]).ravel(), ids.size-1)

    # Convert unit vectors to cartesian basis
    punit = np.array([
        punit[0,:] * np.cos(phi) - punit[1,:] * np.sin(phi),
        punit[0,:] * np.sin(phi) + punit[1,:] * np.cos(phi),
        punit[2,:]])

    i0 = 0
    for i in range(wetted.size):
        i1 = idx[i] + 1
        vec = nvec[:,wetted[i]-1]
        dotprod = punit[0,i0:i1] * vec[0] + punit[1,i0:i1] * vec[1] \
                            + punit[2,i0:i1] * vec[2]
        angles = np.arccos(dotprod) * (180 / np.pi) * unyt.deg
        angles[angles.v > 90] = 180 * unyt.deg - angles[angles.v > 90]

        if weights: #include weights in the deposition
            pdepo[i]  = np.sum(weight[i0:i1])
            edepo[i]  = np.sum(energy[i0:i1]*weight[i0:i1])
            iangle[i] = np.sum(angles * energy[i0:i1] * weight[i0:i1]) \
                / edepo[i]
        else:
            pdepo[i]  = i1 - i0
            edepo[i]  = np.sum(energy[i0:i1])
            iangle[i] = np.mean(angles)
        i0 = i1

    return wetted, area, edepo, pdepo, iangle

def getwall_3dmesh(w_indices=None, p_ids=None):
    """Return 3D mesh representation of 3D wall and associated loads

    Parameters
    ----------
    w_indices : array_like, optional
        List of triangle indices for which the 3D mesh is made.
    p_ids : array_like, optional
        List of particle ids for which the wall loads are calculated.

    Returns
    -------
    wallmesh : Polydata
        Mesh representing the wall.

        The mesh cell data has fields:

        - "pload" particle load in units of prt/m^2 or prt/m^2s,
        - "eload" power/energy load in units of W/m^2 or J/m^2
        - "mload" marker load in units of markers
        - "iangle" angle of incidence (the angle between power flux and
            the surface normal) in deg
        - "label" flag of the wall element
    """
    wallmesh = pv.PolyData(
        *self.wall.noderepresentation(w_indices=w_indices))
    ids, area, eload, pload, iangle = self.getwall_loads(p_ids=p_ids)
    # get the marker load at each tile
    _, _, _, mload, _ = self.getwall_loads(weights=False, p_ids=p_ids)
    ids = ids - 1 # Convert IDs to indices

    w = self.wall.read()
    n_tri, flag = w["nelements"], w["flag"]
    wall_f = np.arange(n_tri, dtype=int)
    if w_indices is not None:
        wall_f = wall_f[w_indices]

    eload_tot = np.zeros((n_tri, )) + np.nan
    eload_tot[ids] = eload/area
    pload_tot = np.zeros((n_tri, )) + np.nan
    pload_tot[ids] = pload/area
    mload_tot = np.zeros((n_tri, )) + np.nan
    mload_tot[ids] = mload
    iangle_tot = np.zeros((n_tri, )) + np.nan
    iangle_tot[ids] = iangle

    wallmesh.cell_data["label"] = flag[wall_f]
    wallmesh.cell_data["pload"] = pload_tot[wall_f]
    wallmesh.cell_data["eload"] = eload_tot[wall_f]
    wallmesh.cell_data["mload"] = mload_tot[wall_f]
    wallmesh.cell_data["iangle"] = iangle_tot[wall_f]
    return wallmesh