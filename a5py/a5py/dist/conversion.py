"""
Routines for converting distribution abscissae to other coordinates.

File: conversion.py
"""

import numpy as np
import copy
import scipy.constants as constants
import itertools

from scipy.interpolate import griddata, RectBivariateSpline


def convert_ppappe_to_Exi(dist, masskg, E_edges=None, xi_edges=None):
    """
    Converts vpa and vpe distribution abscissae to energy and pitch.

    This function operates by looping through all other coordinates except
    ppa and ppe, and at each loop calculates
    f_Exi(E, xi) = f_ppappe(ppa(E_i, xi_i), ppe(E_i, xi_i)) where E_i and
    xi_i are grid points of the new energy-pitch distribution. Interpolation
    is done bilinearly.

    Energy is in electronvolts and pitch is ppa/(ppa^2 + ppe^2)^0.5. The
    transformation is not relativistic.

    Args:
        dist : dict_like <br>
            A ppa-ppe distribution. May hold other dimensions as well.
        masskg : float <br>
            Mass of the species (required for energy conversion) in kg. Note
            that distribution is assumed to consist of markers with equal mass.
        E_edges : array_like or int, optional <br>
            Energy grid edges in the new distribution. If not given,
            linspace(0, Emax, 10) will be used where Emax is
            e*0.5*masskg*max(vpa^2, vpe^2). If an integer is given then
            linspace(0, Emax, E_edges) is used.
        xi_edges : array_like or int, optional <br>
            Pitch grid edges in the new distribution. If not given,
            linspace(-1, 1, 10) will be used. If an integer is given then
            linspace(-1, 1, xi_edges) is used.

    Returns:
       Energy-pitch distribution dictionary whose other dimensions are same as
       in input.
    """

    if E_edges is None:
        pmax2 = np.maximum( dist["ppar_edges"][-1]*dist["ppar_edges"][-1],
                           dist["pperp_edges"][-1]*dist["pperp_edges"][-1] )
        gamma = np.sqrt( 1 +  pmax2 / ( masskg * constants.c) **2)
        Emax = (1/constants.e) * (gamma - 1) * masskg * constants.c**2
        E_edges = np.linspace(0, Emax, 10)
    if xi_edges is None:
        xi_edges = np.linspace(-1, 1, 10)

    if not isinstance(E_edges, np.ndarray):
        pmax2 = np.maximum( dist["ppar_edges"][-1]*dist["ppar_edges"][-1],
                           dist["pperp_edges"][-1]*dist["pperp_edges"][-1] )
        gamma = np.sqrt( 1 +  pmax2 / ( masskg * constants.c) **2)
        Emax = (1/constants.e) * (gamma - 1) * masskg * constants.c**2
        E_edges = np.linspace(0, Emax, E_edges)

    if not isinstance(xi_edges, np.ndarray):
        xi_edges = np.linspace(-1, 1, xi_edges)

    ## Create E-xi distribution ##
    Exidist = copy.deepcopy(dist)

    # Remove vpa and vpe components
    del Exidist["distribution"]
    Exidist["abscissae"].remove("ppar")
    Exidist["abscissae"].remove("pperp")
    for k in list(Exidist):
        if "ppar" in k or "pperp" in k:
            del Exidist[k]

    # Add E and xi abscissae and initialize a new density
    Exidist["abscissae"].insert(3, "energy")
    Exidist["abscissae"].insert(4, "pitch")
    abscissae = Exidist["abscissae"]

    Exidist["energy"]       = (E_edges[0:-1] + E_edges[1:]) / 2
    Exidist["energy_edges"] = E_edges
    Exidist["nenergy"]      = Exidist["energy"].size

    Exidist["pitch"]       = (xi_edges[0:-1] + xi_edges[1:]) / 2
    Exidist["pitch_edges"] = xi_edges
    Exidist["npitch"]      = Exidist["pitch"].size

    dims = []
    for a in abscissae:
        dims.append(Exidist["n" + a])

    Exidist["distribution"]  = np.zeros(tuple(dims))

    # Transform E-xi grid to points in (vpa,vpa) space that are used in
    # interpolation.
    xig, Eg = np.meshgrid(Exidist["pitch"], Exidist["energy"])
    pg   = np.sqrt( ( Eg * constants.e / constants.c + masskg.v * constants.c )**2
                     - (masskg.v * constants.c)**2 )
    ppag = ( xig * pg ).ravel()
    ppeg = (np.sqrt(1 - xig*xig) * pg).ravel()

    # Coordinate transform Jacobian: dvpa dvpe = |jac| dE dxi
    # Jacobian for transform (ppa, ppe) -> (p, xi) is p / sqrt(1-xi^2)
    # because jac = dppa / dp  = xi, dppe / dp  = sqrt(1-xi^2)
    #               dppa / dxi = p,  dppe / dxi = -xi p / sqrt(1-xi^2),
    # and the Jacobian for (p, xi) -> (E, xi) is e E0 / c^2 p when
    # E is in electronvolts. Therefore the combined Jacobian is
    # (e E0 / c^2) / sqrt(1-xi*xi).
    E0  = np.sqrt( (pg*constants.c)**2 + (masskg.v*constants.c**2)**2 )
    jac = (constants.e * E0 / constants.c**2) / np.sqrt(1 - xig*xig)

    # Interpolate.
    ranges = []
    for a in dist["abscissae"]:
        if a != "ppar" and a != "pperp":
            ranges.append(range(dist["n" + a]))

    iE   = Exidist["abscissae"].index("energy")
    ippa = dist["abscissae"].index("ppar")
    for itr in itertools.product(*ranges):

        idx = []
        for i in range(0, ippa):
            idx.append(itr[i])

        idx.append(slice(None))
        idx.append(slice(None))

        for i in range(ippa, len(dist["abscissae"])-2):
            idx.append(itr[i])

        idx = tuple(idx)

        if dist["ppar"].size == 1 and dist["pperp"].size == 1:
            d = np.ones( (Exidist["nenergy"], Exidist["npitch"]) ) \
                * dist["distribution"][idx]
        else:
            f = RectBivariateSpline(
                dist["ppar"], dist["pperp"],
                np.squeeze(dist["distribution"][idx]),
                kx=1, ky=1)
            d = np.reshape(f.ev(ppag, ppeg),
                           (Exidist["nenergy"], Exidist["npitch"]))

        Exidist["distribution"][idx] = d * jac

    return Exidist
