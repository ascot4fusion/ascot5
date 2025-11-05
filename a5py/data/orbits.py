"""Orbit diagnostics.
"""
import numpy as np
import unyt
import a5py.physlib as physlib

import ctypes
from typing import Optional

import unyt
import numpy as np
from numpy.ctypeslib import ndpointer

from a5py import utils
from a5py.libascot import LIBASCOT, DataStruct
from a5py.exceptions import AscotMeltdownError
from a5py.data.access import InputVariant, Leaf, TreeMixin



# pylint: disable=too-few-public-methods
class Struct(DataStruct):
    _fields_ = [
        ("r", ctypes.POINTER(ctypes.c_double)),
        ("z", ctypes.POINTER(ctypes.c_double)),
        ("phi", ctypes.POINTER(ctypes.c_double)),
        ("p1", ctypes.POINTER(ctypes.c_double)),
        ("p2", ctypes.POINTER(ctypes.c_double)),
        ("p3", ctypes.POINTER(ctypes.c_double)),
        ("mileage", ctypes.POINTER(ctypes.c_double)),
        ("stamp", ctypes.POINTER(ctypes.c_double)),
        ("id", ctypes.POINTER(ctypes.c_size_t)),
        ("idx", ctypes.POINTER(ctypes.c_size_t)),
        ("charge", ctypes.POINTER(ctypes.c_int)),
        ("poincare", ctypes.POINTER(ctypes.c_int)),
        ("simmode", ctypes.POINTER(ctypes.c_int)),
    ]




class Orbit():
    """Orbit diagnostics that collect marker phase-space coordinates and related
    quantities at certain points along the marker orbit.
    """

    GYROORBIT     = 1
    GUIDINGCENTER = 2
    FIELDLINE     = 3


    def init(nmrk, npoint):
        orbit = Struct()

        orbit.r = (ctypes.c_double * nmrk * npoint)()
        orbit.z = (ctypes.c_double * nmrk * npoint)()
        orbit.phi = (ctypes.c_double * nmrk * npoint)()
        orbit.p1 = (ctypes.c_double * nmrk * npoint)()
        orbit.p2 = (ctypes.c_double * nmrk * npoint)()
        orbit.p3 = (ctypes.c_double * nmrk * npoint)()
        orbit.mileage = (ctypes.c_double * nmrk * npoint)()
        orbit.stamp = (ctypes.c_double * nmrk)()
        orbit.idx = (ctypes.c_size_t * nmrk)()
        orbit.id = (ctypes.c_size_t * nmrk * npoint)()
        orbit.charge = (ctypes.c_int * nmrk * npoint)()
        orbit.poincare = (ctypes.c_int * nmrk * npoint)()
        orbit.simmode = (ctypes.c_int * nmrk * npoint)()

    def get(self, inistate, endstate, *qnt):
        """Return marker quantity.

        This function accesses the orbit data within the HDF5 file and uses that
        to evaluate the queried quantity. The evaluated quantity at a given
        position corresponds to that mode which was active at the simulation:
        GO simulations return particle phase-space, GC guiding-center
        phase-space and hybrid depends on whether marker was GO or GC at that
        moment when data point was written.

        Parameters
        ----------
        inistate : :class:`State`
            Inistate is needed to evaluate some orbit quantities.
        endstate : :class:`State`
            Endstate is needed to evaluate some orbit quantities.
        *qnt : str
            Names of the quantities.

        Returns
        -------
        *value : array_like
            The quantities as an array ordered by marker ID (major) and mileage
            (minor).
        """
        # Prepare helper variables and functions
        def _val(q, mask=None):
            """Read quantity from HDF5.
            """
            with self as h5:
                if q in h5:
                    qnt = fileapi.read_data(h5, q)
                    if q == "weight":
                        qnt *= unyt.unyt_quantity.from_string("particles/s")
                    return qnt if mask is None else qnt[mask]
            return None

        # Sort using the fact that inistate.get return values ordered by ID
        # and also np.unique returns indices that produce a sorted array.
        mode    = _val("simmode")
        _, idx  = np.unique(_val("ids").v, return_inverse=True)
        mass    = inistate.get("mass")[0][idx]
        time    = inistate.get("time")[0][idx]
        connlen = inistate.get("mileage")[0][idx] - _val("mileage")

        # Only field lines are constant in time
        if not Orbits.FIELDLINE in mode: time = time + _val("mileage")

        def _eval(q, mask=None):
            """Evaluate input quantities at marker position.
            """
            return self._root._ascot.input_eval(
                _val("r", mask=mask), _val("phi",  mask=mask),
                _val("z", mask=mask), time[mask], *[q])

        return Orbits._getactual(mass, time, connlen, mode, _val, _eval, *qnt)

    @staticmethod
    def _getactual(mass, time, totmil, mode, _val, _eval, *qnt):
        """Calculate orbit quantities using the helper functions and data.

        Parameters
        ----------
        mass : array_like, (n,)
            Marker mass repeated at each orbit point.
        time : array_like, (n,)
            Laboratory time at each orbit point.
        totmil : array_like, (n,)
            Total mileage (from endstate) repeated at each orbit point.
        mode : array_like, (n,)
            The simmode value at each orbit point.
        _val : callable
            Function that returns a stored orbit parameter at masked positions.

            ``_val(qnt : str, mask : array_like) -> value``
        _eval : callable
            Function that returns interpolated input quantity at masked
            positions.

            ``_eval(qnt : str, mask : array_like) -> value``
        *qnt : str
            Names of the quantities.

        Returns
        -------
        *value : array_like
            The quantities as an array ordered by marker ID (major) and mileage
            (minor).
        """
        items = [None]*len(qnt)
        def add(q, val):
            if q in qnt:
                i = qnt.index(q)
                if items[i] is None:
                    items[i] = val()
                else:
                    # For hybrid mode GC values are appended by GO
                    items[i].convert_to_base("ascot")
                    v = val().in_base("ascot")
                    items[i] = np.append(v.v, items[i].v) * v.units

        # Some helper quantities as functions so they are evaluated only
        # when needed. Mask is used to separate GOs, GCs, and FLs since
        # they are all stored in the same data.
        bvec = lambda mask : unyt.unyt_array(
            [_val("br", mask), _val("bphi", mask), _val("bz", mask)])
        pvecprt = lambda mask : unyt.unyt_array(
            [_val("pr", mask), _val("pphi", mask), _val("pz", mask)])
        def pvecgc(mask):
            bnorm = np.sqrt(np.sum(bvec(mask)**2,axis=0))
            pnorm = physlib.momentum_muppar(
                mass[mask], _val("mu", mask),
                _val("ppar", mask), bnorm)
            bhat  = bvec(mask) / bnorm
            e1 = np.zeros(bhat.shape)
            e1[2,:] = 1
            e2 = np.cross(bhat.T, e1.T).T
            e1 = e2 / np.sqrt(np.sum(e2**2, axis=0))
            e2 = np.cross(bhat.T, e1.T).T
            perphat = -np.sin(_val("zeta", mask))*e1 \
                - np.cos(_val("zeta", mask))*e2
            return bhat * _val("ppar", mask) \
                + perphat * np.sqrt(pnorm**2 - _val("ppar", mask)**2)
        vvecprt = lambda mask : physlib.velocity_momentum(mass[mask],
                                                          pvecprt(mask))
        vvecgc = lambda mask : physlib.velocity_momentum(mass[mask],
                                                         pvecgc(mask))
        # Common quantities
        add("ids", lambda : _val("ids"))
        add("r", lambda : _val("r"))
        add("z", lambda : _val("z"))
        add("phi", lambda : _val("phi"))
        add("phimod", lambda : np.mod(_val("phi"), 2 * np.pi * unyt.rad))
        add("theta", lambda : _val("theta"))
        add("thetamod", lambda : np.mod(_val("theta"), 2 * np.pi * unyt.rad))
        add("x", lambda : physlib.pol2cart(_val("r"), _val("phi"))[0])
        add("y", lambda : physlib.pol2cart(_val("r"), _val("phi"))[1])
        add("br", lambda : _val("br"))
        add("bz", lambda : _val("bz"))
        add("bphi", lambda : _val("bphi"))
        add("bnorm", lambda : np.sqrt( _val("br")**2 + _val("bphi")**2
                                       + _val("bz")**2 ))
        add("rho", lambda : _val("rho"))
        add("psi", lambda : _eval("psi"))
        add("mileage", lambda : _val("mileage"))
        add("mass", lambda : mass)
        add("charge", lambda : _val("charge"))
        add("weight", lambda : _val("weight"))
        add("time", lambda : time)
        add("pncrid", lambda : _val("pncrid"))
        add("pncrdir", lambda : _val("pncdir"))
        add("connlen", lambda : totmil)

        mask = mode > 0
        if Orbits.GYROORBIT in mode:
            mask = mode == Orbits.GYROORBIT
            add("pr", lambda : _val("pr", mask))
            add("pz", lambda : _val("pz", mask))
            add("pphi", lambda : _val("pphi", mask))
            add("ppar", lambda : physlib.ppar_momentum(pvecprt(mask),
                                                       bvec(mask)))
            if "pperp" in qnt:
                pnorm2 = np.sum( pvecprt(mask)**2, axis=0 )
                pperp2 = physlib.ppar_momentum(pvecprt(mask), bvec(mask))**2
                add("pperp", lambda : np.sqrt(pnorm2 - pperp2))
            add("pnorm", lambda : np.sqrt( np.sum( pvecprt(mask)**2, axis=0 ) ))
            add("vr", lambda : vvecprt(mask)[0,:])
            add("vz", lambda : vvecprt(mask)[2,:])
            add("vphi", lambda : vvecprt(mask)[1,:])
            add("vpar", lambda : physlib.vpar_momentum(
                mass[mask], pvecprt(mask), bvec(mask)))
            if "vperp" in qnt:
                vpar = physlib.vpar_momentum(
                    mass[mask], pvecprt(mask), bvec(mask))
                pnorm = np.sqrt(np.sum(pvecprt(mask)**2, axis=0))
                vnorm = physlib.velocity_momentum(mass[mask], pnorm)
                add("vperp", lambda : np.sqrt(vnorm**2 - vpar**2))
            if "vnorm" in qnt:
                pnorm = np.sqrt(np.sum(pvecprt(mask)**2, axis=0))
                add("vnorm", lambda : physlib.velocity_momentum(mass[mask],
                                                                pnorm))
            add("ekin", lambda : physlib.energy_momentum(mass[mask],
                                                         pvecprt(mask)))
            add("pitch", lambda : physlib.pitch_momentum(pvecprt(mask),
                                                         bvec(mask)))
            add("mu", lambda : physlib.mu_momentum(mass[mask], pvecprt(mask),
                                                   bvec(mask)))
            add("zeta", lambda : mass[mask]*np.nan) # TODO implement properly
            add("ptor", lambda : physlib.torcanangmom_momentum(
                _val("charge", mask), _val("r", mask), pvecprt(mask),
                _eval("psi", mask)))

        if Orbits.GUIDINGCENTER in mode:
            mask = mode == Orbits.GUIDINGCENTER
            add("pr", lambda : pvecgc(mask)[0,:])
            add("pz", lambda : pvecgc(mask)[2,:])
            add("pphi", lambda : pvecgc(mask)[1,:])
            add("ppar", lambda : _val("ppar", mask))
            if "pperp" in qnt:
                pnorm = physlib.momentum_muppar(
                    mass[mask], _val("mu", mask),
                    _val("ppar", mask), bvec(mask))
                add("pperp", lambda : np.sqrt(pnorm**2 - _val("ppar")**2))
            add("pnorm", lambda : physlib.momentum_muppar(
                mass[mask], _val("mu", mask),
                _val("ppar", mask), bvec(mask)))
            add("vr", lambda : vvecgc(mask)[0,:])
            add("vz", lambda : vvecgc(mask)[2,:])
            add("vphi", lambda : vvecgc(mask)[1,:])
            add("vpar", lambda : physlib.vpar_muppar(
                mass[mask], _val("mu", mask), _val("ppar", mask), bvec(mask)))
            if "vperp" in qnt:
                pnorm = physlib.momentum_muppar(
                    mass[mask], _val("mu", mask),
                    _val("ppar", mask), bvec(mask))
                pperp = np.sqrt(pnorm**2 - _val("ppar", mask)**2)
                gamma = physlib.gamma_momentum(mass[mask], pnorm)
                add("vperp", lambda : pperp / (gamma * mass[mask]))
            if "vnorm" in qnt:
                pnorm = physlib.momentum_muppar(
                    mass[mask], _val("mu", mask),
                    _val("ppar", mask), bvec(mask))
                gamma = physlib.gamma_momentum(mass[mask], pnorm)
                add("vnorm", lambda : pnorm / (gamma * mass[mask]))
            add("ekin", lambda : physlib.energy_muppar(
                mass[mask], _val("mu", mask), _val("ppar", mask), bvec(mask)))
            add("pitch", lambda : physlib.pitch_muppar(
                mass[mask], _val("mu", mask), _val("ppar", mask), bvec(mask)))
            add("mu", lambda : _val("mu", mask))
            add("zeta", lambda : _val("zeta", mask))
            add("ptor", lambda : physlib.torcanangmom_ppar(
                _val("charge", mask), _val("r", mask), _val("ppar", mask),
                bvec(mask), _eval("psi", mask)))

        for i in range(len(items)):
            if items[i] is None:
                raise ValueError("Unknown quantity in " + qnt[i])

        # Sort first by IDs and then by mileage
        ids  = _val("ids", mask=mask).v
        mile = _val("mileage", mask=mask).v
        idx  = np.lexsort((mile, ids))
        for i in range(len(items)):
            items[i] = items[i][idx]
            items[i].convert_to_base("ascot")
        return items

    @staticmethod
    def listqnts():
        """List all available quantities.
        """
        out = {
            "r":        "R coordinate",
            "z":        "z coordinate",
            "phi":      "Toroidal coordinate (cumulative)",
            "phimod":   "Toroidal coordinate",
            "theta":    "Poloidal coordinate (cumulative)",
            "thetamod": "Poloidal coordinate",
            "x":        "x coordinate",
            "y":        "y coordinate",
            "pr":       "Momentum R component",
            "pz":       "Momentum z component",
            "pphi":     "Momentum phi component",
            "ppar":     "Momentum component parallel to B",
            "pperp":    "Momentum component perpendicular to B",
            "pnorm":    "Momentum norm",
            "vr":       "Velocity R component",
            "vz":       "Velocity z component",
            "vphi":     "Velocity phi component",
            "vpar":     "Velocity component parallel to B",
            "vperp":    "Velocity component perpendicular to B",
            "vnorm":    "Velocity norm",
            "br":       "Magnetic field R component at the marker position",
            "bz":       "Magnetic field z component at the marker position",
            "bphi":     "Magnetic field phi component at the marker position",
            "bnorm":    "Magnetic field strength at the marker position",
            "ekin":     "Kinetic energy",
            "pitch":    "vpa / vnorm",
            "mu":       "Magnetic moment",
            "zeta":     "Gyroangle",
            "psi":      "Poloidal flux at the marker position",
            "rho":      "Square root of normalized psi",
            "ptor":     "Canonical toroidal angular momentum",
            "mass":     "Mass",
            "charge":   "Charge",
            "time":     "Current laboratory time",
            "mileage":  "Laboratory time elapsed in simulation",
            "weight":   "How many physical particles a marker represents",
            "ids":      "Marker ID",
            "connlen":  "Connection length for lost markers",
            "pncrid":   "Poincaré plane this point corresponds to",
            "pncrdir":  "Direction at which Poincaré plane was crossed",
        }
        return out
