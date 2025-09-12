"""Marker initial and final phase-space positions and related quantities.
"""
import ctypes
import numpy as np
import unyt
import a5py.physlib as physlib

from ..access import Diagnostic


class MarkerState(Diagnostic):
    """Marker initial and final phase-space positions and related quantities.
    """

    class Structure(ctypes.Structure):
        """Python wrapper for the particle state struct in particle.h."""

        _pack_ = 1
        _fields_ = [
            ('r', ctypes.c_double),
            ('phi', ctypes.c_double),
            ('z', ctypes.c_double),
            ('ppar', ctypes.c_double),
            ('mu', ctypes.c_double),
            ('zeta', ctypes.c_double),
            ('rprt', ctypes.c_double),
            ('phiprt', ctypes.c_double),
            ('zprt', ctypes.c_double),
            ('p_r', ctypes.c_double),
            ('p_phi', ctypes.c_double),
            ('p_z', ctypes.c_double),
            ('mass', ctypes.c_double),
            ('charge', ctypes.c_double),
            ('anum', ctypes.c_int32),
            ('znum', ctypes.c_int32),
            ('weight', ctypes.c_double),
            ('time', ctypes.c_double),
            ('mileage', ctypes.c_double),
            ('cputime', ctypes.c_double),
            ('rho', ctypes.c_double),
            ('theta', ctypes.c_double),
            ('id', ctypes.c_int64),
            ('endcond', ctypes.c_int64),
            ('walltile', ctypes.c_int64),
            ('B_r', ctypes.c_double),
            ('B_phi', ctypes.c_double),
            ('B_z', ctypes.c_double),
            ('B_r_dr', ctypes.c_double),
            ('B_phi_dr', ctypes.c_double),
            ('B_z_dr', ctypes.c_double),
            ('B_r_dphi', ctypes.c_double),
            ('B_phi_dphi', ctypes.c_double),
            ('B_z_dphi', ctypes.c_double),
            ('B_r_dz', ctypes.c_double),
            ('B_phi_dz', ctypes.c_double),
            ('B_z_dz', ctypes.c_double),
            ('err', ctypes.c_uint64),
            ]


    def __init__(self, nmrk: int, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._struct_ = ctypes.POINTER(particle_state)()
        self.nmrk = nmrk

    def combine(self):
        """Combines two into one (copy)."""

    def get(self, *qnt, mode="gc"):
        """Return marker quantity.

        This function accesses the state data within the HDF5 file and uses that
        to evaluate the queried quantity.

        Parameters
        ----------
        qnt : str
            Name of the quantities.
        mode : {"prt", "gc"}, optional
            Is the quantity evaluated in particle or guiding center phase-space.

        Returns
        -------
        value : array_like
            The quantities as an array ordered by marker ID.
        """
        def _val(q):
            """Read quantity from HDF5.
            """
            with self as h5:
                if q in h5:
                    return fileapi.read_data(h5, q)
            return None

        def _eval(r, phi, z, t, *q):
            """Evaluate quantity with ascotpy.
            """
            return self._root._ascot.input_eval(r, phi, z, t, *q)

        return State._getactual(mode, _val, _eval, *qnt)

    @staticmethod
    def _getactual(mode, _val, _eval, *qnt):
        """Calculate marker quantities using the helper functions and data.

        Parameters
        ----------
        mode : {"prt", "gc"}
            Phase-space where quantity is evaluated.
        _val : callable
            Function that returns a stored marker parameter.

            ``_val(qnt : str, mask : array_like) -> value``
        _eval : callable
            Function that returns interpolated input quantity at masked
            positions at given position.

            ``_eval(r, phi, z, qnt : str) -> value``
        *qnt : str
            Names of the quantities.

        Returns
        -------
        *value : array_like
            The quantities as an array ordered by marker ID.
        """
        items = [None]*len(qnt)
        def add(q, val):
            if q in qnt:
                items[qnt.index(q)] = val()

        evalprt = lambda *q : _eval(_val("rprt"), _val("phiprt"), _val("zprt"),
                                    _val("time"), *q)
        evalgc  = lambda *q : _eval(_val("r"), _val("phi"), _val("z"),
                                    _val("time"), *q)
        bvecprt = lambda : unyt.unyt_array(evalprt("br", "bphi", "bz"))
        bvecgc  = lambda : unyt.unyt_array([_val("br"), _val("bphi"),
                                            _val("bz")])
        pvecprt = lambda : unyt.unyt_array([_val("prprt"), _val("pphiprt"),
                                _val("pzprt")])
        def pvecgc():
            """Evaluate guiding center momentum vector.
            """
            bnorm = np.sqrt(np.sum(bvecgc()**2,axis=0))
            pnorm = physlib.momentum_muppar(_val("mass"), _val("mu"),
                                            _val("ppar"), bnorm)
            bhat  = bvecgc() / bnorm
            e1 = np.zeros(bhat.shape)
            e1[2,:] = 1
            e2 = np.cross(bhat.T, e1.T).T
            e1 = e2 / np.sqrt(np.sum(e2**2, axis=0))
            e2 = np.cross(bhat.T, e1.T).T
            perphat = -np.sin(_val("zeta"))*e1-np.cos(_val("zeta"))*e2
            return bhat * _val("ppar") \
                + perphat * np.sqrt(pnorm**2 - _val("ppar")**2)
        vvecprt = lambda : physlib.velocity_momentum(_val("mass"), pvecprt())
        vvecgc = lambda : physlib.velocity_momentum(_val("mass"), pvecgc())

        # Dimensionless quantities common for particle and guiding-center
        add("ids", lambda : _val("ids").v)
        add("anum", lambda : _val("anum").v)
        add("znum", lambda : _val("znum").v)
        add("walltile", lambda : _val("walltile").v)
        add("errormsg", lambda : _val("errormsg").v)
        add("errorline", lambda : _val("errorline").v)
        add("errormod", lambda : _val("errormod").v)
        if "endcond" in qnt:
            item = _val("endcond").v
            err = _val("errormsg").v
            item = item << 2
            item[err > 0] = item[err > 0] & State._ABORTED
            item[item==0] = State._NONE
            add("endcond", lambda : item)

        # Dimensional quantities common for particle and guiding-center
        add("mass", lambda : _val("mass"))
        add("charge", lambda : _val("charge"))
        add("time", lambda : _val("time"))
        add("cputime", lambda : _val("cputime"))
        add("mileage", lambda : _val("mileage"))
        add("weight", lambda : _val("weight"))

        # Quantities that have to be evaluated separately (keep this in same
        # order for both and in the list method if new quantities are added)
        if mode == "prt":
            add("r", lambda : _val("rprt"))
            add("z", lambda : _val("zprt"))
            add("phi", lambda : _val("phiprt"))
            add("phimod", lambda : np.mod(_val("phiprt"), 2 * np.pi * unyt.rad))
            if "theta" in qnt:
                axis = evalprt("axisr", "axisz")
                thetaprt = physlib.cart2pol(
                    _val("rprt")-axis[0], _val("zprt")-axis[1])[1]
                thetacum = np.floor(_val("phi") / (2*np.pi * unyt.rad))
                theta = 2*np.pi * thetacum * unyt.rad + thetaprt
                add("theta", lambda : np.mod(theta + np.pi*unyt.rad,
                                             2*np.pi*unyt.rad))
            if "thetamod" in qnt:
                axis = evalprt("axisr", "axisz")
                theta = physlib.cart2pol(
                    _val("rprt")-axis[0], _val("zprt")-axis[1])[1]
                add("thetamod", lambda : np.mod(theta + np.pi*unyt.rad,
                                                2*np.pi*unyt.rad))
            add("x", lambda : physlib.pol2cart(_val("rprt"), _val("phiprt"))[0])
            add("y", lambda : physlib.pol2cart(_val("rprt"), _val("phiprt"))[1])
            add("pr", lambda : _val("prprt"))
            add("pz", lambda : _val("pzprt"))
            add("pphi", lambda : _val("pphiprt"))
            add("ppar", lambda : physlib.ppar_momentum(pvecprt(), bvecprt()))
            if "pperp" in qnt:
                pnorm2 = np.sum( pvecprt()**2, axis=0 )
                ppar2  = physlib.ppar_momentum(pvecprt(), bvecprt())**2
                add("pperp", lambda : np.sqrt(pnorm2 - ppar2))
            if "pnorm" in qnt:
                add("pnorm", lambda : np.sqrt( np.sum( pvecprt()**2, axis=0 ) ))
            add("vr", lambda : vvecprt()[0,:])
            add("vz", lambda : vvecprt()[2,:])
            add("vphi", lambda : vvecprt()[1,:])
            add("vpar", lambda : physlib.vpar_momentum(_val("mass"), pvecprt(),
                                              bvecprt()))
            if "vperp" in qnt:
                vpar = physlib.vpar_momentum(
                    _val("mass"), pvecprt(), bvecprt())
                pnorm = np.sqrt(np.sum(pvecprt()**2, axis=0))
                vnorm = physlib.velocity_momentum(_val("mass"), pnorm)
                add("vperp", lambda : np.sqrt(vnorm**2 - vpar**2))
            if "vnorm" in qnt:
                pnorm = np.sqrt(np.sum(pvecprt()**2, axis=0))
                add("vnorm", lambda : physlib.velocity_momentum(_val("mass"),
                                                                pnorm))
            add("br", lambda : evalprt("br"))
            add("bz", lambda : evalprt("bz"))
            add("bphi", lambda : evalprt("bphi"))
            add("bnorm", lambda : evalprt("bnorm"))
            add("ekin", lambda : physlib.energy_momentum(_val("mass"),
                                                         pvecprt()))
            add("pitch", lambda : physlib.pitch_momentum(pvecprt(), bvecprt()))
            add("mu", lambda : physlib.mu_momentum(_val("mass"), pvecprt(),
                                          bvecprt()))
            add("zeta", lambda : _val("zeta"))
            add("psi", lambda : evalprt("psi"))
            add("rho", lambda : evalprt("rho"))
            add("ptor", lambda : physlib.torcanangmom_momentum(
                _val("charge"), _val("rprt"), pvecprt(), evalprt("psi")))

        if mode == "gc":
            add("r", lambda : _val("r"))
            add("z", lambda : _val("z"))
            add("phi", lambda : _val("phi"))
            add("phimod", lambda : np.mod(_val("phi"), 2*np.pi * unyt.rad))
            add("theta", lambda : _val("theta"))
            add("thetamod", lambda : np.mod(_val("theta"), 2*np.pi * unyt.rad))
            add("x", lambda : physlib.pol2cart(_val("r"), _val("phi"))[0])
            add("y", lambda : physlib.pol2cart(_val("r"), _val("phi"))[1])
            add("pr", lambda : pvecgc()[0,:])
            add("pz", lambda : pvecgc()[2,:])
            add("pphi", lambda : pvecgc()[1,:])
            add("ppar", lambda : _val("ppar"))
            if "pperp" in qnt:
                pnorm = physlib.momentum_muppar(_val("mass"), _val("mu"),
                                                _val("ppar"), bvecgc())
                add("pperp", lambda : np.sqrt(pnorm**2 - _val("ppar")**2))
            add("pnorm",lambda :  physlib.momentum_muppar(
                _val("mass"), _val("mu"), _val("ppar"), bvecgc()))
            add("vr", lambda : vvecgc()[0,:])
            add("vz", lambda : vvecgc()[2,:])
            add("vphi", lambda : vvecgc()[1,:])
            add("vpar", lambda : physlib.vpar_muppar(_val("mass"), _val("mu"),
                                                     _val("ppar"), bvecgc()))
            if "vperp" in qnt:
                pnorm = physlib.momentum_muppar(_val("mass"), _val("mu"),
                                                _val("ppar"), bvecgc())
                pperp = np.sqrt(pnorm**2 - _val("ppar")**2)
                gamma = physlib.gamma_momentum(_val("mass"), pnorm)
                add("vperp", lambda : pperp / (gamma * _val("mass")))
            if "vnorm" in qnt:
                pnorm = physlib.momentum_muppar(_val("mass"), _val("mu"),
                                                _val("ppar"), bvecgc())
                gamma = physlib.gamma_momentum(_val("mass"), pnorm)
                add("vnorm", lambda : pnorm / (gamma * _val("mass")))
            add("br", lambda : _val("br"))
            add("bz", lambda : _val("bz"))
            add("bphi", lambda : _val("bphi"))
            add("bnorm", lambda : np.sqrt(_val("br")**2 + _val("bphi")**2 +
                                 _val("bz")**2))
            add("ekin", lambda : physlib.energy_muppar(_val("mass"), _val("mu"),
                                              _val("ppar"), bvecgc()))
            add("pitch", lambda : physlib.pitch_muppar(_val("mass"), _val("mu"),
                                              _val("ppar"), bvecgc()))
            add("mu", lambda : _val("mu"))
            add("zeta", lambda : _val("zeta"))
            add("psi", lambda : evalgc("psi"))
            add("rho", lambda : _val("rho"))
            add("ptor", lambda : physlib.torcanangmom_ppar(
                _val("charge"), _val("r"), _val("ppar"),
                bvecgc(), evalgc("psi")))

        for i in range(len(items)):
            if items[i] is None:
                raise ValueError("Unknown quantity in " + qnt[i])

        # Order by ID
        idx = _val("ids").argsort()
        for i in range(len(items)):
            items[i] = items[i][idx]
            try:
                items[i].convert_to_base("ascot")
            except AttributeError:
                # Non-dimensional
                pass
        return items

    @staticmethod
    def listqnts():
        """List available quantities.
        """

        from typing import Literal

        Literal["ids"]
        out = {
            "ids":      "Marker ID",
            "anum":     "Mass number",
            "znum":     "Charge number",
            "walltile": "Element ID if the marker hit wall",
            "endcond":  "Bitarray showing active end conditions",
            "errormsg": "Error message if marker was aborted",
            "errorline":"Line where (if) marker was aborted",
            "errormod": "Name of the module where (if) marker was aborted",
            "mass":     "Mass",
            "charge":   "Charge",
            "time":     "Current laboratory time",
            "cputime":  "CPU time elapsed in simulation",
            "mileage":  "Laboratory time elapsed in simulation",
            "weight":   "How many physical particles a marker represents",
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
        }
        return out
