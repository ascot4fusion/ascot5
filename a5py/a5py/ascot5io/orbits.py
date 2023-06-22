"""Orbit data IO module.
"""
import numpy as np
import h5py
import unyt
import a5py.physlib as physlib

from .coreio import fileapi
from .coreio.treedata import DataContainer


class Orbits(DataContainer):
    """Orbit diagnostics that collect marker phase-space coordinates and related
    quantities at certain points along the marker orbit.
    """

    GYROORBIT = 1
    GUIDINGCENTER = 2
    FIELDLINE = 3

    def read_hdf5(self, fn, qid):
        """
        Read orbit output from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing orbit data.
        """

        path = "results/run_" + qid + "/orbit"

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        return out

    @staticmethod
    def write_hdf5(fn, data, desc=None):
        """
        Write state data in HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        desc : str, optional <br>
            Input description.
        """

        gname = "results/" + run + "/orbit"

        N = data["id"].size

        with h5py.File(fn, "a") as f:
            g = f.create_group(gname)

            g.create_dataset("id",   (N,1), data=data["id"],   dtype="i8")
            g.create_dataset("time", (N,1), data=data["time"], dtype="f8")

            if "pncrid" in data:
                g.create_dataset("pncrid", (N,1), data=data["pncrid"], dtype="i4")

            if "charge" in data:
                g.create_dataset("charge", (N,1), data=data["charge"], dtype="i4")

            if "weight" in data:
                g.create_dataset("weight", (N,1), data=data["weight"], dtype="f8")

            if "pr" in data:
                g.create_dataset("r",     (N,1), data=data["r"],     dtype="f8")
                g.create_dataset("phi",   (N,1), data=data["phi"],   dtype="f8")
                g.create_dataset("z",     (N,1), data=data["z"],     dtype="f8")

                g.create_dataset("pr",    (N,1), data=data["pr"],    dtype="f8")
                g.create_dataset("pphi",  (N,1), data=data["pphi"],  dtype="f8")
                g.create_dataset("pz",    (N,1), data=data["pz"],    dtype="f8")

                g.create_dataset("rho",   (N,1), data=data["rho"],   dtype="f8")
                g.create_dataset("theta", (N,1), data=data["theta"], dtype="f8")
                g.create_dataset("br",    (N,1), data=data["br"],    dtype="f8")
                g.create_dataset("bphi",  (N,1), data=data["bphi"],  dtype="f8")
                g.create_dataset("bz",    (N,1), data=data["bz"],    dtype="f8")
            else:
                g.create_dataset("r",     (N,1), data=data["r"],     dtype="f8")
                g.create_dataset("phi",   (N,1), data=data["phi"],   dtype="f8")
                g.create_dataset("z",     (N,1), data=data["z"],     dtype="f8")

                g.create_dataset("rho",   (N,1), data=data["rho"],   dtype="f8")
                g.create_dataset("theta", (N,1), data=data["theta"], dtype="f8")
                g.create_dataset("br",    (N,1), data=data["br"],    dtype="f8")
                g.create_dataset("bphi",  (N,1), data=data["bphi"],  dtype="f8")
                g.create_dataset("bz",    (N,1), data=data["bz"],    dtype="f8")

            if "ppar" in data:
                g.create_dataset("ppar",  (N,1), data=data["ppar"],  dtype="f8")
                g.create_dataset("mu",    (N,1), data=data["mu"],    dtype="f8")
                g.create_dataset("zeta",  (N,1), data=data["zeta"],  dtype="f8")

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
        inistate : State
            Inistate is needed to evaluate some orbit quantities.
        endstate : State
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
                    return fileapi.read_data(h5, q)
            return None

        mode = _val("simmode")
        ids  = inistate.get("ids")[0]
        val     = inistate.get("mass")[0].v
        f       = lambda x: val[ids == x]
        mass = np.array([f(x) for x in _val("ids")]).ravel() * unyt.amu
        val     = inistate.get("time")[0]
        f       = lambda x: val[ids == x]
        time = np.array([f(x) for x in _val("ids")]).ravel() * unyt.s
        val     = inistate.get("mileage")[0]
        f       = lambda x: val[ids == x]
        connlen = np.array([f(x) for x in _val("ids")]).ravel() * unyt.s
        connlen -= _val("mileage")
        if not Orbits.FIELDLINE in mode:
            time = time + _val("mileage")

        def _eval(q, mask=None):
            """Evaluate input quantities at marker position.
            """
            return self._root._ascot.input_eval(
                _val("r", mask=mask), _val("phi",  mask=mask),
                _val("z", mask=mask), time[mask], *[q])

        return Orbits.getactual(mass, time, connlen, mode, _val, _eval, *qnt)

    @staticmethod
    def getactual(mass, time, totmil, mode, _val, _eval, *qnt):
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
            positions along the marker trajectory.

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
                items[qnt.index(q)] = val()

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
        add("rho", lambda : _eval("rho"))
        add("psi", lambda : _eval("psi"))
        add("mileage", lambda : _val("mileage"))
        add("mass", lambda : mass)
        add("charge", lambda : _val("charge"))
        add("weight", lambda : _val("weight"))
        add("time", lambda : time)
        add("pncrid", lambda : _val("pncrid"))
        add("pncrdir", lambda : _val("pncdir"))
        add("connlen", lambda : totmil)

        firstmask = 0; lastmask=0
        if Orbits.GYROORBIT in mode:
            # Record the index of first masked array so that we can later append
            # all other arrays created here with GC mask (in hybrid mode)
            firstmask = len(items)
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
                add("vperp", lambda : np.sqrt(np.sum(pvecprt(mask)**2, axis=0)
                                              - vpar**2))
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
            lastmask = len(items)

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

        if Orbits.GYROORBIT in mode and Orbits.GUIDINGCENTER in mode:
            # Hybrid data, combine masked arrays
            nmasked = lastmask - firstmask
            for i in range(lastmask):
                if i < firstmask: continue
                items[i] = np.append(items[i], items[i+nmasked])
            for i in range(nmasked):
                del items[lastmask+i]

        for i in range(len(items)):
            if items[i] is None:
                raise ValueError("Unknown quantity in " + qnt[i])

        # Sort first by IDs and then by mileage
        ids  = _val("ids").v
        mile = _val("mileage").v
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
