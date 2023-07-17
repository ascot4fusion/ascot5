"""Marker initial and final phase-space positions and related quantities.
"""
import numpy as np
import h5py
import unyt
import a5py.physlib as physlib

from .coreio import fileapi
from .coreio.treedata import DataContainer

class State(DataContainer):
    """Marker initial and final phase-space positions and related quantities.
    """

    # These end conditions are used internally so we define them as class
    # variables. All end conditions are defined as properties so that they
    # can have docstrings.
    _NONE    = 0x1
    _ABORTED = 0x2
    _TLIM    = 0x4
    _EMIN    = 0x8
    _THERM   = 0x10
    _WALL    = 0x20
    _RHOMIN  = 0x40
    _RHOMAX  = 0x80
    _POLMAX  = 0x100
    _TORMAX  = 0x200
    _CPUMAX  = 0x400
    _NEUTR   = 0x800
    _IONIZ   = 0x1000

    @property
    def ABORTED(self):
        """Marker simulation terminated in an error.
        """
        return State._ABORTED

    @property
    def NONE(self):
        """No active end condition meaning the marker hasn't been simulated yet.
        """
        return State._NONE

    @property
    def TLIM(self):
        """Simulation time limit or maximum mileage reached.
        """
        return State._TLIM

    @property
    def EMIN(self):
        """Minimum energy reached.
        """
        return State._EMIN

    @property
    def THERM(self):
        """Local thermal energy reached.
        """
        return State._THERM

    @property
    def WALL():
        """Marker intersected a wall element.
        """
        return State._WALL

    @property
    def RHOMIN(self):
        """Minimum radial coordinate (rho) reached.
        """
        return State._RHOMIN

    @property
    def RHOMAX(self):
        """Maximum radial coordinate (rho) reached.
        """
        return _RHOMAX

    @property
    def POLMAX(self):
        """Maximum poloidal turns reached.
        """
        return State._POLMAX

    @property
    def TORMAX(self):
        """Maximum toroidal turns reached.
        """
        return State._TORMAX

    @property
    def CPUMAX(self):
        """Simulation for this marker exceeded the set CPU time.
        """
        return State._CPUMAX

    @property
    def NEUTR(self):
        """Ion marker neutralized.
        """
        return State._NEUTR

    @property
    def IONIZ(self):
        """Neutral marker ionized.
        """
        return State._IONIZ

    def write_hdf5(self):
        """Write state data in HDF5 file.

        Parameters
        ----------
        fn : str
        Full path to the HDF5 file.
        """

        gname = "results/" + run + "/" + name

        fields = [("mass",      "amu"),
                  ("time",      "s"),
                  ("cputime",   "s"),
                  ("weight",    "markers/s"),
                  ("rprt",      "m"),
                  ("phiprt",    "deg"),
                  ("zprt",      "m"),
                  ("prprt",     "kg*m/s"),
                  ("pphiprt",   "kg*m/s"),
                  ("pzprt",     "kg*m/s"),
                  ("r",         "m"),
                  ("phi",       "deg"),
                  ("z",         "m"),
                  ("mu",        "eV/T"),
                  ("ppar",      "kg*m/s"),
                  ("zeta",      "rad"),
                  ("rho",       "1"),
                  ("theta",     "deg"),
                  ("ids",       "1"),
                  ("walltile",  "1"),
                  ("endcond",   "1"),
                  ("anum",      "1"),
                  ("znum",      "1"),
                  ("charge",    "e"),
                  ("errormsg",  "1"),
                  ("errorline", "1"),
                  ("errormod",  "1"),
                  ("br",        "T"),
                  ("bphi",      "T"),
                  ("bz",        "T")]

        N = data["ids"].size

        with h5py.File(fn, "a") as f:
            g = f.create_group(gname)

            for field in fields:
                ds = g.create_dataset(field[0], (N,1), data=data[field[0]],
                                      dtype="f8")
                ds.attrs["unit"] = field[1]


    def read_hdf5(self, fn, qid, name):
        """
        Read state output from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.
        name : str <br>
            Name of the data to read, e.g. "inistate", "endstate", "distrho5d"

        Returns:
        Dictionary containing input data.
        """

        path = "results/run_" + qid + "/" + name

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        return out

    def read(self):
        """
        Read state data to dictionary.
        """
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid(),
                         self._path.split("/")[-1])

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
        bvecgc = lambda : unyt.unyt_array([_val("br"), _val("bphi"),
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
                vperp = physlib.vpar_momentum(
                    _val("mass"), pvecprt(), bvecprt())
                add("vperp", lambda : np.sum(pvecprt()**2, axis=0))
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
    def endcond_check(bitarr, string):
        """Check if the binary end condition matches the human-readable.

        Parameters
        ----------
        bitarr : :obj:`np.uint32`
            Value of the end condition as it is stored in the HDF5 file.
        string : str
            Human-readable end condition (case-insensitive).

            Markers may have multiple end conditions active simultaneously.
            If just the name of the end condition e.g. "MAXPOL" is passed,
            then all markers with the ``MAXPOL`` end condition are returned.

            If the end cond is preceded by "NOT", e.g. "NOT MAXPOL", then
            markers that don't have that end condition are returned.

            Passing multiple end conditions returns markers that have all listed
            end conditions active, e.g. "MAXPOL MAXTOR" returns markers that
            have both ``MAXPOL`` and ``MAXTOR`` active simultaneously.

        Returns
        -------
        match : bool
            True if the two representations of end conditions match.
        """
        endconds = string.upper().split()
        ec_yes = np.array(0, dtype=np.uint32)
        ec_non = np.array(0, dtype=np.uint32)

        i = 0
        while i < len(endconds):
            ec = endconds[i]

            NOT = False
            if ec == "NOT":
                NOT = True
                i += 1
                ec = endconds[i]

            try:
                ec = getattr(State, "_" + ec)
            except AttributeError:
                raise ValueError("Unknown end condition: " + ec)

            i += 1
            if NOT:
                ec_non = ec_non | ec
            else:
                ec_yes = ec_yes | ec

        return (bitarr & ec_yes) == ec_yes and (bitarr & ec_non) == 0

    @staticmethod
    def endcond_tostring(bitarr):
        """Convert end condition bitarray to human readable.

        Parameters
        ----------
        bitarr : :obj:`np.uint32`
            Value of the end condition as it is stored in the HDF5 file.

        Returns
        -------
        string : str
            End condition in a human readable format.
        """
        endcond = ["NONE", "ABORTED", "TLIM", "EMIN", "THERM", "WALL", "RHOMIN",
                   "RHOMAX", "POLMAX", "TORMAX", "CPUMAX"]
        string = ""
        for ec in endcond:
            if bitarr & getattr(State, "_" + ec):
                string += ec + " and "
        return string[:-5] # remove last "and"

    @staticmethod
    def listqnts():
        """List available quantities.
        """
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
