"""Marker initial and final phase-space positions and related quantities.
"""
import numpy as np
import h5py
import unyt
import a5py.physlib as physlib

from .coreio import fileapi
from .coreio.treedata import DataContainer

class State(DataContainer):
    """State object.
    """

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

    def _endcond_check(self, bitarr, string):
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
                ec = getattr(self, ec)
            except AttributeError:
                raise ValueError("Unknown end condition: " + ec)

            i += 1
            if NOT:
                ec_non = ec_non | ec
            else:
                ec_yes = ec_yes | ec

        return (bitarr & ec_yes) == ec_yes and (bitarr & ec_non) == 0

    @property
    def ABORTED(self):
        """Marker simulation terminated in an error.
        """
        return 0x2

    @property
    def NONE(self):
        """No active end condition meaning the marker hasn't been simulated yet.
        """
        return 0x1

    @property
    def TLIM(self):
        """Simulation time limit or maximum mileage reached.
        """
        return 0x4

    @property
    def EMIN(self):
        """Minimum energy reached.
        """
        return 0x8

    @property
    def THERM(self):
        """Local thermal energy reached.
        """
        return 0x10

    @property
    def WALL():
        """Marker intersected a wall element.
        """
        return 0x20

    @property
    def RHOMIN(self):
        """Minimum radial coordinate (rho) reached.
        """
        return 0x40

    @property
    def RHOMAX(self):
        """Maximum radial coordinate (rho) reached.
        """
        return 0x80

    @property
    def POLMAX(self):
        """Maximum poloidal turns reached.
        """
        return 0x100

    @property
    def TORMAX(self):
        """Maximum toroidal turns reached.
        """
        return 0x200

    @property
    def CPUMAX(self):
        """Simulation for this marker exceeded the set CPU time.
        """
        return 0x400


    def read(self):
        """
        Read state data to dictionary.
        """
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid(),
                         self._path.split("/")[-1])


    def get(self, qnt, mode="gc"):
        """Return marker quantity.

        This function accesses the state data within the HDF5 file and uses that
        to evaluate the queried quantity.

        This method first checks whether the data can be found directly as
        a dataset in the HDF5 file.

        Parameters
        ----------
        qnt : str
            Name of the quantity.
        mode : {"prt", "gc"}, optional
            Is the quantity evaluated in particle or guiding center phase-space.

        Returns
        -------
        value : array_like
            The quantity as an array ordered by marker ID.
        """
        item = None
        #qnt  = physlib.alias.getalias(qnt)

        def read(q):
            """Read quantity from HDF5.
            """
            with self as h5:
                if q in h5:
                    return fileapi.read_data(h5, q)
            return None

        def evalprt(*q):
            """Evaluate input quantity at particle position.
            """
            return self._root._ascot.input_eval(
                read("rprt"), read("phiprt"), read("zprt"),
                read("time"), *q)

        def evalgc(*q):
            """Evaluate input quantity at guiding center position.
            """
            return self._root._ascot.input_eval(
                read("r"), read("phi"), read("z"),
                read("time"), *q)

        def bvecprt():
            """Evaluate particle magnetic field vector.
            """
            return unyt.unyt_array(evalprt("br", "bphi", "bz"))

        def bvecgc():
            """Evaluate guiding center magnetic field vector.
            """
            return unyt.unyt_array([read("br"), read("bphi"),
                                    read("bz")])

        def pvecprt():
            """Evaluate particle momentum vector.
            """
            return unyt.unyt_array([read("prprt"), read("pphiprt"),
                                    read("pzprt")])

        def pvecgc():
            """Evaluate guiding center momentum vector.
            """
            bnorm = np.sqrt(np.sum(bvecgc()**2,axis=0))
            pnorm = physlib.momentum_muppar(read("mass"), read("mu"),
                                            read("ppar"), bnorm)
            bhat  = bvecgc() / bnorm
            e1 = np.zeros(bhat.shape)
            e1[2,:] = 1
            e2 = np.cross(bhat.T, e1.T).T
            e1 = e2 / np.sqrt(np.sum(e2**2, axis=0))
            e2 = np.cross(bhat.T, e1.T).T
            perphat = -np.sin(read("zeta")) * e1 - np.cos(read("zeta")) * e2
            return bhat * read("ppar") \
                + perphat * np.sqrt(pnorm**2 - read("ppar")**2)

        def vvecprt():
            """Evaluate particle velocity vector.
            """
            return physlib.velocity_momentum(read("mass"), pvecprt())

        def vvecgc():
            """Evaluate guiding center velocity vector.
            """
            return physlib.velocity_momentum(read("mass"), pvecgc())

        # Dimensionless quantities common for particle and guiding-center
        if item is None and qnt in ["ids", "anum", "znum", "walltile",
                                    "endcond", "errormsg", "errorline",
                                    "errormod"]:
            item = read(qnt).v
            if qnt == "endcond":
                err = read("errormsg").v
                item = item << 2
                item[err > 0] = item[err > 0] & self.ABORTED
                item[item==0] = self.NONE

        # Dimensional quantities common for particle and guiding-center
        if qnt in ["mass", "charge", "time", "cputime", "mileage", "weight"]:
            item = read(qnt)
            item.convert_to_base("ascot")

        # Quantities that have to be evaluated separately (keep this in same
        # order for both and in the list method if new quantities are added)
        if item is None and mode == "prt":
            if   qnt == "r":
                item = read("rprt")
            elif qnt == "z":
                item = read("zprt")
            elif qnt == "phi":
                item = read("phiprt")
            elif qnt == "phimod":
                item = np.mod(read("phiprt"), 2 * np.pi * unyt.rad)
            elif qnt == "theta":
                axis = evalprt("axisr", "axisz")
                thetaprt = physlib.cart2pol(
                    read("rprt")-axis[0], read("zprt")-axis[1])[1]
                thetacum = np.floor(read("phi") / (2*np.pi * unyt.rad))
                theta = 2*np.pi * thetacum * unyt.rad + thetaprt
                item = np.mod(theta + np.pi*unyt.rad, 2*np.pi*unyt.rad)
            elif qnt == "thetamod":
                axis = evalprt("axisr", "axisz")
                theta = physlib.cart2pol(
                    read("rprt")-axis[0], read("zprt")-axis[1])[1]
                item = np.mod(theta + np.pi*unyt.rad, 2*np.pi*unyt.rad)
            elif qnt == "x":
                item = physlib.pol2cart(read("rprt"), read("phiprt"))[0]
            elif qnt == "y":
                item = physlib.pol2cart(read("rprt"), read("phiprt"))[1]
            elif qnt == "pr":
                item = read("prprt")
            elif qnt == "pz":
                item = read("pzprt")
            elif qnt == "pphi":
                item = read("pphiprt")
            elif qnt == "ppar":
                item = physlib.ppar_momentum(pvecprt(), bvecprt())
            elif qnt == "pperp":
                pnorm2 = np.sum( pvecprt()**2, axis=0 )
                pperp2 = physlib.ppar_momentum(pvecprt(), bvecprt())**2
                item = np.sqrt(pnorm2 - pperp2)
            elif qnt == "pnorm":
                item = np.sqrt( np.sum( pvecprt()**2, axis=0 ) )
            elif qnt == "vr":
                item = vvecprt()[0,:]
            elif qnt == "vz":
                item = vvecprt()[2,:]
            elif qnt == "vphi":
                item = vvecprt()[1,:]
            elif qnt == "vpar":
                item = physlib.vpar_momentum(read("mass"), pvecprt(), bvecprt())
            elif qnt == "vperp":
                vperp = physlib.vpar_momentum(
                    read("mass"), pvecprt(), bvecprt())
                item = np.sum(pvecprt()**2, axis=0)
            elif qnt == "vnorm":
                pnorm = np.sqrt(np.sum(pvecprt()**2, axis=0))
                item = physlib.velocity_momentum(read("mass"), pnorm)
            elif qnt == "br":
                item = evalprt("br")
            elif qnt == "bz":
                item = evalprt("bz")
            elif qnt == "bphi":
                item = evalprt("bphi")
            elif qnt == "bnorm":
                item = evalprt("bnorm")
            elif qnt == "ekin":
                item = physlib.energy_momentum(read("mass"), pvecprt())
            elif qnt == "pitch":
                item = physlib.pitch_momentum(pvecprt(), bvecprt())
            elif qnt == "mu":
                item = physlib.mu_momentum(read("mass"), pvecprt(), bvecprt())
            elif qnt == "zeta":
                item = read("zeta")
            elif qnt == "rho":
                item = evalprt("rho")
            elif qnt == "psi":
                item = evalprt("psi")
            elif qnt == "ptor":
                item = physlib.torcanangmom_momentum(
                    read("charge"), read("rprt"), pvecprt(), evalprt("psi"))
            if item is not None: item.convert_to_base("ascot")

        if item is None and mode == "gc":
            if   qnt == "r":
                item = read("r")
            elif qnt == "z":
                item = read("z")
            elif qnt == "phi":
                item = read("phi")
            elif qnt == "phimod":
                item = np.mod(read("phi"), 2 * np.pi * unyt.rad)
            elif qnt == "theta":
                item = read("theta")
            elif qnt == "thetamod":
                item = np.mod(read("theta"), 2 * np.pi * unyt.rad)
            elif qnt == "x":
                item = physlib.pol2cart(read("r"), read("phi"))[0]
            elif qnt == "y":
                item = physlib.pol2cart(read("r"), read("phi"))[1]
            elif qnt == "pr":
                item = pvecgc()[0,:]
            elif qnt == "pz":
                item = pvecgc()[2,:]
            elif qnt == "pphi":
                item = pvecgc()[1,:]
            elif qnt == "ppar":
                item = read("ppar")
            elif qnt == "pperp":
                pnorm = physlib.momentum_muppar(read("mass"), read("mu"),
                                                read("ppar"), bvecgc())
                item = np.sqrt(pnorm**2 - read("ppar")**2)
            elif qnt == "pnorm":
                item = physlib.momentum_muppar(read("mass"), read("mu"),
                                               read("ppar"), bvecgc())
            elif qnt == "vr":
                item = vvecgc()[0,:]
            elif qnt == "vz":
                item = vvecgc()[2,:]
            elif qnt == "vphi":
                item = vvecgc()[1,:]
            elif qnt == "vpar":
                item = physlib.vpar_muppar(read("mass"), read("mu"),
                                           read("ppar"), bvecgc())
            elif qnt == "vperp":
                pnorm = physlib.momentum_muppar(read("mass"), read("mu"),
                                                read("ppar"), bvecgc())
                pperp = np.sqrt(pnorm**2 - read("ppar")**2)
                gamma = physlib.gamma_momentum(read("mass"), pnorm)
                item = pperp / (gamma * read("mass"))
            elif qnt == "vnorm":
                pnorm = physlib.momentum_muppar(read("mass"), read("mu"),
                                                read("ppar"), bvecgc())
                gamma = physlib.gamma_momentum(read("mass"), pnorm)
                item = pnorm / (gamma * read("mass"))
            elif qnt == "br":
                item = read("br")
            elif qnt == "bz":
                item = read("bz")
            elif qnt == "bphi":
                item = read("bphi")
            elif qnt == "bnorm":
                item = np.sqrt(read("br")**2 + read("bphi")**2 + read("bz")**2)
            elif qnt == "ekin":
                item = physlib.energy_muppar(read("mass"), read("mu"),
                                             read("ppar"), bvecgc())
            elif qnt == "pitch":
                item = physlib.pitch_muppar(read("mass"), read("mu"),
                                            read("ppar"), bvecgc())
            elif qnt == "mu":
                item = read("mu")
            elif qnt == "zeta":
                item = read("zeta")
            elif qnt == "rho":
                item = read("rho")
            elif qnt == "psi":
                item = evalgc("psi")
            elif qnt == "ptor":
                item = physlib.torcanangmom_ppar(
                    read("charge"), read("r"), read("ppar"), bvecgc(),
                    evalgc("psi"))
            if item is not None: item.convert_to_base("ascot")

        if item is None:
            raise ValueError("Unknown quantity " + qnt)

        # Order by ID and return.
        idx  = (read("ids").v).argsort()
        return item[idx]

    def list(self):
        """List available quantities.
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
            "cputime":  "CPU time elapsed in simulation",
            "mileage":  "Laboratory time elapsed in simulation",
            "weight":   "How many physical particles a marker represents",
            "ids":      "Marker ID",
            "anum":     "Mass number",
            "znum":     "Charge number",
            "walltile": "Element ID if the marker hit wall",
            "endcond":  "Bitarray showing active end conditions",
            "errormsg": "Error message if marker was aborted",
            "errorline":"Line where (if) marker was aborted",
            "errormod": "Name of the module where (if) marker was aborted",
        }
        return out
