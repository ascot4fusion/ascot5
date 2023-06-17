"""State HDF5 IO module.
"""

import numpy as np
import h5py
import unyt

import a5py.physlib as physlib

from a5py.physlib.alias import getalias

from ._iohelpers.treedata import DataContainer
from ._iohelpers.fileapi import read_data

def write_hdf5(fn, run, name, data):
    """Write state data in HDF5 file.

    Parameters
    ----------
    fn : `str`
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


def read_hdf5(fn, qid, name):
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


class State(DataContainer):
    """State object.
    """

    def __init__(self, root, hdf5):
        """
        Initialize state object from given HDF5 file to given RunNode.
        """
        super().__init__(root, hdf5)


    def _endcond_check(self, bitarr, string):
        """Check if the binary end condition matches the human-readable.

        Parameters
        ----------
        bitarr : `np.uint32`
            Value of the end condition as it is stored in the HDF5 file.
        string : `str`
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
        match : `bool`
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


    def write(self, fn, run, name, data=None):
        """
        Write state data to HDF5 file.
        """
        if data is None:
            data = self.read()

        write_hdf5(fn, run, name, data)


    def __getitem__(self, key):
        """Return marker quantity.

        This function accesses the state data within the HDF5 file and uses that
        to evaluate the queried quantity.

        This method first checks whether the data can be found directly as
        a dataset in the HDF5 file.

        Parameters
        ----------
        key : `str`
            Name of the quantity.

        Returns
        -------
        value : np.array
            The quantity as an array ordered by marker ID.
        """

        # Wrapper for read data which opens the HDF5 file
        def read_dataw(data):
            with self as h5:
                data = read_data(h5, data)
            return data

        # Helper function that returns guiding center magnetic field vector
        def getbvec():
            return unyt.T * np.array(
                [read_dataw("br"),
                 read_dataw("bphi"),
                 read_dataw("bz")])

        # Helper function that returns particle momentum vector
        def getpvecprt():
            return unyt.kg * unyt.m / unyt.s * np.array(
                [read_dataw("prprt"),
                 read_dataw("pphiprt"),
                 read_dataw("pzprt")])

        # Helper function that evaluates ascotpy at guiding center position
        def evalapy(quantity):
            return self._root._ascot.input_eval(
                R   = read_dataw("r"),
                phi = read_dataw("phi").to("rad"),
                z   = read_dataw("z"),
                t   = read_dataw("time"),
                quantity = quantity
            )

        # Helper function that evaluates ascotpy at particle position
        def evalapyprt(quantity):
            return self._root._ascot.input_eval(
                R   = read_dataw("rprt"),
                phi = read_dataw("phiprt").to("rad"),
                z   = read_dataw("zprt"),
                t   = read_dataw("time"),
                quantity= quantity
            )

        # Helper function that returns particle magnetic field vector
        def getbvecprt():
            return unyt.T * np.array(
                [evalapyprt("br"),
                 evalapyprt("bphi"),
                 evalapyprt("bz")])

        # Get alias
        key  = getalias(key)
        item = None

        ## See if the field can be read directly and without conversions ##
        with self as h5:
            h5keys = list(h5.keys())
        if key in h5keys:
            item = read_dataw(key)

        if item is not None:
            pass

        ## Coordinates ##
        elif key  == "x":
            item = physlib.xcoord(
                r   = read_dataw("r"),
                phi = read_dataw("phi")
            )
        elif key == "xprt":
            item = physlib.xcoord(
                r   = read_dataw("rprt"),
                phi = read_dataw("phiprt")
            )
        elif key  == "y":
            item = physlib.ycoord(
                r   = read_dataw("r"),
                phi = read_dataw("phi")
            )
        elif key == "yprt":
            item = physlib.ycoord(
                r   = read_dataw("rprt"),
                phi = read_dataw("phiprt")
            )
        elif key == "phimod":
            item = np.mod(read_dataw("phi"), 2 * np.pi * unyt.rad)

        elif key == "phimodprt":
            item = np.mod(read_dataw("phiprt"), 2 * np.pi * unyt.rad)

        elif key == "thetamod":
            item = np.mod(read_dataw("theta"), 2 * np.pi * unyt.rad)

        ## Energy, gamma, and pitch ##
        elif key == "energy":
            item = physlib.energy_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "energyprt":
            item = physlib.energy_momentum(
                m = read_dataw("mass"),
                p = getpvecprt()
            )
        elif key == "gamma":
            item = physlib.gamma_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "gammaprt":
            item = physlib.gamma_momentum(
                m = read_dataw("mass"),
                p = getpvecprt()
            )
        elif key == "pitch":
            item = physlib.pitch_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "pitchprt":
            item = pitch_momentum(
                p = getpvecprt(),
                b = getbvecprt()
            )

        ## Velocity and momentum components, norms and mu ##
        elif key == "vpar":
            item = physlib.vpar_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "vparprt":
            item = physlib.vpar_momentum(
                m = read_dataw("mass"),
                p = getpvecprt(),
                b = getbvecprt()
            )

        elif key == "pparprt":
            item = physlib.ppar_momentum(
                p = getpvecprt(),
                b = getbvecprt()
            )

        elif key == "pnorm":
            item = physlib.momentum_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "pnormprt":
            item = getpvecprt()
            item = np.sqrt( np.sum( item**2, axis=1 ) )

        elif key == "vnorm":
            item = physlib.velocity_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "vnormprt":
            item = getpvecprt()
            item = np.sqrt( np.sum( item**2, axis=1 ) )
            item = physlib.velocity_momentum(
                m = read_dataw("mass"),
                p = item
            )
        elif key == "vrprt":
            item = physlib.velocity_momentum(
                m = read_dataw("mass"),
                p = getpvecprt()
            )[0,:]
        elif key == "vphiprt":
            item = physlib.velocity_momentum(
                m = read_dataw("mass"),
                p = getpvecprt()
            )[1,:]
        elif key == "vzprt":
            item = physlib.velocity_momentum(
                m = read_dataw("mass"),
                p = getpvecprt()
            )[2,:]
        elif key == "muprt":
            item = physlib.mu_momentum(
                m = read_dataw("mass"),
                p = getpvecprt(),
                b = getbvecprt()
            )

        elif key == "muprt":
            item = getbvec()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        ## Background quantities ##
        elif key == "bnorm":
            item = getbvec()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        elif key == "bnormprt":
            item = getbvecprt()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        elif key == "psi":
            item = evalapy("psi") * unyt.Wb

        elif key == "psiprt":
            item = evalapyprt("psi") * unyt.Wb

        elif key == "rhoprt":
            item = evalapyprt("rho") * unyt.dimensionless

        elif key == "ptor":
            item = physlib.torcanangmom_ppar(
                q    = read_dataw("charge"),
                r    = read_dataw("r"),
                ppar = read_dataw("ppar"),
                b    = getbvec(),
                psi  = evalapy("psi") * unyt.Wb
            )

        elif key == "ptorprt":
            item = physlib.torcanangmom_momentum(
                q   = read_dataw("charge"),
                r   = read_dataw("r"),
                p   = getpvecprt(),
                psi = evalapyprt("psi") * unyt.Wb
            )

        ## Boozer and MHD parameters ##
        elif key == "psi(bzr)":
            item = evalapy("psi (bzr)") * unyt.dimensionless

        elif key == "psi(bzr)prt":
            item = evalapyprt("psi (bzr)") * unyt.dimensionless

        elif key == "theta(bzr)":
            item = evalapy("theta") * unyt.rad

        elif key == "theta(bzr)prt":
            item = evalapyprt("theta") * unyt.rad

        elif key == "phi(bzr)":
            item = evalapy("zeta") * unyt.rad

        elif key == "phi(bzr)prt":
            item = evalapyprt("zeta") * unyt.rad

        elif key == "db/b(mhd)":
            item = evalapy("db/b") * unyt.T

        elif key == "db/b(mhd)prt":
            item = evalapyprt("db/b") * unyt.dimensionless

        elif key == "mhdepot":
            item = evalapy("phi") * unyt.V

        elif key == "mhdepotprt":
            item = evalapyprt("phi") * unyt.V

        elif key == "mhdalpha":
            item = evalapy("alpha") * unyt.m

        elif key == "mhdalphaprt":
            item = evalapyprt("alpha") * unyt.m

        if item is None:
            raise ValueError("Invalid query: " + key)

        # Strip units from fields to which they do not belong
        if key in ["ids", "endcond", "errormsg", "errorline", "errormod",
                   "walltile", "anum", "znum"]:
            item = item.v
        else:
            # Convert to ascot unit system.
            item.convert_to_base("ascot")

        # Dissect endcondition
        if key == "endcond":
            err = read_dataw("errormsg")
            item = item << 2
            item[err > 0] = item[err > 0] & endcondmod.getbin("aborted")
            item[item==0] = endcondmod.getbin("none")

        # Order by ID and return.
        idx  = (read_dataw("ids").v).argsort()
        return item[idx]
