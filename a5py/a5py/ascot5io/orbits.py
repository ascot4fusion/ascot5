"""
Orbit data IO module.

File: orbits.py
"""
import numpy as np
import h5py
import unyt

import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.endcond as endcondmod
import a5py.marker.plot as plot
import a5py.physlib as physlib

from a5py.physlib.alias import getalias

from a5py.ascot5io.ascot5data import AscotData
from a5py.ascot5io.ascot5file import read_data

def read_hdf5(fn, qid):
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


def write_hdf5(fn, run, data, desc=None):
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


class Orbits(AscotData):
    """
    Object representing orbit data.
    """

    def __init__(self, hdf5, runnode):
        """
        Initialize orbit object from given HDF5 file to given RunNode.
        """
        self._runnode = runnode
        super().__init__(hdf5)


    def read(self):
        """
        Read orbit data to dictionary.
        """
        return read_hdf5(self._file, self.get_qid())


    def write(self, fn, run, data=None):
        """
        Write orbit data to HDF5 file.
        """
        if data is None:
            data = self.read()

        write_hdf5(fn, run, data)


    def _eval(self, key, simmode, a5):
        """
        Evaluate a quantity.
        """
        with self as h5:
            mask = read_data(h5, "simmode").v == simmode

        # Wrapper for read data which opens the HDF5 file
        def read_dataw(data):
            with self as h5:
                data = read_data(h5, data)
            return data[mask]

        # Helper function that returns magnetic field vector
        def getbvec():
            return unyt.T * np.array(
                [read_dataw("br"),
                 read_dataw("bphi"),
                 read_dataw("bz")])

        # Helper function that returns particle momentum vector
        def getpvecprt():
            return unyt.kg * unyt.m / unyt.s * np.array(
                [read_dataw("pr"),
                 read_dataw("pphi"),
                 read_dataw("pz")])

        # Helper function that evaluates ascotpy at marker position
        def evalapy(quantity):
            return a5.evaluate(
                R   = read_dataw("r"),
                phi = read_dataw("phi").to("rad"),
                z   = read_dataw("z"),
                t   = read_dataw("mileage") + initime(),
                quantity = quantity
            )

        # Gets mass and initime for each datapoint from the inistate
        def mass():
            m = self._read_from_inistate(
                "mass", read_dataw("ids")).ravel() * unyt.amu
            return m[mask]

        def initime():
            m = self._read_from_inistate(
                "time", read_dataw("ids")).ravel() * unyt.s
            return m[mask]

        item = None
        with self as h5:
            h5keys = list(h5.keys())

        # See if the field can be read directly and without conversions
        if key in h5keys:
            item = read_dataw(key)

        # Mass and time are read from inistate
        elif key == "mass":
            item = mass()
        elif key == "time":
            item = initime() + read_dataw("mileage")
        elif key == "mileagerev":
            item = read_dataw("mileage")
            item = item[0] - item
        else:
            if simmode == 1:
                key += "prt"
            elif simmode == 2:
                key += "gc"
            else:
                key += "fl"

        # New quantities can be added here. Use suffix "gc" for quantities that
        # can be evaluated from guiding-center data, "prt" for particle data
        # and "fl" for field line data.
        if item is not None:
            pass

        ## Coordinates ##
        if key in ["xgc", "xprt", "xfl"]:
            item = physlib.xcoord(
                r   = read_dataw("r"),
                phi = read_dataw("phi")
            )
        elif key in ["ygc", "yprt", "yfl"]:
            item = physlib.ycoord(
                r   = read_dataw("r"),
                phi = read_dataw("phi")
            )
        elif key in ["phimodgc", "phimodprt", "phimodfl"]:
            item = np.mod(read_dataw("phi"), 2 * np.pi * unyt.rad)

        elif key in ["thetamodgc", "thetamodprt", "thetamodfl"]:
            item = np.mod(read_dataw("theta"), 2 * np.pi * unyt.rad)

        ## Energy, gamma, and pitch ##
        elif key == "energygc":
            item = physlib.energy_muppar(
                m    = mass(),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "energyprt":
            item = physlib.energy_momentum(
                m = mass(),
                p = getpvecprt()
            )
        elif key == "gammagc":
            item = physlib.gamma_muppar(
                m    = mass(),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "gammaprt":
            item = physlib.gamma_momentum(
                m = mass(),
                p = getpvecprt()
            )
        elif key == "pitchgc":
            item = physlib.pitch_muppar(
                m    = mass(),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "pitchprt":
            item = pitch_momentum(
                p = getpvecprt(),
                b = getbvec()
            )

        ## Velocity and momentum components, norms and mu ##
        elif key == "vpargc":
            item = physlib.vpar_muppar(
                m    = mass(),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "vparprt":
            item = physlib.vpar_momentum(
                m = mass(),
                p = getpvecprt(),
                b = getbvec()
            )
        elif key == "pparprt":
            item = physlib.ppar_momentum(
                p = getpvecprt(),
                b = getbvec()
            )
        elif key == "pnormgc":
            item = physlib.momentum_muppar(
                m    = mass(),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "pnormprt":
            item = getpvecprt()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        elif key == "vnormgc":
            item = physlib.velocity_muppar(
                m    = mass(),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "vnormprt":
            item = getpvecprt()
            item = np.sqrt( np.sum( item**2, axis=0 ) )
            item = physlib.velocity_momentum(
                m = mass(),
                p = item
            )
        elif key == "vrprt":
            item = physlib.velocity_momentum(
                m = mass(),
                p = getpvecprt()
            )[0,:]
        elif key == "vphiprt":
            item = physlib.velocity_momentum(
                m = mass(),
                p = getpvecprt()
            )[1,:]
        elif key == "vzprt":
            item = physlib.velocity_momentum(
                m = mass(),
                p = getpvecprt()
            )[2,:]
        elif key == "muprt":
            item = physlib.mu_momentum(
                m = mass(),
                p = getpvecprt(),
                b = getbvec()
            )

        ## Background quantities ##
        elif key in ["bnormgc", "bnormprt", "bnormfl"]:
            item = getbvec()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        elif key in ["psigc", "psiprt", "psifl"]:
            a5.init(bfield=True)
            item = evalapy("psi") * unyt.Wb
            a5.free(bfield=True)

        elif key in ["rhogc", "rhoprt", "rhofl"]:
            a5.init(bfield=True)
            item = evalapy("rho") * unyt.dimensionless
            a5.free(bfield=True)

        elif key == "ptorgc":
            a5.init(bfield=True)
            item = physlib.torcanangmom_ppar(
                q    = read_dataw("charge"),
                r    = read_dataw("r"),
                ppar = read_dataw("ppar"),
                b    = getbvec(),
                psi  = evalapy("psi") * unyt.Wb
            )
            a5.free(bfield=True)

        elif key == "ptorprt":
            a5.init(bfield=True)
            item = physlib.torcanangmom_momentum(
                q   = read_dataw("charge"),
                r   = read_dataw("r"),
                p   = getpvecprt(),
                psi = evalapy("psi") * unyt.Wb
            )
            a5.free(bfield=True)

        ## Boozer and MHD parameters ##
        elif key in ["psi(bzr)gc", "psi(bzr)prt", "psi(bzr)fl"]:
            a5.init(bfield=True, boozer=True)
            item = evalapy("psi (bzr)") * unyt.dimensionless
            a5.free(bfield=True, boozer=True)

        elif key in ["theta(bzr)gc", "theta(bzr)prt", "theta(bzr)fl"]:
            a5.init(bfield=True, boozer=True)
            item = evalapy("theta") * unyt.rad
            a5.free(bfield=True, boozer=True)

        elif key in ["phi(bzr)gc", "phi(bzr)prt", "phi(bzr)fl"]:
            a5.init(bfield=True, boozer=True)
            item = evalapy("zeta") * unyt.rad
            a5.free(bfield=True, boozer=True)

        elif key in ["db/b(mhd)gc", "db/b(mhd)prt", "db/b(mhd)fl"]:
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapy("db/b") * unyt.dimensionless
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key in ["mhdepotgc", "mhdepotprt", "mhdepotfl"]:
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapy("phi") * unyt.V
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key in ["mhdalphagc", "mhdalphaprt", "mhdalphafl"]:
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapy("alpha") * unyt.m
            a5.free(bfield=True, boozer=True, mhd=True)

        return item

    
    def __getitem__(self, key):
        """
        Return queried quantity.

        The quantity is returned as a single numpy array ordered by id and time.

        The quantity will be in the same picture as the stored data is, e.g. if
        one queries for magnetic momentum mu, mu will be in guiding center
        picture if the data corresponds to guiding center data, particle picture
        for particle data, and none for field line data.

        Args:
            key : str <br>
                Name of the quantity.
        Returns:
            The quantity ordered by id and time.
        """

        # Init ascotpy object (no data is initialized yet)
        from a5py.ascotpy import Ascotpy
        a5 = Ascotpy(self._file)

        # Wrapper for read data which opens the HDF5 file
        def read_dataw(data):
            with self as h5:
                data = read_data(h5, data)
            return data

        # Get alias
        key  = getalias(key)

        # Get simmode in each data point which is used as a mask
        simmode = read_dataw("simmode")
        sids    = np.unique(simmode)

        item = self._eval(key, sids[0], a5)
        item.convert_to_base("ascot")

        # Deal with the hybrid case
        if sids.size > 1:
            item0 = np.zeros(simmode.shape) * item.unit_quantity
            item0[simmode == sids[0]] = item
            item = self._eval(key, sids[1], a5)
            item.convert_to_base("ascot")

            item0[simmode == sids[1]] = item
            item = item0


        if item is None:
            raise Exception("Invalid query: " + key)

        # Strip units from fields to which they do not belong
        if key in ["ids"]:
            item = item.v

        # Dissect endcondition
        if key == "endcond":
            err = read_dataw("errormsg")
            item = item << 2
            item[err > 0] = item[err > 0] & endcondmod.getbin("aborted")
            item[item==0] = endcondmod.getbin("none")

        # Order by id and mileage
        ids  = read_dataw("ids").v
        mile = read_dataw("mileage")
        idx  = np.lexsort((mile, ids))

        return item[idx]


    def get(self, key, ids=None, endcond=None, pncrid=None):
        """
        Same as __getitem__ but with option to filter which points are returned.

        Args:
            key : str <br>
                Name of the quantity (see alias.py for a complete list).
            ids : int, array_like, optional <br>
                Id or a list of ids whose data points are returned.
            endcond : str, array_like, optional <br>
                Endcond of those  markers which are returned.
            pncrid : str, array_like, optional <br>
                Poincare ID of those  markers which are returned.
        Returns:
            The quantity.
        """
        val = self[key]

        idx = np.ones(val.shape, dtype=bool)

        if endcond is not None:
            ec = self._read_from_endstate("endcond", self["ids"]).ravel()
            idx = np.logical_and( idx, ec == endcondmod.getbin(endcond) )

        if pncrid is not None:
            idx = np.logical_and(idx, self["pncrid"] == pncrid)

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self["id"], ids))

        val = val[idx]

        return val


    def get_datatype(self):
        """
        Get string describing what data this object contains (prt, gc, fl).
        """
        with self as h5:
            simmode = np.unique(read_data(h5, "simmode"))

        if 1 in simmode and 2 in simmode:
            return "hybrid"
        if 1 in simmode:
            return "prt"
        if 2 in simmode:
            return "gc"
        if 3 in simmode:
            return "fl"


    def plot(self, x=None, y=None, z=None, endcond=None, pncrid=None,
             equal=False, log=False, axes=None, ids=None, **kwargs):
        """
        Plot orbits as a continuous line.
        """
        ids0 = ids
        ids = self.get("id", endcond=endcond, ids=ids0)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, endcond=endcond, pncrid=pncrid, ids=ids0)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, pncrid=pncrid, ids=ids0)

        zc = None
        if z is not None:
            zc = self.get(z, endcond=endcond, pncrid=pncrid, ids=ids0)

        if isinstance(log, tuple):
            if log[0]:
                xc = np.log10(np.absolute(xc))
            if log[1]:
                yc = np.log10(np.absolute(yc))
            if z is not None and log[2]:
                zc = np.log10(np.absolute(zc))
        elif log:
            xc = np.log10(np.absolute(xc))
            yc = np.log10(np.absolute(yc))
            if z is not None:
                zc = np.log10(np.absolute(zc))

        axes = plot.plot_line(x=xc, y=yc, z=zc, ids=ids, equal=equal,
                              xlabel=x, ylabel=y, zlabel=z, axes=axes, **kwargs)

        return axes


    def scatter(self, x=None, y=None, z=None, c=None, endcond=None, pncrid=None,
                sepid=False, equal=False, log=False, axes=None, prune=1, ids=None,
                markersize=5, **kwargs):
        """
        Make scatter plot.
        """
        ids = self.get("id", ids=ids, endcond=endcond, pncrid=pncrid)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, ids=ids, endcond=endcond, pncrid=pncrid)

        yc = None
        if y is not None:
            yc = self.get(y, ids=ids, endcond=endcond, pncrid=pncrid)

        zc = None
        if z is not None:
            zc = self.get(z, ids=ids, endcond=endcond, pncrid=pncrid)

        cc = None
        if c is not None:
            cc = self.get(c, ids=ids, endcond=endcond, pncrid=pncrid)

        if not sepid:
            ids = None

        if isinstance(log, tuple):
            if log[0]:
                xc = np.log10(np.absolute(xc))
            if log[1]:
                yc = np.log10(np.absolute(yc))
            if z is not None and log[2]:
                zc = np.log10(np.absolute(zc))
            if c is not None and log[3]:
                cc = np.log10(np.absolute(cc))
        elif log:
            xc = np.log10(np.absolute(xc))
            yc = np.log10(np.absolute(yc))
            if z is not None:
                zc = np.log10(np.absolute(zc))
            if c is not None:
                cc = np.log10(np.absolute(cc))

        axes = plot.plot_scatter(x=xc, y=yc, z=zc, c=cc, equal=equal,
                                 ids=ids, xlabel=x, ylabel=y, zlabel=z,
                                 axes=axes, prune=prune, s=markersize, **kwargs)

        return axes


    def poincare(self, *args, log=False, endcond=None, equal=False,
                 prune=1, ids=None, markersize=5, axes=None):
        """
        Make a Poincare plot.
        """

        z   = None
        sepid = False
        if len(args) == 1:
            x = "rho"
            y = "phimod"
            pncrid = args[0]
            sepid = True

        elif len(args) == 3:
            x = args[0]
            y = args[1]
            pncrid = args[2]
            sepid = True

        if len(args) == 4:
            x = args[0]
            y = args[1]
            z = args[2]
            pncrid = args[3]
            if log:
                log = (0, 0, 1, 0)

        if x == "R" and y == "z":
            equal = True

        axes = self.scatter(x=x, y=y, c=z, pncrid=pncrid, endcond=endcond,
                            prune=prune, sepid=sepid, log=log, equal=equal,
                            axes=axes, markersize=markersize, ids=ids)

        if y == "phimod":
            axes.set_ylim([0, 360])

        return axes


    def _read_from_inistate(self, key, ids):
        isval = self._runnode.inistate[key]
        isid  = self._runnode.inistate["ids"]
        f     = lambda x: isval[isid == x]
        return np.array([f(x) for x in ids])


    def _read_from_endstate(self, key, ids):
        esval = self._runnode.endstate[key]
        esid  = self._runnode.endstate["ids"]
        f     = lambda x: esval[esid == x]
        return np.array([f(x) for x in ids])
