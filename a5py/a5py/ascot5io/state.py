"""
State HDF5 IO module.

File: state.py
"""

import numpy as np
import h5py
import unyt
import warnings

from a5py.marker.alias import get as alias
import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.plot as plot
import a5py.marker.endcond as endcondmod
import a5py.physlib as physlib

from a5py.physlib.alias import getalias

from a5py.ascot5io.ascot5data import AscotData
from a5py.ascot5io.ascot5file import read_data

def write_hdf5(fn, run, name, data):
    """
    Write state data in HDF5 file.

    Args:
        fn : str <br>
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


class State(AscotData):
    """
    State object.
    """

    def __init__(self, hdf5, runnode):
        """
        Initialize state object from given HDF5 file to given RunNode.
        """
        self._runnode = runnode
        super().__init__(hdf5)


    def read(self):
        """
        Read state data to dictionary.
        """
        return read_hdf5(self._file, self.get_qid(), self._path.split("/")[-1])


    def write(self, fn, run, name, data=None):
        """
        Write state data to HDF5 file.
        """
        if data is None:
            data = self.read()

        write_hdf5(fn, run, name, data)


    def __getitem__(self, key):
        """
        Return queried marker quantity.

        This function accesses the state data within the HDF5 file and uses that
        to evaluate the queried quantity.

        This method first checks whether the data can be found directly as
        a dataset in the HDF5 file.

        Args:
            key : str <br>
                Name of the quantity.
        Returns:
            The quantity as an array ordered by marker ID.
        """
        # Init ascotpy object (no data is initialized yet)
        from a5py.ascotpy import Ascotpy
        a5 = Ascotpy(self._file)

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
            return a5.evaluate(
                R   = read_dataw("r"),
                phi = read_dataw("phi").to("rad"),
                z   = read_dataw("z"),
                t   = read_dataw("time"),
                quantity = quantity
            )

        # Helper function that evaluates ascotpy at particle position
        def evalapyprt(quantity):
            return a5.evaluate(
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
            a5.init(bfield=True)
            item = pitch_momentum(
                p = getpvecprt(),
                b = getbvecprt()
            )
            a5.free(bfield=True)

        ## Velocity and momentum components, norms and mu ##
        elif key == "vpar":
            item = physlib.vpar_muppar(
                m    = read_dataw("mass"),
                mu   = read_dataw("mu"),
                ppar = read_dataw("ppar"),
                b    = getbvec()
            )
        elif key == "vparprt":
            a5.init(bfield=True)
            item = physlib.vpar_momentum(
                m = read_dataw("mass"),
                p = getpvecprt(),
                b = getbvecprt()
            )
            a5.free(bfield=True)

        elif key == "pparprt":
            a5.init(bfield=True)
            item = physlib.ppar_momentum(
                p = getpvecprt(),
                b = getbvecprt()
            )
            a5.free(bfield=True)

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
            a5.init(bfield=True)
            item = physlib.mu_momentum(
                m = read_dataw("mass"),
                p = getpvecprt(),
                b = getbvecprt()
            )
            a5.free(bfield=True)

        elif key == "muprt":
            a5.init(bfield=True)
            item = getbvec()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        ## Background quantities ##
        elif key == "bnorm":
            item = getbvec()
            item = np.sqrt( np.sum( item**2, axis=0 ) )

        elif key == "bnormprt":
            a5.init(bfield=True)
            item = getbvecprt()
            item = np.sqrt( np.sum( item**2, axis=0 ) )
            a5.free(bfield=True)

        elif key == "psi":
            a5.init(bfield=True)
            item = evalapy("psi") * unyt.Wb
            a5.free(bfield=True)

        elif key == "psiprt":
            a5.init(bfield=True)
            item = evalapyprt("psi") * unyt.Wb
            a5.free(bfield=True)

        elif key == "rhoprt":
            a5.init(bfield=True)
            item = evalapyprt("rho") * unyt.dimensionless
            a5.free(bfield=True)

        elif key == "ptor":
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
                psi = evalapyprt("psi") * unyt.Wb
            )
            a5.free(bfield=True)

        ## Boozer and MHD parameters ##
        elif key == "psi(bzr)":
            a5.init(bfield=True, boozer=True)
            item = evalapy("psi (bzr)") * unyt.dimensionless
            a5.free(bfield=True, boozer=True)

        elif key == "psi(bzr)prt":
            a5.init(bfield=True, boozer=True)
            item = evalapyprt("psi (bzr)") * unyt.dimensionless
            a5.free(bfield=True, boozer=True)

        elif key == "theta(bzr)":
            a5.init(bfield=True, boozer=True)
            item = evalapy("theta") * unyt.rad
            a5.free(bfield=True, boozer=True)

        elif key == "theta(bzr)prt":
            a5.init(bfield=True, boozer=True)
            item = evalapyprt("theta") * unyt.rad
            a5.free(bfield=True, boozer=True)

        elif key == "phi(bzr)":
            a5.init(bfield=True, boozer=True)
            item = evalapy("zeta") * unyt.rad
            a5.free(bfield=True, boozer=True)

        elif key == "phi(bzr)prt":
            a5.init(bfield=True, boozer=True)
            item = evalapyprt("zeta") * unyt.rad
            a5.free(bfield=True, boozer=True)

        elif key == "db/b(mhd)":
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapy("db/b") * unyt.T
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key == "db/b(mhd)prt":
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapyprt("db/b") * unyt.dimensionless
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key == "mhdepot":
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapy("phi") * unyt.V
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key == "mhdepotprt":
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapyprt("phi") * unyt.V
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key == "mhdalpha":
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapy("alpha") * unyt.m
            a5.free(bfield=True, boozer=True, mhd=True)

        elif key == "mhdalphaprt":
            a5.init(bfield=True, boozer=True, mhd=True)
            item = evalapyprt("alpha") * unyt.m
            a5.free(bfield=True, boozer=True, mhd=True)

        if item is None:
            raise Exception("Invalid query: " + key)

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
            if hasattr(self._runnode, "endstate"):
                ec = self._runnode.endstate["endcond"]
            else:
                ec = self["endcond"]

            idx = np.logical_and( idx, ec==endcondmod.getbin(endcond) )

        if pncrid is not None:
            idx = np.logical_and(idx, self["pncrid"] == pncrid)

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self["id"], ids))

        val = val[idx]

        return val


    def listendconds(self):
        ec, counts = np.unique(self["endcond"], return_counts=True)
        endcond = []
        for e in ec:
            endcond.append(endcondmod.getname(e))

        return (endcond, counts)


    def scatter(self, x=None, y=None, z=None, c=None, endcond=None, pncrid=None,
                equal=False, log=False, axes=None):
        """
        Make scatter plot.

        """
        ids = self.get("id", endcond=endcond, pncrid=pncrid)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, endcond=endcond, pncrid=pncrid)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, pncrid=pncrid)

        zc = None
        if z is not None:
            zc = self.get(z, endcond=endcond, pncrid=pncrid)

        cc = None
        if c is not None:
            cc = self.get(c, endcond=endcond, pncrid=pncrid)

        if isinstance(log, tuple):
            if log[0]:
                xc = np.log10(np.absolute(xc))
            if log[1]:
                yc = np.log10(np.absolute(yc))
            if z is not None and log[2]:
                zc = np.log10(np.absolute(zc))
            if c is not None and log[2]:
                cc = np.log10(np.absolute(cc))
        elif log:
            xc = np.log10(np.absolute(xc))
            yc = np.log10(np.absolute(yc))
            if z is not None:
                zc = np.log10(np.absolute(zc))
            if c is not None:
                cc = np.log10(np.absolute(cc))

        plot.plot_scatter(x=xc, y=yc, z=zc, c=cc, equal=equal,
                          xlabel=x, ylabel=y, zlabel=z, axes=axes)


    def histogram(self, x, y=None, xbins=None, ybins=None, weight=False,
                  logx=False, logy=False, logscale=False, endcond=None,
                  axes=None):
        """
        Make histogram plot.
        """
        if y is None:
            # 1D plot

            if logx:
                xbins = np.linspace(np.log10(xbins[0]), np.log10(xbins[1]),
                                    xbins[2])
            else:
                xbins = np.linspace(xbins[0], xbins[1], xbins[2])

            weights=None
            if endcond is not None or not hasattr(self._runnode, "endstate"):
                xc = self.get(x, endcond=endcond)
                if weight:
                    weights = self.get("weight", endcond=endcond)
                if logx:
                    xc = np.log10(np.absolute(xc))
            else:
                xc = []
                if weight:
                    weights = []

                ecs, count = self._runnode.endstate.listendconds()
                for ec in ecs:
                    xc0 = self.get(x, endcond=ec)
                    if weight:
                        weights.append(self.get("weight", endcond=ec))
                    if logx:
                        xc0 = np.log10(np.absolute(xc0))

                    xc.append(xc0)

            plot.plot_histogram(x=xc, y=None, xbins=xbins, ybins=ybins,
                                weights=weights, logscale=logscale,
                                xlabel=x, ylabel=y, axes=axes)

        else:
            # 2D plot
            xc = self.get(x, endcond=endcond)
            yc = self.get(y, endcond=endcond)

            if logx:
                xc = np.log10(np.absolute(xc))
                xbins = np.linspace(np.log10(xbins[0]), np.log10(xbins[1]),
                                    xbins[2])
            else:
                xbins = np.linspace(xbins[0], xbins[1], xbins[2])

            if logy:
                yc = np.log10(np.absolute(yc))
                ybins = np.linspace(np.log10(ybins[0]), np.log10(ybins[1]),
                                    ybins[2])
            else:
                ybins = np.linspace(ybins[0], ybins[1], ybins[2])

            weights=None
            if weight:
                weights = self.get("weight", endcond=endcond)

            plot.plot_histogram(x=xc, y=yc, xbins=xbins, ybins=ybins,
                                weights=weights, logscale=logscale,
                                xlabel=x, ylabel=y, axes=axes)
