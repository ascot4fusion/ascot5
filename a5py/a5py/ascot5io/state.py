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

    N = data["id"].size

    with h5py.File(fn, "a") as f:
        g = f.create_group(gname)

        g.create_dataset("mass",      (N,1), data=data["mass"],      dtype="f8")
        g.create_dataset("time",      (N,1), data=data["time"],      dtype="f8")
        g.create_dataset("cputime",   (N,1), data=data["cputime"],   dtype="f8")
        g.create_dataset("weight",    (N,1), data=data["weight"],    dtype="f8")

        g.create_dataset("rprt",      (N,1), data=data["rprt"],      dtype="f8")
        g.create_dataset("phiprt",    (N,1), data=data["phiprt"],    dtype="f8")
        g.create_dataset("zprt",      (N,1), data=data["zprt"],      dtype="f8")
        g.create_dataset("prprt",     (N,1), data=data["prprt"],     dtype="f8")
        g.create_dataset("pphiprt",   (N,1), data=data["pphiprt"],   dtype="f8")
        g.create_dataset("pzprt",     (N,1), data=data["pzprt"],     dtype="f8")

        g.create_dataset("rgc",       (N,1), data=data["rgc"],       dtype="f8")
        g.create_dataset("phigc",     (N,1), data=data["phigc"],     dtype="f8")
        g.create_dataset("zgc",       (N,1), data=data["zgc"],       dtype="f8")
        g.create_dataset("mugc",      (N,1), data=data["mugc"],      dtype="f8")
        g.create_dataset("ppargc",    (N,1), data=data["ppargc"],    dtype="f8")
        g.create_dataset("zetagc",    (N,1), data=data["zetagc"],    dtype="f8")

        g.create_dataset("rhogc",     (N,1), data=data["rhogc"],     dtype="f8")
        g.create_dataset("thetagc",   (N,1), data=data["thetagc"],   dtype="f8")

        g.create_dataset("ids",       (N,1), data=data["ids"],      dtype="i8")
        g.create_dataset("walltile",  (N,1), data=data["walltile"],  dtype="i8")
        g.create_dataset("endcond",   (N,1), data=data["endcond"],   dtype="i8")
        g.create_dataset("anum",      (N,1), data=data["anum"],      dtype="i4")
        g.create_dataset("znum",      (N,1), data=data["znum"],      dtype="i4")
        g.create_dataset("charge",    (N,1), data=data["charge"],    dtype="i4")
        g.create_dataset("errormsg",  (N,1), data=data["errormsg"],  dtype="i4")
        g.create_dataset("errorline", (N,1), data=data["errorline"], dtype="i4")
        g.create_dataset("errormod",  (N,1), data=data["errormod"],  dtype="i4")

        g.create_dataset("brgc",      (N,1), data=data["brgc"],      dtype="f8")
        g.create_dataset("bphigc",    (N,1), data=data["bphigc"],    dtype="f8")
        g.create_dataset("bzgc",      (N,1), data=data["bzgc"],      dtype="f8")


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
                 read_dataw("bz")]).T

        # Helper function that returns particle momentum vector
        def getpvecprt():
            return unyt.kg * unyt.m / unyt.s * np.array(
                [read_dataw("prprt"),
                 read_dataw("pphiprt"),
                 read_dataw("pzprt")]).T

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
                 evalapyprt("bz")]).T

        # Get alias
        key  = getalias(key)
        item = None

        # See if the field can be read directly and without conversions
        with self as h5:
            h5keys = list(h5.keys())

        if key in h5keys:
            item = read_dataw(key)
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
            item = 0
        elif key == "thetamod":
            item = 0
        elif key == "phimodprt":
            item = 0
        elif key == "thetamodprt":
            item = 0
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
            item = 0
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
        elif key == "muprt":
            item = 0

        # Dissect endcondition
        if key == "endcond":
            err = read_dataw("errormsg")
            item = item << 2
            item[err > 0] = item[err > 0] & endcondmod.getbin("aborted")
            item[item==0] = endcondmod.getbin("none")

        # Convert to ascot unit system.
        item.convert_to_base("ascot")

        # Order by ID and return.
        idx  = read_dataw("ids").argsort()
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

            idx = np.logical_and( idx, ec == endcondmod.getbin(endcond) )

        if pncrid is not None:
            idx = np.logical_and(idx, self["pncrid"] == pncrid)

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
            xc = self.get(x, endcond=endcond, pncrid=pncrid, SI=False)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, pncrid=pncrid, SI=False)

        zc = None
        if z is not None:
            zc = self.get(z, endcond=endcond, pncrid=pncrid, SI=False)

        cc = None
        if c is not None:
            cc = self.get(c, endcond=endcond, pncrid=pncrid, SI=False)

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
                xc = self.get(x, endcond=endcond, SI=False)
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
                    xc0 = self.get(x, endcond=ec, SI=False)
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
            xc = self.get(x, endcond=endcond, SI=False)
            yc = self.get(y, endcond=endcond, SI=False)

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
