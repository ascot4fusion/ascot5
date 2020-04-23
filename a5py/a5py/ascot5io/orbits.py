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

        if "prprt" in data:
            g.create_dataset("rprt",    (N,1), data=data["rprt"],    dtype="f8")
            g.create_dataset("phiprt",  (N,1), data=data["phiprt"],  dtype="f8")
            g.create_dataset("zprt",    (N,1), data=data["zprt"],    dtype="f8")

            g.create_dataset("prprt",   (N,1), data=data["prprt"],   dtype="f8")
            g.create_dataset("pphiprt", (N,1), data=data["pphiprt"], dtype="f8")
            g.create_dataset("pzprt",   (N,1), data=data["pzprt"],   dtype="f8")

            g.create_dataset("rhoprt",   (N,1), data=data["rhoprt"],   dtype="f8")
            g.create_dataset("thetaprt", (N,1), data=data["thetaprt"], dtype="f8")
            g.create_dataset("brprt",    (N,1), data=data["brprt"],    dtype="f8")
            g.create_dataset("bphiprt",  (N,1), data=data["bphiprt"],  dtype="f8")
            g.create_dataset("bzprt",    (N,1), data=data["bzprt"],    dtype="f8")
        else:
            g.create_dataset("rgc",    (N,1), data=data["rgc"],    dtype="f8")
            g.create_dataset("phigc",  (N,1), data=data["phigc"],  dtype="f8")
            g.create_dataset("zgc",    (N,1), data=data["zgc"],    dtype="f8")

            g.create_dataset("rhogc",   (N,1), data=data["rhogc"],   dtype="f8")
            g.create_dataset("thetagc", (N,1), data=data["thetagc"], dtype="f8")
            g.create_dataset("brgc",    (N,1), data=data["brgc"],    dtype="f8")
            g.create_dataset("bphigc",  (N,1), data=data["bphigc"],  dtype="f8")
            g.create_dataset("bzgc",    (N,1), data=data["bzgc"],    dtype="f8")

        if "ppargc" in data:
            g.create_dataset("ppargc", (N,1), data=data["ppargc"], dtype="f8")
            g.create_dataset("mugc",   (N,1), data=data["mugc"],   dtype="f8")
            g.create_dataset("zetagc", (N,1), data=data["zetagc"], dtype="f8")


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

        # Helper function that returns magnetic field vector
        def getbvec():
            return unyt.T * np.array(
                [read_dataw("br"),
                 read_dataw("bphi"),
                 read_dataw("bz")]).T

        # Helper function that returns particle momentum vector
        def getpvecprt():
            return unyt.kg * unyt.m / unyt.s * np.array(
                [read_dataw("pr"),
                 read_dataw("pphi"),
                 read_dataw("pz")]).T

        # Helper function that evaluates ascotpy at marker position
        def evalapy(quantity):
            return a5.evaluate(
                R   = read_dataw("r"),
                phi = read_dataw("phi").to("rad"),
                z   = read_dataw("z"),
                t   = read_dataw("time"),
                quantity = quantity
            )

        # Gets mass for each datapoint from the endstate
        def mass():
            return self._read_from_inistate(
                "mass", read_dataw("ids")).ravel() * unyt.amu

        # Get alias and data field names in HDF5
        key  = getalias(key)
        item = None
        with self as h5:
            h5keys = list(h5.keys())

        # See if the field can be read directly and without conversions
        if key in h5keys:
            item = read_dataw(key)
        else:
            # See what type of data is stored
            if "mu" in h5keys:
                key += "gc"
            elif "charge" in h5keys:
                key += "prt"
            else:
                key += "fl"

        # Evaluate the quantity if direct read could not be done
        # New quantities can be added here. Use suffix "gc" for quantities that
        # can be evaluated from guiding-center data, "prt" for particle data
        # and "fl" for field line data.
        if item is not None:
            pass
        elif key in ["xgc", "xprt", "xfl"]:
            item = physlib.xcoord(
                r   = read_dataw("r"),
                phi = read_dataw("phi")
            )
        elif key in ["ygc", "yprt", "yfl"]:
            item = physlib.ycoord(
                r   = read_dataw("r"),
                phi = read_dataw("phi")
            )
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

        elif key == "pnormprt":
            item = getpvecprt()
            item = np.sqrt( np.sum( item**2, axis=1 ) )

        elif key == "vnormprt":
            item = getpvecprt()
            item = np.sqrt( np.sum( item**2, axis=1 ) )
            item = physlib.velocity_momentum(
                m = mass(),
                p = item
            )
        elif key == "vrprt":
            item = physlib.velocity_momentum(
                m = mass(),
                p = getpvecprt()
            )[:,0]
        elif key == "vphiprt":
            item = physlib.velocity_momentum(
                m = mass(),
                p = getpvecprt()
            )[:,1]
        elif key == "vzprt":
            item = physlib.velocity_momentum(
                m = mass(),
                p = getpvecprt()
            )[:,2]
        elif key == "muprt":
            item = physlib.mu_momentum(
                m = mass(),
                p = getpvecprt(),
                b = getbvec()
            )
        elif key in ["bnormgc", "bnormprt", "bnormfl"]:
            item = getbvec()
            item = np.sqrt( np.sum( item**2, axis=1 ) )

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


        # Dissect endcondition
        if key == "endcond":
            err = read_dataw("errormsg")
            item = item << 2
            item[err > 0] = item[err > 0] & endcondmod.getbin("aborted")
            item[item==0] = endcondmod.getbin("none")

        # Convert to ascot unit system.
        item.convert_to_base("ascot")

        # Order by id and time
        ids  = read_dataw("ids")
        time = read_dataw("time")
        idx  = np.lexsort((time, ids))

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
            with self as h5:
                ec = self._read_from_endstate("endcond", h5).ravel()

            idx = np.logical_and( idx, ec == endcondmod.getbin(endcond) )

        if pncrid is not None:
            idx = np.logical_and(idx, self["pncrid"] == pncrid)

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self["id"], ids))

        val = val[idx]

        return val


    def plot(self, x=None, y=None, z=None, endcond=None, pncrid=None,
             equal=False, log=False, axes=None):
        """
        Plot orbits as a continuous line.
        """
        ids = self.get("id", endcond=endcond)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, endcond=endcond, pncrid=pncrid, SI=False)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, pncrid=pncrid, SI=False)

        zc = None
        if z is not None:
            zc = self.get(z, endcond=endcond, pncrid=pncrid, SI=False)

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

        plot.plot_line(x=xc, y=yc, z=zc, ids=ids, equal=equal,
                       xlabel=x, ylabel=y, zlabel=z, axes=axes)


    def scatter(self, x=None, y=None, z=None, c=None, endcond=None, pncrid=None,
                sepid=False, equal=False, log=False, axes=None, prune=1, ids=None,
                markersize=5, **kwargs):
        """
        Make scatter plot.
        """
        ids = self.get("id", ids=ids, endcond=endcond, pncrid=pncrid)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, ids=ids, endcond=endcond, pncrid=pncrid, SI=False)

        yc = None
        if y is not None:
            yc = self.get(y, ids=ids, endcond=endcond, pncrid=pncrid, SI=False)

        zc = None
        if z is not None:
            zc = self.get(z, ids=ids, endcond=endcond, pncrid=pncrid, SI=False)

        cc = None
        if c is not None:
            cc = self.get(c, ids=ids, endcond=endcond, pncrid=pncrid, SI=False)

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

        plot.plot_scatter(x=xc, y=yc, z=zc, c=cc, equal=equal,
                          ids=ids, xlabel=x, ylabel=y, zlabel=z, axes=axes,
                          prune=prune, s=markersize, **kwargs)


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

        self.scatter(x=x, y=y, c=z, pncrid=pncrid, endcond=endcond, prune=prune,
                     sepid=sepid, log=log, equal=equal, axes=axes,
                     markersize=markersize, ids=ids)

        if y == "phimod":
            axes.set_ylim([0, 360])



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
