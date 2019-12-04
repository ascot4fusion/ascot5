"""
State HDF5 IO module.

File: state.py
"""

import numpy as np
import h5py
import warnings

from a5py.marker.alias import get as alias
import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.plot as plot
import a5py.marker.endcond as endcondmod

from a5py.marker.alias import get as alias

from a5py.ascot5io.ascot5data import AscotData

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
        g.create_dataset("vr",        (N,1), data=data["vr"],        dtype="f8")
        g.create_dataset("vphi",      (N,1), data=data["vphi"],      dtype="f8")
        g.create_dataset("vz",        (N,1), data=data["vz"],        dtype="f8")

        g.create_dataset("r",         (N,1), data=data["r"],         dtype="f8")
        g.create_dataset("phi",       (N,1), data=data["phi"],       dtype="f8")
        g.create_dataset("z",         (N,1), data=data["z"],         dtype="f8")
        g.create_dataset("mu",        (N,1), data=data["mu"],        dtype="f8")
        g.create_dataset("vpar",      (N,1), data=data["vpar"],      dtype="f8")
        g.create_dataset("zeta",      (N,1), data=data["zeta"],      dtype="f8")

        g.create_dataset("rho",       (N,1), data=data["rho"],       dtype="f8")
        g.create_dataset("theta",     (N,1), data=data["theta"],     dtype="f8")

        g.create_dataset("id",        (N,1), data=data["id"],        dtype="i8")
        g.create_dataset("walltile",  (N,1), data=data["walltile"],  dtype="i8")
        g.create_dataset("endcond",   (N,1), data=data["endcond"],   dtype="i8")
        g.create_dataset("anum",      (N,1), data=data["anum"],      dtype="i4")
        g.create_dataset("znum",      (N,1), data=data["znum"],      dtype="i4")
        g.create_dataset("charge",    (N,1), data=data["charge"],    dtype="i4")
        g.create_dataset("errormsg",  (N,1), data=data["errormsg"],  dtype="i4")
        g.create_dataset("errorline", (N,1), data=data["errorline"], dtype="i4")
        g.create_dataset("errormod",  (N,1), data=data["errormod"],  dtype="i4")

        g.create_dataset("br",        (N,1), data=data["br"],        dtype="f8")
        g.create_dataset("bphi",      (N,1), data=data["bphi"],      dtype="f8")
        g.create_dataset("bz",        (N,1), data=data["bz"],        dtype="f8")
        g.create_dataset("brdr",      (N,1), data=data["brdr"],      dtype="f8")
        g.create_dataset("brdphi",    (N,1), data=data["brdphi"],    dtype="f8")
        g.create_dataset("brdz",      (N,1), data=data["brdz"],      dtype="f8")
        g.create_dataset("bphidr",    (N,1), data=data["bphidr"],    dtype="f8")
        g.create_dataset("bphidphi",  (N,1), data=data["bphidphi"],  dtype="f8")
        g.create_dataset("bphidz",    (N,1), data=data["bphidz"],    dtype="f8")
        g.create_dataset("bzdr",      (N,1), data=data["bzdr"],      dtype="f8")
        g.create_dataset("bzdphi",    (N,1), data=data["bzdphi"],    dtype="f8")
        g.create_dataset("bzdz",      (N,1), data=data["bzdz"],      dtype="f8")


def read_hdf5(fn, qid, name):
    """
    Read state output from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

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
        Return queried quantity.

        The quantity is returned as a single numpy array ordered by id and time.
        Internally, this method first sees if the quantity can be read directly
        from HDF5. If not, then it tries to see if it is present in endstate and
        can be copied from there (e.g. mass). If not, then quantity is evaluated
        by first determining if the stored orbit type is field line (has no
        charge), guiding center (has magnetic moment), or particle.

        Args:
            key : str <br>
                Name of the quantity (see alias.py for a complete list).
        Returns:
            The quantity in SI units ordered by id and time.
        """
        with self as h5:

            key  = alias(key)
            item = None

            # See if the field can be read directly and without conversions
            h5keys = list(h5.keys())
            if key in h5keys:
                item = h5[key][:]

                # Unit conversions
                if key == "charge":
                    f    = lambda x: interpret.charge_C(x)
                    item = np.array([f(x) for x in item]).ravel()
                if key == "mu":
                    f    = lambda x: interpret.energy_J(x)
                    item = np.array([f(x) for x in item]).ravel()
                if key == "phi":
                    item = item * np.pi/180
                if key == "mass":
                    f    = lambda x: interpret.mass_kg(x)
                    item = np.array([f(x) for x in item]).ravel()
                if key == "endcond":
                    err = h5["errormsg"][:]
                    item = item << 2
                    item[err > 0] = item[err > 0] & endcondmod.getbin("aborted")
                    item[item==0] = endcondmod.getbin("none")

            if item is None:

                # Convert guiding-center quantities to SI units
                f      = lambda x: interpret.mass_kg(x)
                mass   = np.array([f(x) for x in h5["mass"][:]]).ravel()
                f      = lambda x: interpret.charge_C(x)
                charge = np.array([f(x) for x in h5["charge"][:]]).ravel()
                f      = lambda x: interpret.energy_J(x)
                mu     = np.array([f(x) for x in h5["mu"][:]]).ravel()
                phi    = h5["phi"][:] * np.pi/180

                try:
                    item = marker.eval_guidingcenter(
                        key, mass=mass, charge=charge,
                        R=h5["r"][:], phi=phi, z=h5["z"][:],
                        mu=mu, vpar=h5["vpar"][:],
                        theta=h5["theta"][:],
                        BR=h5["br"][:], Bphi=h5["bphi"][:],
                        Bz=h5["bz"][:])
                except ValueError:
                    pass

                if item is None:
                    phi = h5["phiprt"][:] * np.pi/180

                    # The magnetic field values in the HDF5 file are for the
                    # guiding center position. Those need to be evaluated in
                    # particle position using Ascotpy.
                    from a5py.ascotpy import Ascotpy
                    a5 = None
                    try:
                        a5 = Ascotpy(self._file)
                        a5.init(bfield=self._runnode.bfield.get_qid())

                        br   = a5.evaluate(h5["rprt"][:], phi=phi,
                                           z=h5["zprt"][:],
                                           t=h5["time"][:], quantity="br")
                        bphi = a5.evaluate(h5["rprt"][:], phi=phi,
                                           z=h5["zprt"][:],
                                           t=h5["time"][:], quantity="bphi")
                        bz   = a5.evaluate(h5["rprt"][:], phi=phi,
                                           z=h5["zprt"][:],
                                           t=h5["time"][:], quantity="bz")
                    except:
                        warnings.warn("Ascotpy initialization failed.\n" +
                                      "Using magnetic field values at gc " +
                                      "position to evaluate particle data.")
                        br   = h5["br"][:]
                        bphi = h5["bphi"][:]
                        bz   = h5["bz"][:]

                    if key == "brprt":
                        item = br
                    elif key == "bphiprt":
                        item = bphi
                    elif key == "bzprt":
                        item = bz
                    else:
                        item = marker.eval_particle(
                            key, mass=mass, charge=charge,
                            R=h5["rprt"][:], phi=phi, z=h5["zprt"][:],
                            vR=h5["vr"][:], vphi=h5["vphi"][:], vz=h5["vz"][:],
                            BR=br, Bphi=bphi, Bz=bz)

                    if a5 is not None:
                        a5.free(bfield=True)

            # Order by id
            ids  = h5["id"][:]
            time = h5["time"][:]
            idx  = np.lexsort((time, ids))

            return item[idx]


    def get(self, key, ids=None, endcond=None, pncrid=None, SI=True):
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
            SI : bool, optional <br>
                Whether to return quantity in SI units or Ascot units.
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

        if not SI:
            key = alias(key)

            if key in ["energy", "mu"]:
                f      = lambda x: interpret.energy_eV(x)
                val   = np.array([f(x) for x in val]).ravel()

            if key in ["phi", "phimod"]:
                val = val*180/np.pi

            if key in ["mass"]:
                f      = lambda x: interpret.mass_amu(x)
                val   = np.array([f(x) for x in val]).ravel()

            if key in ["charge"]:
                f      = lambda x: interpret.charge_e(x)
                val   = np.array([f(x) for x in val]).ravel()

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
