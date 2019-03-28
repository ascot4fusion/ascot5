"""
State HDF5 IO module.

File: state.py
"""

import numpy as np
import h5py

from a5py.marker.alias import get as alias
import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.plot as plot
from a5py.marker.alias import get as alias
from a5py.marker.endcond import Endcond
from a5py.marker.endcond import endconds

from a5py.ascot5io.ascot5data import AscotData

def read_hdf5(fn, qid, name):
    """
    Read all or specified states.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    read : string list, optional
        Which states are read e.g. "inistate". Default is all.

    Returns
    -------

    Dictionary storing the states that were read.
    """

    path = "results/run_" + qid + "/" + name

    with h5py.File(fn,"r") as f:
        out = {}

        for field in f[path]:
            out[field]           = f[path][field][:]
            out[field + "_unit"] = f[path][field].attrs["unit"]

        # TODO Parse endconditions.

        # Find number of markers and check that no markers share same id
        # (which they shouldn't).
        out["N"] = np.unique(out["id"]).size
        if out["N"] != out["id"].size:
            print("Warning: Markers don't have unique Id.")

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
        Read orbit data to dictionary.
        """
        return read_hdf5(self._file, self.get_qid(), self._path.split("/")[-1])


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
                    item = marker.eval_particle(
                        key, mass=mass, charge=charge,
                        R=h5["r"][:], phi=phi, z=h5["z"][:],
                        vR=h5["vr"][:], vphi=h5["vphi"][:], vz=h5["vz"][:],
                        BR=h5["br"][:], Bphi=h5["bphi"][:], Bz=h5["bz"][:])

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
                Wheter to return quantity in SI units or Ascot units.
        Returns:
            The quantity.
        """
        val = self[key]

        idx = np.ones(val.shape, dtype=bool)

        if endcond is not None and hasattr(self._runnode, "endstate"):
            ec = self._runnode.endstate["endcond"]
            er = self._runnode.endstate["errormsg"]

            endcondlist = [Endcond(ec[i], er[i]) for i in range(ec.size)]
            for i in range(idx.size):
                idx[i] = np.logical_and(idx[i], endcondlist[i] == endcond)

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


    def histogram(self, x=None, y=None, xbins=None, ybins=None, weight=False,
                  logx=False, logy=False, logscale=False, endcond=None,
                  axes=None):
        """
        Make histogram plot.

        """
        if endcond is None and y is None:
            # Repeat 1D histogram for all endconds so we get a stacked histogram
            for ec in endconds.keys():
                self.histogram(x, y=None, xbins=xbins, ybins=None,
                               weight=weight, logx=logx, logy=False,
                               logscale=logscale, endcond=ec, axes=axes)

            return

        ids = self.get("id", endcond=endcond)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, endcond=endcond, SI=False)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, SI=False)

        if logx:
            xc = np.log10(np.absolute(xc))
            xbins = np.linspace(np.log10(xbins[0]), np.log10(xbins[1]),
                                xbins[2])
        else:
            xbins = np.linspace(xbins[0], xbins[1], xbins[2])

        if logy and y is not None:
            yc = np.log10(np.absolute(yc))
            ybins = np.linspace(np.log10(ybins[0]), np.log10(ybins[1]),
                                ybins[2])
        elif y is not None:
            ybins = np.linspace(ybins[0], ybins[1], ybins[2])

        weights=None
        if weight:
            weights = self.get("weight", endcond=endcond)

        plot.plot_histogram(x=xc, y=yc, xbins=xbins, ybins=ybins,
                            weights=weights, logscale=logscale,
                            xlabel=x, ylabel=y, axes=axes)
