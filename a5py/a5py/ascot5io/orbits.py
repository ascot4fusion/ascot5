"""
Orbit data IO module.

File: orbits.py
"""
import numpy as np
import h5py
import scipy.constants as constants

from a5py.ascot5io.ascot5data import AscotData
import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.endcond as endcondmod
import a5py.marker.plot as plot
from a5py.marker.alias import get as alias

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
        g.create_dataset("r",    (N,1), data=data["r"],    dtype="f8")
        g.create_dataset("phi",  (N,1), data=data["phi"],  dtype="f8")
        g.create_dataset("z",    (N,1), data=data["z"],    dtype="f8")
        g.create_dataset("time", (N,1), data=data["time"], dtype="f8")

        g.create_dataset("rho",   (N,1), data=data["rho"],   dtype="f8")
        g.create_dataset("theta", (N,1), data=data["theta"], dtype="f8")
        g.create_dataset("br",    (N,1), data=data["br"],    dtype="f8")
        g.create_dataset("bphi",  (N,1), data=data["bphi"],  dtype="f8")
        g.create_dataset("bz",    (N,1), data=data["bz"],    dtype="f8")

        if "pncrid" in data:
            g.create_dataset("pncrid", (N,1), data=data["pncrid"], dtype="i4")

        if "charge" in data:
            g.create_dataset("charge", (N,1), data=data["charge"], dtype="i4")

        if "weight" in data:
            g.create_dataset("weight", (N,1), data=data["weight"], dtype="f8")

        if "vr" in data:
            g.create_dataset("vr",   (N,1), data=data["vr"],   dtype="f8")
            g.create_dataset("vphi", (N,1), data=data["vphi"], dtype="f8")
            g.create_dataset("vz",   (N,1), data=data["vz"],   dtype="f8")

        if "vpar" in data:
            g.create_dataset("vpar", (N,1), data=data["vpar"], dtype="f8")
            g.create_dataset("mu",   (N,1), data=data["mu"],   dtype="f8")
            g.create_dataset("zeta", (N,1), data=data["zeta"], dtype="f8")


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

            # Is it modulus of poloidal angle.
            if (item is None) and (key == "thetamod"):
                item = (np.mod(h5["theta"][:] + 360, 360) - 180) * np.pi/180

            # See if it is supposed to be read from inistate instead.
            if (item is None) and (key == "mass"):
                item = self._read_from_inistate("mass", h5)
                f    = lambda x: interpret.mass_kg(x)
                mass = np.array([f(x) for x in mass]).ravel()

            if (item is None) and ("mu" in h5keys):
                # HDF5 contains guiding center data
                mass = self._read_from_inistate("mass", h5)

                # Convert guiding-center quantities to SI units
                f      = lambda x: interpret.mass_kg(x)
                mass   = np.array([f(x) for x in mass]).ravel()
                f      = lambda x: interpret.charge_C(x)
                charge = np.array([f(x) for x in h5["charge"][:]]).ravel()
                f      = lambda x: interpret.energy_J(x)
                mu     = np.array([f(x) for x in h5["mu"][:]]).ravel()
                phi    = h5["phi"][:] * np.pi/180

                item = marker.eval_guidingcenter(
                    key, mass=mass, charge=charge,
                    R=h5["r"][:], phi=phi, z=h5["z"][:],
                    mu=mu, vpar=h5["vpar"][:],
                    theta=h5["theta"][:],
                    BR=h5["br"][:], Bphi=h5["bphi"][:],
                    Bz=h5["bz"][:])

            if (item is None) and ("charge" in h5keys):
                # HDF5 contains particle data
                mass = self._read_from_inistate("mass", h5)

                # Convert particle quantities to SI units
                f      = lambda x: interpret.mass_kg(x)
                mass   = np.array([f(x) for x in mass]).ravel()
                f      = lambda x: interpret.charge_C(x)
                charge = np.array([f(x) for x in h5["charge"][:]]).ravel()
                phi    = h5["phi"][:] * np.pi/180

                item = marker.eval_particle(
                    key, mass=mass, charge=charge,
                    R=h5["r"][:], phi=phi, z=h5["z"][:],
                    vR=h5["vr"][:], vphi=h5["vphi"][:], vz=h5["vz"][:],
                    BR=h5["br"][:], Bphi=h5["bphi"][:], Bz=h5["bz"][:])

            if item is None:
                # HDF5 contains field line data

                # Convert particle quantities to SI units
                phi    = h5["phi"][:] * np.pi/180

                # All physical field-line quantities can be get like this.
                item = marker.eval_particle(key, R=h5["r"][:], phi=phi,
                                            z=h5["z"][:], BR=h5["br"][:],
                                            Bphi=h5["bphi"][:],
                                            Bz=h5["bz"][:])

            # Order by id and time
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

        if endcond is not None:
            with self as h5:
                ec = self._read_from_endstate("endcond", h5).ravel()

            idx = np.logical_and( idx, ec == endcondmod.getbin(endcond) )

        if pncrid is not None:
            idx = np.logical_and(idx, self["pncrid"] == pncrid)

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self["id"], ids))

        val = val[idx]

        if not SI:
            key = alias(key)

            if key in ["energy", "mu"]:
                f      = lambda x: interpret.energy_eV(x)
                val   = np.array([f(x) for x in val]).ravel()

            if key in ["phi", "phimod", "theta", "thetamod"]:
                val = val*180/np.pi

            if key in ["mass"]:
                f      = lambda x: interpret.mass_amu(x)
                val   = np.array([f(x) for x in val]).ravel()

            if key in ["charge"]:
                f      = lambda x: interpret.charge_e(x)
                val   = np.array([f(x) for x in val]).ravel()

        return val


    def plot(self, x=None, y=None, z=None, endcond=None, pncrid=None,
             equal=False, log=False, axes=None, ids=None, **kwargs):
        """
        Plot orbits as a continuous line.
        """
        ids0 = ids
        ids = self.get("id", endcond=endcond, ids=ids0)

        xc = np.linspace(0, ids.size, ids.size)
        if x is not None:
            xc = self.get(x, endcond=endcond, pncrid=pncrid, SI=False, ids=ids0)

        yc = None
        if y is not None:
            yc = self.get(y, endcond=endcond, pncrid=pncrid, SI=False, ids=ids0)

        zc = None
        if z is not None:
            zc = self.get(z, endcond=endcond, pncrid=pncrid, SI=False, ids=ids0)

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
                       xlabel=x, ylabel=y, zlabel=z, axes=axes, **kwargs)


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



    def _read_from_inistate(self, key, h5):
        isval = self._runnode.inistate[key]
        isid  = self._runnode.inistate["id"]
        f     = lambda x: isval[isid == x]
        return np.array([f(x) for x in h5["id"][:]])


    def _read_from_endstate(self, key, h5):
        esval = self._runnode.endstate[key]
        esid  = self._runnode.endstate["id"]
        f     = lambda x: esval[esid == x]
        return np.array([f(x) for x in h5["id"][:]])
