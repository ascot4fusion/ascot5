"""
Orbits HDF5 IO module.

File: orbits.py
"""
import numpy as np
import h5py
import scipy.constants as constants

from a5py.ascot5io.ascot5data import AscotData
import a5py.marker.interpret as interpret
import a5py.marker as marker
import a5py.marker.plot as plot

def read_hdf5(fn, qid):
    """
    Read orbits.

    Parameters
    ----------

    fn : str
        Full path to HDF5 file.
    qid : str
        qid of the run these orbit data correspond to.

    Returns
    -------

    Dictionary storing the orbits that were read.
    """

    with h5py.File(fn,"r") as f:
        orbits = f["/results/run-"+qid+"/orbits"]

        out = {}

        # Read data from file.
        for field in orbits:
            out[field]           = orbits[field][:]
            out[field + "_unit"] = orbits[field].attrs["unit"]

        # Find how many markers we have and their ids.
        out["N"]        = (np.unique(orbits["id"][:])).size
        out["uniqueId"] = np.unique(orbits["id"][:])

        # Sort fields by id (major) and time (minor), both ascending.
        if out["N"] > 0:
            ind = np.lexsort((out["time"], out["id"]))
            for field in orbits:
                if field[-4:] != "unit":
                    out[field] = out[field][ind]

    return out

class Orbits(AscotData):

    def __init__(self, hdf5, runnode):
        self._runnode = runnode
        super().__init__(hdf5)

    def read(self):
        return read_hdf5(self._file, self.get_qid())

    def __getitem__(self, key):
        with self as h5:

            unit = None
            if isinstance(key, tuple):
                unit = key[1]
                key  = key[0]

            # Remove underscores, spaces and uppercase from key and hdf5 keys.
            def cleankey(k):
                k = k.lower().replace("_", "")
                return k.replace(" ", "")

            key    = cleankey(key)
            h5keys = [cleankey(x) for x in h5.keys()]
            if unit is not None:
                unit   = cleankey(unit)

            # See if the field can be read directly and without conversions
            if key in h5keys and key not in ["m", "mass", "q", "charge" ,"phi"]:
                dirtykey = list(h5.keys())[h5keys.index(key)]
                item = h5[dirtykey][:]

            elif "mu" in h5keys:
                # This contains guiding center data

                # Take mass from inistate
                inistatemass = self._runnode.inistate["mass"]
                inistateid   = self._runnode.inistate["id"]
                f    = lambda x: inistatemass[inistateid == x]
                mass = np.array([f(x) for x in h5["id"][:]])

                # Convert to SI units
                f      = lambda x: interpret.mass_kg(x)
                mass   = np.array([f(x) for x in mass]).ravel()
                f      = lambda x: interpret.charge_C(x)
                charge = np.array([f(x) for x in h5["charge"][:]]).ravel()
                phi    = h5["phi"][:] * np.pi/180

                if key in ["m", "mass"]:
                    item = mass
                elif key in ["q", "charge"]:
                    item = charge
                elif key is "phi":
                    item = phi
                else:
                    item = marker.eval_guidingcenter(
                        key, mass=mass, charge=charge,
                        R=h5["R"][:], phi=h5["phi"][:], z=h5["z"][:],
                        mu=h5["mu"][:], vpar=h5["vpar"][:],
                        theta=h5["theta"][:],
                        BR=h5["B_R"][:], Bphi=h5["B_phi"][:],
                        Bz=h5["B_z"][:])

            elif "charge" in h5keys:
                # This contains particle data

                # Take mass from inistate
                inistatemass = self._runnode.inistate["mass"]
                inistateid   = self._runnode.inistate["id"]
                f    = lambda x: inistatemass[inistateid == x]
                mass = np.array([f(x) for x in h5["id"][:]])

                # Convert to SI units
                f      = lambda x: interpret.mass_kg(x)
                mass   = np.array([f(x) for x in mass]).ravel()
                f      = lambda x: interpret.charge_C(x)
                charge = np.array([f(x) for x in h5["charge"][:]]).ravel()
                phi    = h5["phi"][:] * np.pi/180

                if key in ["m", "mass"]:
                    item = mass
                elif key in ["q", "charge"]:
                    item = charge
                elif key is "phi":
                    item = phi
                else:
                    item = marker.eval_particle(
                        key, mass=mass, charge=charge,
                        R=h5["R"][:], phi=h5["phi"][:], z=h5["z"][:],
                        vR=h5["v_R"][:], vphi=h5["v_phi"][:], vz=h5["v_z"][:],
                        BR=h5["B_R"][:], Bphi=h5["B_phi"][:], Bz=h5["B_z"][:])

            else:
                # This contains field line data

                if key is "phi":
                    item = h5["phi"][:] * np.pi/180
                else:
                    # We only need Bnorm which we can get e.g. from this method
                    item = marker.eval_particle(key, BR=h5["B_R"][:],
                                                Bphi=h5["B_phi"][:],
                                                Bz=h5["B_z"][:])

            # Unit conversions
            if unit is not None:
                unit = unit.lower()
                if unit == "ev":
                    f    = lambda x: interpret.energy_eV(x)
                    item = np.array([f(x) for x in item])
                elif unit == "e":
                    f    = lambda x: interpret.charge_e(x)
                    item = np.array([f(x) for x in item])
                elif unit == "amu":
                    f    = lambda x: interpret.mass_amu(x)
                    item = np.array([f(x) for x in item])
                elif unit == "deg":
                    f    = lambda x: x * 180 / np.pi
                    item = np.array([f(x) for x in item])

            # Order by id and time
            ids  = h5["id"][:]
            time = h5["time"][:]
            idx  = np.lexsort((time, ids))

            return item[idx]

    def plot_orbit(self, x, y, equal=False):
        x = self[x]
        y = self[y]
        ids = self["id"]

        plot.plot_orbit(x, y, ids, equal=equal)
