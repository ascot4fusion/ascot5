"""Import data from Novatron.
"""
import numpy as np
import h5py
import unyt

class ImportNovatron():

    def novatron_bfield(self, filename=None):
        with h5py.File(filename) as h5:
            r = h5["r"][:]
            z = h5["z"][:]
            psi = h5["psi"][:]

        z = np.hstack((-z[::-1], z[1:]))
        psi = np.hstack((np.fliplr(psi[:, 1:]), psi))

        b2ds = {
            "rmin": r[0],
            "rmax": r[-1],
            "nr": r.size,
            "zmin": z[0],
            "zmax": z[-1],
            "nz": z.size,
            "psi": psi,
            "br": psi*0,
            "bphi": psi*0,
            "bz": psi*0,
            "axisr": 0,
            "axisz": 0,
            "psi0": np.amin(psi),
            "psi1": np.amax(psi),
        }

        return "B_2DS", b2ds

    def novatron_efield(self, filename=None):
        with h5py.File(filename) as h5:
            r = h5["r"][:]
            z = h5["z"][:]
            phi = h5["phi"][:]

        z = np.hstack((-z[::-1], z[1:]))
        phi = np.hstack((np.fliplr(phi[:, 1:]), phi))

        e2ds = {
            "rmin": r[0],
            "rmax": r[-1],
            "nr": r.size,
            "zmin": z[0],
            "zmax": z[-1],
            "nz": z.size,
            "vpot": phi,
        }
        return "E_2DS", e2ds

    def novatron_plasma(self, filename=None):
        with h5py.File(filename) as h5:
            r = h5["r"][:]
            z = h5["z"][:]
            #psi = h5["psi_p"][:]