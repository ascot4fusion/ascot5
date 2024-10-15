"""Atomic reaction data HDF5 IO.

The data consists either of cross-sections or rate coefficients.
"""
import h5py
import numpy as np

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class Asigma_loc(DataGroup):
    """Local atomic data.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        fn   = self._root._ascot.file_getpath()
        path = self._path

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]
                if key in ["nreac"]:
                    out[key] = int(out[key])

        return out

    @staticmethod
    def write_hdf5(fn, nreac, z1, a1, z2, a2, reactype, nenergy, energymin,
                   energymax, ndensity, densitymin, densitymax, ntemperature,
                   temperaturemin, temperaturemax, sigma, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Path to hdf5 file.
        nreac : int
            Number of available atomic reactions.
        z1 : array_like (nreac,1)
            Atomic number of test particle.
        a1 : array_like (nreac,1)
            Atomic mass number of test particle.
        z2 : array_like (nreac,1)
            Atomic number of bulk particle.
        a2 : array_like (nreac,1)
            Atomic mass number of bulk particle.
        reactype : array_like (nreac,1)
            Type of atomic reaction.
        nenergy : array_like (nreac,1)
            Number of energy grid points.
        energymin : array_like (nreac,1)
            Energy grid minimum edge [eV].
        energymax : array_like (nreac,1)
            Energy grid maximum edge [eV].
        ndensity : array_like (nreac,1)
            Number of density grid points.
        densitymin : array_like (nreac,1)
            Density grid minimum edge [m^-3].
        densitymax : array_like (nreac,1)
            Density grid maximum edge [m^-3].
        ntemperature : array_like (nreac,1)
            Number of temperature grid points.
        temperaturemin : array_like (nreac,1)
            Temperature grid minimum edge [eV].
        temperaturemax : array_like (nreac,1)
            Temperature grid maximum edge [eV].
        sigma : array_like (1,sum(nenergy[i]*ndensity[i]*ntemperature[i]))
            Reaction cross-section or other probability data [cm^2 or other].
        desc : str, optional
            Input description.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.

        Raises
        ------
        ValueError
            If inputs were not consistent.
        """
        if z1.size != nreac:
            raise ValueError("Invalid number of reactions.")
        n = np.zeros(nreac, dtype=int)
        ntot = 0
        for i in range(0, nreac):
            n[i] = nenergy[i] * ndensity[i] * ntemperature[i]
            ntot += n[i]
        if sigma.shape != (1,ntot):
            raise ValueError("Invalid size for sigma.")
        parent = "asigma"
        group  = "asigma_loc"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nreac',          (1,1), data=nreac,      dtype='i4')
            g.create_dataset('z1',         (nreac,1), data=z1,         dtype='i4')
            g.create_dataset('a1',         (nreac,1), data=a1,         dtype='i4')
            g.create_dataset('z2',         (nreac,1), data=z2,         dtype='i4')
            g.create_dataset('a2',         (nreac,1), data=a2,         dtype='i4')
            g.create_dataset('reactype',   (nreac,1), data=reactype,   dtype='i4')
            g.create_dataset('nenergy',    (nreac,1), data=nenergy,    dtype='i4')
            g.create_dataset('energymin',  (nreac,1), data=energymin,  dtype='f8')
            g.create_dataset('energymax',  (nreac,1), data=energymax,  dtype='f8')
            g.create_dataset('ndensity',   (nreac,1), data=ndensity,   dtype='i4')
            g.create_dataset('densitymin', (nreac,1), data=densitymin, dtype='f8')
            g.create_dataset('densitymax', (nreac,1), data=densitymax, dtype='f8')
            g.create_dataset('ntemperature',   (nreac,1), data=ntemperature,
                             dtype='i4')
            g.create_dataset('temperaturemin', (nreac,1), data=temperaturemin,
                             dtype='f8')
            g.create_dataset('temperaturemax', (nreac,1), data=temperaturemax,
                             dtype='f8')
            g.create_dataset('sigma',      (1,ntot),  data=sigma,      dtype='f8')

        return gname

    @staticmethod
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        N_reac    = 1
        z_1       = 1 + np.zeros(N_reac, dtype=int)
        a_1       = 1 + np.zeros(N_reac, dtype=int)
        z_2       = 1 + np.zeros(N_reac, dtype=int)
        a_2       = 1 + np.zeros(N_reac, dtype=int)
        reac_type = 7 + np.zeros(N_reac, dtype=int)
        N_E       = 3    + np.zeros(N_reac, dtype=int)
        E_min     = 1e3  + np.zeros(N_reac, dtype=float)
        E_max     = 1e4  + np.zeros(N_reac, dtype=float)
        N_n       = 4    + np.zeros(N_reac, dtype=int)
        n_min     = 1e18 + np.zeros(N_reac, dtype=float)
        n_max     = 1e20 + np.zeros(N_reac, dtype=float)
        N_T       = 5    + np.zeros(N_reac, dtype=int)
        T_min     = 1e3  + np.zeros(N_reac, dtype=float)
        T_max     = 1e4  + np.zeros(N_reac, dtype=float)
        sigma = np.zeros((1,3*4*5))
        return {"nreac":N_reac, "z1":z_1, "a1":a_1, "z2":z_2, "a2":a_2,
                "reactype":reac_type, "nenergy":N_E, "energymin":E_min,
                "energymax":E_max, "ndensity":N_n, "densitymin":n_min,
                "densitymax":n_max, "ntemperature":N_T, "temperaturemin":T_min,
                "temperaturemax":T_max, "sigma":sigma}
