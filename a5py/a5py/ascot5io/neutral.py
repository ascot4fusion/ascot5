"""Input representing plasma neutral bakcground.

Neutral data is used in simulations with atomic reactions enabled.
"""
import numpy as np
import h5py

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class N0_1D(DataGroup):
    """Constant-on-flux-surfaces neutral profile.
    """

    def read_hdf5(self):
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

        out["density"]     = np.transpose(out["density"],     (1,0))
        out["temperature"] = np.transpose(out["temperature"], (1,0))

        return out

    @staticmethod
    def write_hdf5(fn, rhomin, rhomax, nrho,
               nspecies, anum, znum, density, temperature, maxwellian=1,
               desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        rhomin : float
            Minimum value in rho grid [1].
        rhomax : float
            Maximum value in rho grid [1].
        nrho : int
            Number of rho grid points.
        nspecies : int
            Number of neutral species.
        anum : array_like (nspecies,1)
            Neutral species' atomic mass number.
        znum array_like (nspecies,1)
            Neutral species' charge number.
        density array_like (nrho,nspecies)
            Neutral species-wise density [m^-3].
        temperature array_like (nrho,nspecies)
            Neutral species-wise temperature [eV].
        maxwellian array_like (nspecies,1)
            Whether species distribution is Maxwellian (1) of monoenergetic (0)
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
        if density.shape != (nrho,nspecies):
            raise ValueError("Density has invalid shape.")
        if temperature.size != (nrho,nspecies):
            raise ValueError("Temperature has invalid shape.")
        if anum.size != nspecies or znum.size != nspecies:
            raise ValueError("Anum or Znum has invalid shape.")
        if maxwellian != 1 and maxwellian.size != nspecies:
            raise ValueError("Failed to interpret maxwellian parameter.")

        parent = "neutral"
        group  = "N0_1D"
        gname  = ""

        # Transpose n0 from (rho, spec) to (spec, rho)
        density     = np.transpose(density)

        if maxwellian == 1:
            maxwellian = np.ones( (int(nspecies),1) )

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("rhomin",   (1,), data=rhomin,   dtype="f8")
            g.create_dataset("rhomax",   (1,), data=rhomax,   dtype="f8")
            g.create_dataset("nrho",     (1,), data=nrho,     dtype="i4")
            g.create_dataset("nspecies", (1,), data=nspecies, dtype="i4")

            g.create_dataset("anum",        (nspecies,),  data=anum,
                             dtype="i4")
            g.create_dataset("znum",        (nspecies,),  data=znum,
                             dtype="i4")
            g.create_dataset("maxwellian",  (nspecies,),  data=maxwellian,
                             dtype="i4")

            g.create_dataset("density",     (nspecies,nrho), data=density,
                             dtype="f8")
            g.create_dataset("temperature", (nspecies,nrho), data=temperature,
                             dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        N0_1D.write_hdf5(
            fn=fn, rhomin=0, rhomax=2, nrho=100, nspecies=1,
            anum=np.array([1]), znum=np.array([1]),
            density=5e16*np.ones( (100, 2) ),
            temperature=1e3*np.ones( (100, 2) ), desc="DUMMY")

class N0_3D(DataGroup):
    """Non-axisymmetric neutral data.

    This input represents neutral density and temperature that can have 3D
    profile. The data is interpolated on an uniform cylindrical grid using
    linear interpolation. The neutral temperature can either be Maxwellian
    or mono-energetic.
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

        out["density"]     = np.transpose(out["density"],     (3,1,2,0))
        out["temperature"] = np.transpose(out["temperature"], (3,1,2,0))

        return out

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
               nspecies, anum, znum, density, temperature, maxwellian=1,
               desc=None):
        """Write input data to the HDF5 file.

        The toroidal angle phi is treated as a periodic coordinate, meaning
        ``A(phi=phimin) == A(phi=phimax)``. However, the phi grid, where input
        arrays are tabulated, is ``linspace(phimin, phimax, nphi+1)[:-1]``
        to avoid storing duplicate data.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        rmin : float
            Minimum value in R grid [m].
        rmax : float
            Maximum value in R grid [m].
        nr : int
            Number of R grid points.
        zmin : float
            Minimum value in z grid [m].
        zmax : float
            Maximum value in z grid [m].
        nz : int
            Number of z grid points.
        phimin : float
            Beginning of the toroidal period [deg].
        phimax : float
            End of the toroidal period [deg].
        nphi : int
            Number of phi grid points.
        nspecies : int
            Number of neutral species.
        anum : array_like (nspecies,1)
            Neutral species' atomic mass number.
        znum : array_like (nspecies,1)
            Neutral species' charge number.
        density : array_like (nr,nphi,nz,nspecies)
            Neutral species-wise density [m^-3].
        temperature : array_like (nr,nphi,nz,nspecies)
            Neutral species-wise temperature [eV].
        maxwellian : array_like (nspecies,1)
            Whether species distribution is Maxwellian (1) of monoenergetic (0)
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
        if density.shape != (nr,nphi,nz,nspecies):
            raise ValueError("Density has invalid shape.")
        if temperature.shape != (nr,nphi,nz,nspecies):
            raise ValueError("Temperature has invalid shape.")
        if anum.size != nspecies or znum.size != nspecies:
            raise ValueError("Anum or Znum has invalid shape.")
        if maxwellian != 1 and maxwellian.size != nspecies:
            raise ValueError("Maxwellian has invalid shape.")

        parent = "neutral"
        group  = "N0_3D"
        gname  = ""

        # Transpose n0 and t0 from (r, phi, z, spec) to (spec, phi, z, r)
        density     = np.transpose(density,(3,1,2,0))
        temperature = np.transpose(temperature,(3,1,2,0))

        if maxwellian == 1:
            maxwellian = np.ones( (int(nspecies),1) )

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("rmin",     (1,), data=rmin,     dtype="f8")
            g.create_dataset("rmax",     (1,), data=rmax,     dtype="f8")
            g.create_dataset("nr",       (1,), data=nr,       dtype="i4")
            g.create_dataset("phimin",   (1,), data=phimin,   dtype="f8")
            g.create_dataset("phimax",   (1,), data=phimax,   dtype="f8")
            g.create_dataset("nphi",     (1,), data=nphi,     dtype="i4")
            g.create_dataset("zmin",     (1,), data=zmin,     dtype="f8")
            g.create_dataset("zmax",     (1,), data=zmax,     dtype="f8")
            g.create_dataset("nz",       (1,), data=nz,       dtype="i4")
            g.create_dataset("nspecies", (1,), data=nspecies, dtype="i4")
            g.create_dataset("anum",        (nspecies,),
                             data=anum, dtype="i4")
            g.create_dataset("znum",        (nspecies,),
                             data=znum, dtype="i4")
            g.create_dataset("maxwellian",  (nspecies,),
                             data=maxwellian, dtype="i4")
            g.create_dataset("density",     (nspecies,nphi,nz,nr),
                             data=density, dtype="f8")
            g.create_dataset("temperature", (nspecies,nphi,nz,nr),
                             data=temperature, dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        N0_3D.write_hdf5(
            fn=fn, rmin=0, rmax=100, nr=3, zmin=-100, zmax=100, nz=3, phimin=0,
            phimax=360, nphi=3, nspecies=1, anum=np.array([1]),
            znum=np.array([1]), density=np.ones((3,3,3,1)),
            temperature=np.ones((3,3,3,1)), desc="DUMMY")
