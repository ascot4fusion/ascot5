"""Input representing plasma neutral bakcground.

Neutral data is used in simulations with atomic reactions enabled.
"""
import numpy as np
import h5py

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

class N0_3D(DataGroup):
    """Non-axisymmetric neutral data.

    This input represents neutral density and temperature that can have 3D
    profile. The data is interpolated on an uniform cylindrical grid using
    linear interpolation. The neutral temperature can either be Maxwellian
    or mono-energetic.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    def write_dummy(self, fn):
        return write_hdf5_dummy(fn)

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
               nspecies, anum, znum, density, temperature, maxwellian=1,
               desc=None):
        """Write 3D neutral input in HDF5 file.

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
            Minimum value in phi grid [deg].
        phimax : float
            Maximum value in phi grid [deg].
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

            g.create_dataset("anum",        (nspecies,),  data=anum,
                             dtype="i4")
            g.create_dataset("znum",        (nspecies,),  data=znum,
                             dtype="i4")
            g.create_dataset("maxwellian",  (nspecies,),  data=maxwellian,
                             dtype="i4")

            g.create_dataset("density",     (nspecies,nphi,nz,nr), data=density,
                             dtype="f8")
            g.create_dataset("temperature", (nspecies,nphi,nz,nr), data=temperature,
                             dtype="f8")

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read 3D neutral input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "neutral/N0_3D_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["density"]     = np.transpose(out["density"],     (3,1,2,0))
        out["temperature"] = np.transpose(out["temperature"], (3,1,2,0))

        return out

    @staticmethod
    def write_hdf5_dummy(fn, desc="Dummy"):
        N0Rmin = 0
        N0Rmax = 100
        N0nR   = 4
        N0zmin = -100
        N0zmax = 100
        N0nz   = 3
        N0pmin = 0
        N0pmax = 2*np.pi
        N0np   = 2
        N0spec = 1
        N0anum = np.array([1])
        N0znum = np.array([1])
        N0dens = np.ones( (N0nR,N0np,N0nz,N0spec) )
        N0temp = np.ones( (N0nR,N0np,N0nz,N0spec) )
        write_hdf5(fn,
                   N0Rmin, N0Rmax, N0nR,
                   N0zmin, N0zmax, N0nz,
                   N0pmin, N0pmax, N0np,
                   N0spec, N0anum, N0znum,
                   N0dens, N0temp,
                   desc=desc)
