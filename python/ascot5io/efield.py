"""Electric field input.

Electric field is present in all simulations and it is automatically accounted
for when orbit-following is enabled. To "disable" electric field, set the values
to zero everywhere. The easiest and most computationally effective way to
accomplish this is to use E_TC input.
"""
import h5py
import numpy as np

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

class E_TC(DataGroup):
    """Uniform electric field in Cartesian basis.

    This input fixes the electric field vector, in Cartesian basis, so that the
    field has same value and direction everywhere. This input is meant for
    testing purposes or for disabling electric field in simulations. To disable
    electric field, use :meth:`write_hdf5_dummy`.
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

        return out

    @staticmethod
    def write_hdf5(fn, exyz, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        exyz : array_like (3,1)
            Electric field value in cartesian coordinates [V/m].
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
        if exyz.shape != (3,) and exyz.shape != (3,1):
            raise ValueError("Exyz has wrong shape.")

        parent = "efield"
        group  = "E_TC"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("exyz", (3,1), data=exyz, dtype="f8")

        return gname


    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return E_TC.write_hdf5(fn=fn, exyz=np.array([0,0,0]), desc="DUMMY")

class E_3D(DataGroup):
    """3D electric field that is linearly interpolated.

    This input tabulates the electric-field components on an uniform cylindrical
    grid, and then uses trilinear interpolation to calculate the values during
    the simulation.
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

        out["er"]   = np.transpose(out["er"],   (2,1,0))
        out["ephi"] = np.transpose(out["ephi"], (2,1,0))
        out["ez"]   = np.transpose(out["ez"],   (2,1,0))
        return out

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
                   er, ephi, ez, desc=None):
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
        er : array_like (nr,nphi,nz)
            Electric field R component [V/m].
        ephi : array_like (nr,nphi,nz)
            Electric field phi component [V/m].
        ez : array_like (nr,nphi,nz)
            Electric field z component [V/m].
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
        if er.shape   != (nr,nphi,nz):
            raise ValueError("ER has an inconsinstent shape.")
        if ephi.shape != (nr,nphi,nz):
            raise ValueError("Ephi has an inconsinstent shape.")
        if ez.shape   != (nr,nphi,nz):
            raise ValueError("Ez has an inconsinstent shape.")

        parent = "efield"
        group  = "E_3D"
        gname  = ""

        er   = np.transpose(er,(2,1,0))
        ephi = np.transpose(ephi,(2,1,0))
        ez   = np.transpose(ez,(2,1,0))

        # Create a group for this input.
        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("rmin",          (1,),  data=rmin,   dtype="f8")
            g.create_dataset("rmax",          (1,),  data=rmax,   dtype="f8")
            g.create_dataset("nr",            (1,),  data=nr,     dtype="i4")
            g.create_dataset("phimin",        (1,),  data=phimin, dtype="f8")
            g.create_dataset("phimax",        (1,),  data=phimax, dtype="f8")
            g.create_dataset("nphi",          (1,),  data=nphi,   dtype="i4")
            g.create_dataset("zmin",          (1,),  data=zmin,   dtype="f8")
            g.create_dataset("zmax",          (1,),  data=zmax,   dtype="f8")
            g.create_dataset("nz",            (1,),  data=nz,     dtype="i4")
            g.create_dataset("er",  (nz, nphi, nr),  data=er,     dtype="f8")
            g.create_dataset("ephi",(nz, nphi, nr),  data=ephi,   dtype="f8")
            g.create_dataset("ez",  (nz, nphi, nr),  data=ez,     dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return E_3D.write_hdf5(fn=fn, rmin=1, rmax=10, nr=3, zmin=-10, zmax=10,
                               nz=3, phimin=0, phimax=360, nphi=3,
                               er=np.zeros((3,3,3)), ephi=np.zeros((3,3,3)),
                               ez=np.zeros((3,3,3)), desc="DUMMY")

class E_3DS(DataGroup):
    """3D electric field interpolated with cubic splines.

    This input tabulates the electric-field components on an uniform cylindrical
    grid, and then uses spline interpolation to calculate the values during
    the simulation. Slower and more memory-intensive than :class:`E_3D`, but
    potentially more accurate.
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

        out["er"]   = np.transpose(out["er"],   (2,1,0))
        out["ephi"] = np.transpose(out["ephi"], (2,1,0))
        out["ez"]   = np.transpose(out["ez"],   (2,1,0))
        return out

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
                   er, ephi, ez, desc=None):
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
        er : array_like (nr,nphi,nz)
            Electric field R component [V/m].
        ephi : array_like (nr,nphi,nz)
            Electric field phi component [V/m].
        ez : array_like (nr,nphi,nz)
            Electric field z component [V/m].
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
        if er.shape   != (nr,nphi,nz):
            raise ValueError("ER has an inconsinstent shape.")
        if ephi.shape != (nr,nphi,nz):
            raise ValueError("Ephi has an inconsinstent shape.")
        if ez.shape   != (nr,nphi,nz):
            raise ValueError("Ez has an inconsinstent shape.")

        parent = "efield"
        group  = "E_3DS"
        gname  = ""

        er   = np.transpose(er,(2,1,0))
        ephi = np.transpose(ephi,(2,1,0))
        ez   = np.transpose(ez,(2,1,0))

        # Create a group for this input.
        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("rmin",          (1,),  data=rmin,   dtype="f8")
            g.create_dataset("rmax",          (1,),  data=rmax,   dtype="f8")
            g.create_dataset("nr",            (1,),  data=nr,     dtype="i4")
            g.create_dataset("phimin",        (1,),  data=phimin, dtype="f8")
            g.create_dataset("phimax",        (1,),  data=phimax, dtype="f8")
            g.create_dataset("nphi",          (1,),  data=nphi,   dtype="i4")
            g.create_dataset("zmin",          (1,),  data=zmin,   dtype="f8")
            g.create_dataset("zmax",          (1,),  data=zmax,   dtype="f8")
            g.create_dataset("nz",            (1,),  data=nz,     dtype="i4")
            g.create_dataset("er",  (nz, nphi, nr),  data=er,     dtype="f8")
            g.create_dataset("ephi",(nz, nphi, nr),  data=ephi,   dtype="f8")
            g.create_dataset("ez",  (nz, nphi, nr),  data=ez,     dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return E_3DS.write_hdf5(fn=fn, rmin=1, rmax=10, nr=3, zmin=-10, zmax=10,
                                nz=3, phimin=0, phimax=360, nphi=3,
                                er=np.zeros((3,3,3)), ephi=np.zeros((3,3,3)),
                                ez=np.zeros((3,3,3)), desc="DUMMY")

class E_3DST(DataGroup):
    """Time-dependent 3D electric field interpolated with cubic splines.

    This input tabulates the electric-field components on an uniform
    (R, phi, z, t) grid, and then uses spline interpolation to calculate
    the values during the simulation.
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

        out["er"]   = np.transpose(out["er"],   (3,2,1,0))
        out["ephi"] = np.transpose(out["ephi"], (3,2,1,0))
        out["ez"]   = np.transpose(out["ez"],   (3,2,1,0))
        return out

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
                   tmin, tmax, nt, er, ephi, ez, desc=None):
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
        tmin : float
            Minimum value in t grid [s].
        tmax : float
            Maximum value in t grid [s].
        nt : int
            Number of t grid points.
        er : array_like (nr,nphi,nz,nt)
            Electric field R component [V/m].
        ephi : array_like (nr,nphi,nz,nt)
            Electric field phi component [V/m].
        ez : array_like (nr,nphi,nz,nt)
            Electric field z component [V/m].
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
        if er.shape   != (nr,nphi,nz,nt):
            raise ValueError("ER has an inconsinstent shape.")
        if ephi.shape != (nr,nphi,nz,nt):
            raise ValueError("Ephi has an inconsinstent shape.")
        if ez.shape   != (nr,nphi,nz,nt):
            raise ValueError("Ez has an inconsinstent shape.")

        parent = "efield"
        group  = "E_3DST"
        gname  = ""

        er   = np.transpose(er,(3,2,1,0))
        ephi = np.transpose(ephi,(3,2,1,0))
        ez   = np.transpose(ez,(3,2,1,0))

        # Create a group for this input.
        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("rmin",              (1,), data=rmin,   dtype="f8")
            g.create_dataset("rmax",              (1,), data=rmax,   dtype="f8")
            g.create_dataset("nr",                (1,), data=nr,     dtype="i4")
            g.create_dataset("phimin",            (1,), data=phimin, dtype="f8")
            g.create_dataset("phimax",            (1,), data=phimax, dtype="f8")
            g.create_dataset("nphi",              (1,), data=nphi,   dtype="i4")
            g.create_dataset("zmin",              (1,), data=zmin,   dtype="f8")
            g.create_dataset("zmax",              (1,), data=zmax,   dtype="f8")
            g.create_dataset("nz",                (1,), data=nz,     dtype="i4")
            g.create_dataset("tmin",              (1,), data=tmin,   dtype="f8")
            g.create_dataset("tmax",              (1,), data=tmax,   dtype="f8")
            g.create_dataset("nt",                (1,), data=nt,     dtype="i4")
            g.create_dataset("er",  (nt, nz, nphi, nr), data=er,     dtype="f8")
            g.create_dataset("ephi",(nt, nz, nphi, nr), data=ephi,   dtype="f8")
            g.create_dataset("ez",  (nt, nz, nphi, nr), data=ez,     dtype="f8")

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return E_3DST.write_hdf5(
            fn=fn, rmin=1, rmax=10, nr=3, zmin=-10, zmax=10,
            nz=3, phimin=0, phimax=360, nphi=3, tmin=0, tmax=1,
            nt=3, er=np.zeros((3,3,3,3)),
            ephi=np.zeros((3,3,3,3)), ez=np.zeros((3,3,3,3)),
            desc="DUMMY")

class E_1DS(DataGroup):
    """One-dimensional electric field interpolated with cubic splines.

    This input tabulates the gradient of the electric field potential with
    respect to minor radius on 1D (radial) grid which is then interpolated with
    splines during the simulation.
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
                if key == "nrho":
                    out[key] = int(out[key])

        return out

    @staticmethod
    def write_hdf5(fn, nrho, rhomin, rhomax, dvdrho, reff, desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        nrho : int
            Number of rho slots in data.
        rhomin : float
            Minimum rho value.
        rhomax : float
            Maximum rho value.
        dvdrho : array_like (nrho,1)
            Derivative of electric potential with respect to minor radius [V/m].

            If ``reff = 1 m``, this is essentially equal to ``dv/dr``.
        reff : float
            Effective minor radius of the plasma used to convert ``dv/drho`` to
            SI units as ``drho/dr=1/reff`` [m].
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
        if dvdrho.shape != (nrho,) and dvdrho.shape != (nrho,1):
            raise ValueError("Input dv/drho has a wrong shape.")

        parent = "efield"
        group  = "E_1DS"
        gname  = ""
        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc)
            gname = g.name.split("/")[-1]

            g.create_dataset('nrho',   (1,1),     data=nrho,   dtype='i8')
            g.create_dataset('rhomin', (1,1),     data=rhomin, dtype='f8')
            g.create_dataset('rhomax', (1,1),     data=rhomax, dtype='f8')
            g.create_dataset('dvdrho', (nrho,1),  data=dvdrho, dtype='f8')
            g.create_dataset('reff',   (1,1),     data=reff,   dtype='f8')

        return gname

    @staticmethod
    def write_hdf5_dummy(fn):
        """Write dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        This dummy input sets electric field to zero everywhere.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.

        Returns
        -------
        name : str
            Name, i.e. "<type>_<qid>", of the new input that was written.
        """
        return E_1DS.write_hdf5(fn=fn, nrho=3, rhomin=0, rhomax=1,
                                dvdrho=np.zeros((3,)), reff=1, desc="DUMMY")
