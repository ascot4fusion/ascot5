"""Magnetic field input.

Present in every simulation.
"""
import h5py
import numpy as np

from ._iohelpers.fileapi import add_group
from ._iohelpers.treedata import DataGroup

import a5py.preprocessing.analyticequilibrium as psifun

class B_TC(DataGroup):
    """
    Object representing B_TC data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())


    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    @staticmethod
    def write_hdf5(fn, bxyz, jacobian, rhoval, psival=None, axisr=1, axisz=0,
                   desc=None):
        """Write trivial cartesian magnetic field input to HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        bxyz : array_like (3,1)
            Magnetic field in cartesian coordinates at origo.
        jacobian : array_like (3,3)
            Magnetic field Jacobian, jacobian[i,j] = dB_i/dx_j
        rhoval: float
            Constant rho value.
        psival: float, optional
            Constant psi value. If None, same as rhoval.
        axisr: float, optional
            Magnetic axis R coordinate.
        axisz: real, optional
            Magnetic axis z coordinate.
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
        if bxyz.shape != (3,1): raise ValueError("Invalid shape for Bxyz.")

        parent = "bfield"
        group  = "B_TC"
        gname  = ""

        if psival is None:
            psival = rhoval

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("bxyz",     (3,1), data=bxyz,     dtype="f8")
            g.create_dataset("jacobian", (3,3), data=jacobian, dtype="f8")
            g.create_dataset("rhoval",   (1,),  data=rhoval,   dtype="f8")
            g.create_dataset("psival",   (1,),  data=psival,   dtype="f8")
            g.create_dataset("axisr",    (1,),  data=axisr,    dtype="f8")
            g.create_dataset("axisz",    (1,),  data=axisz,    dtype="f8")

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read Cartesian magnetic field input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "bfield/B_TC_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        return out

    @staticmethod
    def write_hdf5_dummy(fn, desc="Dummy"):
        """
        Write dummy data.
        """
        B   = np.array([1,0,3])
        jac = np.array([ [0,0,0.01], [0,0,0], [0,0,0] ])
        return write_hdf5(fn, B, jacobian=jac,
                          rhoval=0.5, psival=1.5, axisr=6, axisz=0, desc=desc)

class B_GS(DataGroup):
    """
    Object representing B_GS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    @staticmethod
    def write_hdf5(fn, r0, z0, bphi0, psimult, coefficients, psi0=None, psi1=0,
                   raxis=None, zaxis=None, nripple=0, a0=2, alpha0=2,
                   delta0=0.05, desc=None):
        """Write analytical tokamak magnetic field input in HDF5 file.

        If psi0 is None, magnetic axis location is found numerically and psi0
        value is evaluated at that point. The analytical field is defined by
        assuming midplane is at z=0, but it is moved to z=z0 here.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        r0 : float
            Major radius R coordinate [m].
        z0 : float
            Distance by which midplane is moved from z=0 plane [m].
        bphi0 : float
            Toroidal field at axis [T].
        psimult : float
            Scaling factor for psi.
        coefficients : array_like (13,1)
            Coefficients defining psi [c0, c1, ..., c11, A].
        psi0 : float, optional
            Poloidal flux at magnetic axis.
        psi1 : float, optional
            Poloidal flux at the separatrix.
        raxis : float, optional
            Magnetic axis R coordinate [m].
        zaxis : float, optional
            Magnetic axis z coordinate [m].
        nripple : float, optional
            Number of TF coils.
        a0 : float, optional
            Minor radius. [m]
        alpha0 : float, optional
            Ripple penetration.
        delta0 : float, optional
            Ripple (maximum) strength.
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

        if coefficients.size != 13:
            raise ValueError("Coefficients invalid.")

        parent = "bfield"
        group  = "B_GS"
        gname  = ""

        c = coefficients # For shorter notation.

        # Search for magnetic axis psi
        if psi0 is None:
            x = psifun.find_axis(r0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                                 c[7], c[8], c[9], c[10], c[11], c[12])
            psi0 = psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                               c[5], c[6], c[7], c[8], c[9], c[10], c[11],
                               c[12]) * psimult * 1.001
            raxis = x[0]*r0
            zaxis = x[1]*r0

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("r0",           (1,),   data=r0,           dtype='f8')
            g.create_dataset("z0",           (1,),   data=z0,           dtype='f8')
            g.create_dataset("raxis",        (1,),   data=raxis,        dtype='f8')
            g.create_dataset("zaxis",        (1,),   data=zaxis,        dtype='f8')
            g.create_dataset("bphi0",        (1,),   data=bphi0,        dtype='f8')
            g.create_dataset("psi0",         (1,),   data=psi0,         dtype='f8')
            g.create_dataset("psi1",         (1,),   data=psi1,         dtype='f8')
            g.create_dataset("psimult",      (1,),   data=psimult,      dtype='f8')
            g.create_dataset("coefficients", (13,1), data=coefficients, dtype='f8')
            g.create_dataset("nripple",      (1,),   data=nripple,      dtype='i8')
            g.create_dataset("a0",           (1,),   data=a0,           dtype='f8')
            g.create_dataset("alpha0",       (1,),   data=alpha0,       dtype='f8')
            g.create_dataset("delta0",       (1,),   data=delta0,       dtype='f8')

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read analytic magnetic field input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "bfield/B_GS_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        return out

    @staticmethod
    def write_hdf5_dummy(fn, kind="GS", desc="Dummy"):
        r0      = 6.2
        z0      = 0
        bphi0   = 5.3
        psimult = 200

        # ITER-like but circular equilibrium
        coefficients = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                                  6.200e-03, -1.205e-03, -3.701e-05,  0,
                                  0,          0,          0,          0,    -0.155])

        if kind == "GS":
            return write_hdf5(fn, r0, z0, bphi0, psimult, coefficients,
                              nripple=1, a0=2, alpha0=2, delta0=0.05, desc=desc)

        if kind == "2DS":
            R = (1,  6, 50)
            z = (-4, 4, 100)
            return write_hdf5_B_2DS(fn, r0, z0, bphi0, psimult, coefficients,
                                    R[0], R[1], R[2], z[0], z[1], z[2], psi0=None,
                                    desc=desc)

        if kind == "3DS":
            R   = (1,  6, 50)
            z   = (-4, 4, 100)
            phi = (0, 360, 100)
            return write_hdf5_B_3DS(fn, r0, z0, bphi0, psimult, coefficients,
                                    1, 2, 2, 0.05,
                                    R[0], R[1], R[2], z[0], z[1], z[2],
                                    phi[0], phi[1], phi[2], psi0=None, desc=desc)

        @staticmethod
        def write_hdf5_B_2DS(fn, r0, z0, bphi0, psimult, coefficients,
                             rmin, rmax, nr, zmin, zmax, nz, psi0=None, desc=None):
            """
            Write analytical tokamak magnetic field as a 2D field input in HDF5 file.

            Args:
            fn : str <br>
                Full path to the HDF5 file.
            r0 : float <br>
                Major axis R coordinate [m].
            z0 : float <br>
                Major axis z coordinate [m].
            bphi0 : float <br>
                Toroidal field at axis.
            psimult : real
                Scaling factor for psi.
            coefficients : real 13 x 1 numpy array
                Coefficients defining psi.
            rmin, rmax, zmin, zmax : real
                Edges of the uniform Rz-grid.
            nr, nz : int
                Number of Rz-grid points.
            psi0 : float, optional <br>
                Poloidal flux at magnetic axis.
            desc : str, optional <br>
                Input description.
            """

            rgrid = np.linspace(rmin, rmax, nr)
            zgrid = np.linspace(zmin, zmax, nz)

            zg, rg = np.meshgrid(zgrid, rgrid);

            c = coefficients # For shorter notation.
            psirz = psimult * psifun.psi0(
                rg/r0, zg/r0, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8],
                c[9], c[10], c[11], c[12])

            br = np.zeros((nr,nz))
            bz = np.zeros((nr,nz))
            bphi = ( r0 / rg ) * bphi0

            # search for magnetic axis if not given
            if psi0 == None:
                x = psifun.find_axis(r0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                                     c[7], c[8], c[9], c[10], c[11], c[12])
                psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                                            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
                r0 = x[0]*r0
                z0 = x[1]*r0

                psi1 = 0

            return {"rmin" : rmin, "rmax" : rmax, "nr" : nr, "zmin" : zmin,
                    "zmax" : zmax, "nz" : nz, "r0" : r0, "z0" : z0, "psirz" : psirz,
                    "psi0" : psi0, "psi1" : psi1, "br" : br, "bphi" : bphi, "bz" : bz}

        @staticmethod
        def write_hdf5_B_3DS(fn, r0, z0, bphi0, psimult, coefficients,
                             nripple, a0, alpha0, delta0,
                             rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
                             psi0=None, raxis=None, zaxis=None, desc=None):
            """
            Write analytical tokamak magnetic field as a 3D field input in HDF5 file.

            Parameters
            ----------
            fn : str
                Full path to the HDF5 file.
            r0, z0 : real
                Major axis Rz-coordinates.
            bphi0 : real
                Toroidal field at axis.
            psimult : real
                Scaling factor for psi.
            coefficients : real 13 x 1 numpy array
                Coefficients defining psi.
            nripple : real
                Number of TF coils.
            a0 : real
                Minor radius [m].
            alpha0 : real
                Ripple penetration.
            delta0 : real
                Ripple strength.
            rmin, rmax, zmin, zmax, phimin, phimax : real
                Edges of the uniform Rphiz-grid.
            nr, nz, nphi : int
                Number of Rphiz-grid points.
            psi0 : float, optional
                Poloidal flux at magnetic axis.
            raxis : float, optional
                Magnetic axis R coordinate [m].
            zaxis : float, optional
                Magnetic axis z coordinate [m].
            desc : str, optional
                Input description.
            """

            rgrid = np.linspace(rmin, rmax, nr)
            zgrid = np.linspace(zmin, zmax, nz)

            zg, rg = np.meshgrid(zgrid, rgrid);

            c = coefficients # For shorter notation.
            psirz = psimult*psifun.psi0(rg/r0,zg/r0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                        c[7],c[8],c[9],c[10],c[11],c[12])

            br = np.zeros((nr, nphi, nz))
            bz = np.zeros((nr, nphi, nz))
            bphi = np.zeros((nr, nphi, nz))

            # Ripple
            radius = np.sqrt( ( rg - r0 ) * ( rg - r0 ) + ( zg - z0 ) * ( zg - z0 ))
            theta = np.arctan2( zg - z0, rg - r0 )
            delta = delta0 * np.exp(-0.5*theta*theta) * np.power( radius / a0, alpha0 )

            phigrid = np.linspace(phimin,phimax,nphi+1)
            phigrid = phigrid[0:-1]

            for i in range(0,nphi):
                bphi[:,i,:] = ((r0/rg)*bphi0 * ( 1 + delta * np.cos(Nripple * phigrid[i]*2*np.pi/360) ))

            # search for magnetic axis if not given
            if psi0 == None:
                x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                                     c[7], c[8], c[9], c[10], c[11], c[12])
                psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                                            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
                raxis = x[0]*R0
                zaxis = x[1]*R0

            psi1 = 0

            return {"rmin" : rmin, "rmax" : rmax, "nr" : nr, "zmin" : zmin,
                    "zmax" : zmax, "nz" : nz, "phimin" : phimin, "phimax" : phimax,
                    "nphi" : nphi, "raxis" : raxis, "zaxis" : zaxis, "psirz" : psirz,
                    "psi0" : psi0, "psi1" : psi1, "br" : br, "bphi" : bphi, "bz" : bz}

class B_2DS(DataGroup):
    """
    Object representing B_2DS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz,
                   axisr, axisz, psi, psi0, psi1,
                   br, bphi, bz, desc=None):
        """Write 2DS magnetic field input in HDF5 file.

        Note that br and bz should not include the equilibrium component of the
        magnetic field as that is calculated from psi by ASCOT5 during the
        simulation.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        rmin : float
            R grid min edge [m].
        rmax : float
            R grid max edge [m].
        nr : int
            Number of R grid points.
        zmin : float
            z grid min edge [m].
        zmax : float
            z grid max edge [m].
        nz : int
            Number of z grid points.
        axisr : float
            Magnetic axis R coordinate [m].
        axisz : float
            Magnetic axis z coordinate [m].
        psi0 : float
            On-axis poloidal flux value [Vs/m].
        psi1 : float
            Separatrix poloidal flux value [Vs/m].
        psi : array_like (nr, nz)
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (nr,nz)
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
        bphi : array_like (nr,nz)
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nr,nz)
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
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

        if psi.shape  != (nr,nz): raise ValueError("Inconsistent shape for psi.")
        if br.shape   != (nr,nz): raise ValueError("Inconsistent shape for br.")
        if bphi.shape != (nr,nz): raise ValueError("Inconsistent shape for bphi.")
        if bz.shape   != (nr,nz): raise ValueError("Inconsistent shape for bz.")

        psi  = np.transpose(psi)
        br   = np.transpose(br)
        bphi = np.transpose(bphi)
        bz   = np.transpose(bz)

        parent = "bfield"
        group  = "B_2DS"
        gname  = ""

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("rmin",  (1,), data=rmin,  dtype="f8")
            g.create_dataset("rmax",  (1,), data=rmax,  dtype="f8")
            g.create_dataset("nr",    (1,), data=nr,    dtype="i4")
            g.create_dataset("zmin",  (1,), data=zmin,  dtype="f8")
            g.create_dataset("zmax",  (1,), data=zmax,  dtype="f8")
            g.create_dataset("nz",    (1,), data=nz,    dtype="i4")
            g.create_dataset("axisr", (1,), data=axisr, dtype="f8")
            g.create_dataset("axisz", (1,), data=axisz, dtype="f8")
            g.create_dataset("psi0",  (1,), data=psi0,  dtype="f8")
            g.create_dataset("psi1",  (1,), data=psi1,  dtype="f8")

            g.create_dataset("psi",  (nz, nr), data=psi,  dtype="f8")
            g.create_dataset("br",   (nz, nr), data=br,   dtype="f8")
            g.create_dataset("bphi", (nz, nr), data=bphi, dtype="f8")
            g.create_dataset("bz",   (nz, nr), data=bz,   dtype="f8")

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read 2D magnetic field input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "bfield/B_2DS_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["psi"]  = np.transpose(out["psi"])
        out["br"]   = np.transpose(out["br"])
        out["bphi"] = np.transpose(out["bphi"])
        out["bz"]   = np.transpose(out["bz"])
        return out

    @staticmethod
    def write_hdf5_dummy(fn, desc="Dummy"):
        import a5py.ascot5io.B_GS as B_GS
        return B_GS.write_hdf5_dummy(fn, kind="2DS", desc=desc)

class B_3DS(DataGroup):
    """
    Object representing B_3DS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    @staticmethod
    def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                   b_phimin, b_phimax, b_nphi,
                   axisr, axisz, psi, psi0, psi1, br, bphi, bz,
                   psi_rmin=None, psi_rmax=None, psi_nr=None,
                   psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
        """Write 3DS magnetic field input in HDF5 file.

        Note that br and bz should not include the equilibrium component of the
        magnetic field as that is calculated from psi by ASCOT5 during the
        simulation.

        It is possible to use different Rz grids for psi and magnetic field
        components by giving Rz grid for psi separately.

        The toroidal angle phi is treated as a periodic coordinate meaning that
        B(phi) = B(phi + n*(b_phimax - b_phimin)). Do note that to avoid
        duplicate data, the last points in phi axis in B data are not at
        b_phimax, i.e. br[:,-1,:] != BR(phi=b_phimax).

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        b_rmin : float
            Magnetic field data R grid min edge [m].
        b_rmax : float
            Magnetic field data R grid max edge [m].
        b_nr : int
            Number of R grid points in magnetic field data.
        b_zmin : float
            Magnetic field data z grid min edge [m].
        b_zmax : float
            Magnetic field data z grid max edge [m].
        b_nz : int
            Number of z grid points in magnetic field data.
        b_phimin : float
            Magnetic field data phi grid min edge [deg].
        b_phimax : float
            Magnetic field data phi grid max edge [deg].
        b_nphi : int
            Number of phi grid points in magnetic field data.
        axisr : float
            Magnetic axis R coordinate [m].
        axisz : float
            Magnetic axis z coordinate [m].
        psi0 : float
            On-axis poloidal flux value [Vs/m].
        psi1 : float
            Separatrix poloidal flux value [Vs/m].
        psi : array_like (nr, nz)
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (nr,nphi,nz)
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
        bphi : array_like (nr,nphi,nz)
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nr,nphi,nz)
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
        psi_rmin : float, optional
            Psi data R grid min edge [m].
        psi_rmax : float, optional
            Psi data R grid max edge [m].
        psi_nr : int, optional
            Number of R grid points in psi data.
        psi_zmin : float, optional
            Psi data z grid min edge [m].
        psi_zmax : float, optional
            Psi data z grid max edge [m].
        psi_nz : int, optional
            Number of z grid points in psi data.
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
        parent = "bfield"
        group  = "B_3DS"
        gname  = ""

        # Define psigrid to be same as Bgrid if not stated otherwise.
        if(psi_rmin is None or psi_rmax is None or psi_nr is None or
           psi_zmin is None or psi_zmax is None or psi_nz is None):
            psi_rmin = b_rmin
            psi_rmax = b_rmax
            psi_nr   = b_nr
            psi_zmin = b_zmin
            psi_zmax = b_zmax
            psi_nz   = b_nz

        if psi.shape  != (psi_nr,psi_nz):
            raise ValueError("Inconsistent shape foor psi.")
        if br.shape   != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape foor br.")
        if bphi.shape != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape foor bphi.")
        if bz.shape   != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape foor bz.")

        psi  = np.transpose(psi)
        br   = np.transpose(br,   (2,1,0))
        bphi = np.transpose(bphi, (2,1,0))
        bz   = np.transpose(bz,   (2,1,0))

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("b_rmin",   (1,), data=b_rmin,   dtype="f8")
            g.create_dataset("b_rmax",   (1,), data=b_rmax,   dtype="f8")
            g.create_dataset("b_nr",     (1,), data=b_nr,     dtype="i4")
            g.create_dataset("b_phimin", (1,), data=b_phimin, dtype="f8")
            g.create_dataset("b_phimax", (1,), data=b_phimax, dtype="f8")
            g.create_dataset("b_nphi",   (1,), data=b_nphi,   dtype="i4")
            g.create_dataset("b_zmin",   (1,), data=b_zmin,   dtype="f8")
            g.create_dataset("b_zmax",   (1,), data=b_zmax,   dtype="f8")
            g.create_dataset("b_nz",     (1,), data=b_nz,     dtype="i4")
            g.create_dataset("psi_rmin", (1,), data=psi_rmin, dtype="f8")
            g.create_dataset("psi_rmax", (1,), data=psi_rmax, dtype="f8")
            g.create_dataset("psi_nr",   (1,), data=psi_nr,   dtype="i4")
            g.create_dataset("psi_zmin", (1,), data=psi_zmin, dtype="f8")
            g.create_dataset("psi_zmax", (1,), data=psi_zmax, dtype="f8")
            g.create_dataset("psi_nz",   (1,), data=psi_nz,   dtype="i4")
            g.create_dataset("axisr",    (1,), data=axisr,    dtype="f8")
            g.create_dataset("axisz",    (1,), data=axisz,    dtype="f8")
            g.create_dataset("psi0",     (1,), data=psi0,     dtype="f8")
            g.create_dataset("psi1",     (1,), data=psi1,     dtype="f8")

            g.create_dataset("psi",  (psi_nz, psi_nr),     data=psi,  dtype="f8")
            g.create_dataset("br",   (b_nz, b_nphi, b_nr), data=br,   dtype="f8")
            g.create_dataset("bphi", (b_nz, b_nphi, b_nr), data=bphi, dtype="f8")
            g.create_dataset("bz",   (b_nz, b_nphi, b_nr), data=bz,   dtype="f8")

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read 3D magnetic field input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "bfield/B_3DS_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["psi"]  = np.transpose(out["psi"])
        out["br"]   = np.transpose(out["br"],   (2,1,0))
        out["bphi"] = np.transpose(out["bphi"], (2,1,0))
        out["bz"]   = np.transpose(out["bz"],   (2,1,0))
        return out

    @staticmethod
    def write_hdf5_dummy(fn, desc="Dummy"):
        import a5py.ascot5io.B_GS as B_GS
        return B_GS.write_hdf5_dummy(fn, kind="3DS", desc=desc)

    @staticmethod
    def write_hdf5_B_2DS(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                         b_phimin, b_phimax, b_nphi,
                         axisr, axisz, psi, psi0, psi1, br, bphi, bz,
                         desc=None):
        """
        Convert the 3D data to 2D.
        """
        bphi = bphi=np.mean(bphi,axis=1)
        B_2DS.write_hdf5(fn, rmin=b_rmin, rmax=b_rmax,
                         zmin=b_zmin, zmax=b_zmax, nr=b_nr,
                         nz=b_nz, axisr=axisr, axisz=axisz,
                         psi0=psi0, psi1=psi1, psi=psi,
                         br=bphi*0, bphi=bphi, bz=bphi*0)

class B_3DST(DataGroup):
    """
    Object representing B_3DS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    @staticmethod
    def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                   b_phimin, b_phimax, b_nphi, b_tmin, b_tmax, b_nt,
                   axisr, axisz, psi, psi0, psi1, br, bphi, bz,
                   psi_rmin=None, psi_rmax=None, psi_nr=None,
                   psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
        """
        Write 3DST magnetic field input in HDF5 file.

        Note that br and bz should not include the equilibrium component of the
        magnetic field as that is calculated from psi by ASCOT5 during the
        simulation.

        It is possible to use different Rz grids for psi and magnetic field
        components by giving Rz grid for psi separately.

        The toroidal angle phi is treated as a periodic coordinate meaning that
        B(phi) = B(phi + n*(b_phimax - b_phimin)). Do note that to avoid
        duplicate data, the last points in phi axis in B data are not at
        b_phimax, i.e. br[:,-1,:,:] != BR(phi=b_phimax).

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        b_rmin : float
            Magnetic field data R grid min edge [m].
        b_rmax : float
            Magnetic field data R grid max edge [m].
        b_nr : int
            Number of R grid points in magnetic field data.
        b_zmin : float
            Magnetic field data z grid min edge [m].
        b_zmax : float
            Magnetic field data z grid max edge [m].
        b_nz : int
            Number of z grid points in magnetic field data.
        b_phimin : float
            Magnetic field data phi grid min edge [deg].
        b_phimax : float
            Magnetic field data phi grid max edge [deg].
        b_nphi : int
            Number of phi grid points in magnetic field data.
        b_tmin : float
            Magnetic field data time grid min edge [s].
        b_tmax : float
            Magnetic field data time grid max edge [s].
        b_nt : int
            Number of t grid points in magnetic field data.
        axisr : float
            Magnetic axis R coordinate [m].
        axisz : float
            Magnetic axis z coordinate [m].
        psi0 : float
            On-axis poloidal flux value [Vs/m].
        psi1 : float
            Separatrix poloidal flux value [Vs/m].
        psi : array_like (nr, nz)
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (nr,nphi,nz,nt)
            Magnetic field R component (excl. equilibrium comp.).
        bphi : array_like (nr,nphi,nz,nt)
            Magnetic field phi component on Rz grid [T].
        bz : array_like (nr,nphi,nz,nt)
            Magnetic field z component (excl. equilibrium comp.).
        psi_rmin : float, optional
            Psi data R grid min edge [m].
        psi_rmax : float, optional
            Psi data R grid max edge [m].
        psi_nr : int, optional
            Number of R grid points in psi data.
        psi_zmin : float, optional
            Psi data z grid min edge [m].
        psi_zmax : float, optional
            Psi data z grid max edge [m].
        psi_nz : int, optional
            Number of z grid points in psi data.
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

        parent = "bfield"
        group  = "B_3DST"
        gname  = ""

        # Define psigrid to be same as Bgrid if not stated otherwise.
        if(psi_rmin is None or psi_rmax is None or psi_nr is None or
           psi_zmin is None or psi_zmax is None or psi_nz is None):
            psi_rmin = b_rmin
            psi_rmax = b_rmax
            psi_nr   = b_nr
            psi_zmin = b_zmin
            psi_zmax = b_zmax
            psi_nz   = b_nz

        if psi.shape  != (psi_nr,psi_nz):
            raise ValueError("Inconsistent shape for psi.")
        if br.shape   != (b_nr,b_nphi,b_nz,b_nt):
            raise ValueError("Inconsistent shape for br.")
        if bphi.shape != (b_nr,b_nphi,b_nz,b_nt):
            raise ValueError("Inconsistent shape for bphi.")
        if bz.shape   != (b_nr,b_nphi,b_nz,b_nt):
            raise ValueError("Inconsistent shape for bz.")

        psi  = np.transpose(psi)
        br   = np.transpose(br,   (3,2,1,0))
        bphi = np.transpose(bphi, (3,2,1,0))
        bz   = np.transpose(bz,   (3,2,1,0))

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("b_rmin",   (1,), data=b_rmin,   dtype="f8")
            g.create_dataset("b_rmax",   (1,), data=b_rmax,   dtype="f8")
            g.create_dataset("b_nr",     (1,), data=b_nr,     dtype="i4")
            g.create_dataset("b_phimin", (1,), data=b_phimin, dtype="f8")
            g.create_dataset("b_phimax", (1,), data=b_phimax, dtype="f8")
            g.create_dataset("b_nphi",   (1,), data=b_nphi,   dtype="i4")
            g.create_dataset("b_zmin",   (1,), data=b_zmin,   dtype="f8")
            g.create_dataset("b_zmax",   (1,), data=b_zmax,   dtype="f8")
            g.create_dataset("b_nz",     (1,), data=b_nz,     dtype="i4")
            g.create_dataset("b_tmin",   (1,), data=b_tmin,   dtype="f8")
            g.create_dataset("b_tmax",   (1,), data=b_tmax,   dtype="f8")
            g.create_dataset("b_nt",     (1,), data=b_nt,     dtype="i4")
            g.create_dataset("psi_rmin", (1,), data=psi_rmin, dtype="f8")
            g.create_dataset("psi_rmax", (1,), data=psi_rmax, dtype="f8")
            g.create_dataset("psi_nr",   (1,), data=psi_nr,   dtype="i4")
            g.create_dataset("psi_zmin", (1,), data=psi_zmin, dtype="f8")
            g.create_dataset("psi_zmax", (1,), data=psi_zmax, dtype="f8")
            g.create_dataset("psi_nz",   (1,), data=psi_nz,   dtype="i4")
            g.create_dataset("axisr",    (1,), data=axisr,    dtype="f8")
            g.create_dataset("axisz",    (1,), data=axisz,    dtype="f8")
            g.create_dataset("psi0",     (1,), data=psi0,     dtype="f8")
            g.create_dataset("psi1",     (1,), data=psi1,     dtype="f8")

            g.create_dataset("psi",  (psi_nz, psi_nr),     data=psi,  dtype="f8")
            g.create_dataset("br",   (b_nt, b_nz, b_nphi, b_nr), data=br,   dtype="f8")
            g.create_dataset("bphi", (b_nt, b_nz, b_nphi, b_nr), data=bphi, dtype="f8")
            g.create_dataset("bz",   (b_nt, b_nz, b_nphi, b_nr), data=bz,   dtype="f8")
        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read time-dependent 3D magnetic field input from HDF5 file.
        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.
        Returns:
        Dictionary containing input data.
        """

        path = "bfield" + "/B_3DST_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["psi"]  = np.transpose(out["psi"])
        out["br"]   = np.transpose(out["br"],   (3,2,1,0))
        out["bphi"] = np.transpose(out["bphi"], (3,2,1,0))
        out["bz"]   = np.transpose(out["bz"],   (3,2,1,0))
        return out

class B_STS(DataGroup):
    """
    Object representing B_STS data.
    """

    def read(self):
        return read_hdf5(self._root._ascot.file_getpath(), self.get_qid())

    def write(self, fn, data=None):
        if data is None:
            data = self.read()

        return write_hdf5(fn, **data)

    @staticmethod
    def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                   b_phimin, b_phimax, b_nphi, psi0, psi1,
                   br, bphi, bz, psi,
                   axis_phimin, axis_phimax, axis_nphi, axisr, axisz,
                   psi_rmin=None, psi_rmax=None, psi_nr=None,
                   psi_zmin=None, psi_zmax=None, psi_nz=None,
                   psi_phimin=None, psi_phimax=None, psi_nphi=None,
                   desc=None):
        """Write stellarator magnetic field input in HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        b_rmin : float
            Magnetic field data R grid min edge [m].
        b_rmax : float
            Magnetic field data R grid max edge [m].
        b_nr : int
            Number of R grid points in magnetic field data.
        b_zmin : float
            Magnetic field data z grid min edge [m].
        b_zmax : float
            Magnetic field data z grid max edge [m].
        b_nz : int
            Number of z grid points in magnetic field data.
        b_phimin : float
            Magnetic field data phi grid min edge [deg].
        b_phimax : float
            Magnetic field data phi grid max edge [deg].
        b_nphi : int
            Number of phi grid points in magnetic field data.
        axis_phimin : float
            Magnetic axis phi grid min value [deg].
        axis_phimax : float
            Magnetic axis phi grid max value [deg].
        axis_nphi : float
            Number of points in magnetic axis phi grid.
        axisr : float
            Magnetic axis R coordinates on axis phi grid [m].
        axisz : float
            Magnetic axis z coordinates on axis phi grid [m].
        psi0 : float
            On-axis poloidal flux value [Vs/m].
        psi1 : float
            Separatrix poloidal flux value [Vs/m].
        psi : array_like (psi_nr,psi_nphi,psi_nz)
            Poloidal flux values on the Rz grid [Vs/m].
        br : array_like (b_nr,b_nphi,b_nz)
            Magnetic field R component (excl. equilibrium comp.) on Rz grid [T].
        bphi : array_like (b_nr,b_nphi,b_nz)
            Magnetic field phi component on Rz grid [T].
        bz : array_like (b_nr,b_nphi,b_nz)
            Magnetic field z component (excl. equilibrium comp.) onRz grid [T].
        psi_rmin : float, optional
            Psi data R grid min edge [m].
        psi_rmax : float, optional
            Psi data R grid max edge [m].
        psi_nr : int, optional
            Number of R grid points in psi data.
        psi_zmin : float, optional
            Psi data z grid min edge [m].
        psi_zmax : float, optional
            Psi data z grid max edge [m].
        psi_nz : int, optional
            Number of z grid points in psi data.
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

        parent = "bfield"
        group  = "B_STS"
        gname  = ""

        # Define psigrid to be same as Bgrid if not stated otherwise.
        if(psi_rmin is None or psi_rmax is None or psi_nr is None
           or psi_phimin is None or psi_phimax is None or psi_nphi is None
           or psi_zmin is None or psi_zmax is None or psi_nz is None):
            psi_rmin   = b_rmin
            psi_rmax   = b_rmax
            psi_nr     = b_nr
            psi_phimin = b_phimin
            psi_phimax = b_phimax
            psi_nphi   = b_nphi
            psi_zmin   = b_zmin
            psi_zmax   = b_zmax
            psi_nz     = b_nz

        if psi.shape  != (psi_nr,psi_nphi,psi_nz):
            raise ValueError("Inconsistent shape foor psi.")
        if br.shape   != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape foor br.")
        if bphi.shape != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape foor bphi.")
        if bz.shape   != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape foor bz.")

        psi  = np.transpose(psi,  (2,1,0))
        br   = np.transpose(br,   (2,1,0))
        bphi = np.transpose(bphi, (2,1,0))
        bz   = np.transpose(bz,   (2,1,0))

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("b_rmin",      (1,),     data=b_rmin,      dtype="f8")
            g.create_dataset("b_rmax",      (1,),     data=b_rmax,      dtype="f8")
            g.create_dataset("b_nr",        (1,),     data=b_nr,        dtype="i4")
            g.create_dataset("b_phimin",    (1,),     data=b_phimin,    dtype="f8")
            g.create_dataset("b_phimax",    (1,),     data=b_phimax,    dtype="f8")
            g.create_dataset("b_nphi",      (1,),     data=b_nphi,      dtype="i4")
            g.create_dataset("b_zmin",      (1,),     data=b_zmin,      dtype="f8")
            g.create_dataset("b_zmax",      (1,),     data=b_zmax,      dtype="f8")
            g.create_dataset("b_nz",        (1,),     data=b_nz,        dtype="i4")
            g.create_dataset("psi_rmin",    (1,),     data=psi_rmin,    dtype="f8")
            g.create_dataset("psi_rmax",    (1,),     data=psi_rmax,    dtype="f8")
            g.create_dataset("psi_nr",      (1,),     data=psi_nr,      dtype="i4")
            g.create_dataset("psi_phimin",  (1,),     data=psi_phimin,  dtype="f8")
            g.create_dataset("psi_phimax",  (1,),     data=psi_phimax,  dtype="f8")
            g.create_dataset("psi_nphi",    (1,),     data=psi_nphi,    dtype="i4")
            g.create_dataset("psi_zmin",    (1,),     data=psi_zmin,    dtype="f8")
            g.create_dataset("psi_zmax",    (1,),     data=psi_zmax,    dtype="f8")
            g.create_dataset("psi_nz",      (1,),     data=psi_nz,      dtype="i4")
            g.create_dataset("axis_phimin", (1,),     data=axis_phimin, dtype="f8")
            g.create_dataset("axis_phimax", (1,),     data=axis_phimax, dtype="f8")
            g.create_dataset("axis_nphi",   (1,),     data=axis_nphi,   dtype="i4")
            g.create_dataset("psi0",        (1,),     data=psi0,        dtype="f8")
            g.create_dataset("psi1",        (1,),     data=psi1,        dtype="f8")

            g.create_dataset("axisr", (axis_nphi,), data=axisr, dtype="f8")
            g.create_dataset("axisz", (axis_nphi,), data=axisz, dtype="f8")

            g.create_dataset("br",         (b_nz,b_nphi,b_nr),       data=br,
                             dtype="f8")
            g.create_dataset("bphi",       (b_nz,b_nphi,b_nr),       data=bphi,
                             dtype="f8")
            g.create_dataset("bz",         (b_nz,b_nphi,b_nr),       data=bz,
                             dtype="f8")
            g.create_dataset("psi",        (psi_nz,psi_nphi,psi_nr), data=psi,
                             dtype="f8")

        return gname

    @staticmethod
    def read_hdf5(fn, qid):
        """
        Read stellarator magnetic field input from HDF5 file.

        Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            QID of the data to be read.

        Returns:
        Dictionary containing input data.
        """

        path = "bfield/B_STS_" + qid

        out = {}
        with h5py.File(fn,"r") as f:
            for key in f[path]:
                out[key] = f[path][key][:]

        out["psi"]  = np.transpose(out["psi"],  (2,1,0))
        out["br"]   = np.transpose(out["br"],   (2,1,0))
        out["bphi"] = np.transpose(out["bphi"], (2,1,0))
        out["bz"]   = np.transpose(out["bz"],   (2,1,0))
        return out
