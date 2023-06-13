"""Analytic tokamak field HDF5 IO.
"""
import numpy as np
import h5py
import random
import datetime
import a5py.preprocessing.analyticequilibrium as psifun
import a5py.ascot5io.B_2DS as B_2DS
import a5py.ascot5io.B_3DS as B_3DS

from ._iohelpers.fileapi import add_group

from ._iohelpers.treedata import DataGroup

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
