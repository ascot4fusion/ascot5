import numpy as np
import unyt

from scipy.interpolate import interpn
from scipy.integrate   import quad

from a5py.physlib import aeq

from .template import InputTemplate


class TransmuteAnalyticalBfield(InputTemplate):

    def convert_B_GS(rmin, rmax, nr, zmin, zmax, nz, **kwargs):
        """Convert :class:`B_GS` input to `B_2DS` input.

        Parameters
        ----------
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
        **kwargs
            Arguments passed to :meth:`B_GS.write_hdf5` excluding ``fn`` and
            ``desc``.

        Returns
        -------
        out : dict
            :class:`B_GS` converted as an input for :meth:`write_hdf5`.
        """
        rgrid = np.linspace(rmin, rmax, nr)
        zgrid = np.linspace(zmin, zmax, nz)
        zg, rg = np.meshgrid(zgrid, rgrid)

        c = kwargs["coefficients"] # For shorter notation.
        psirz = kwargs["psimult"] * psifun.psi0(
            rg/kwargs["r0"], zg/kwargs["r0"], c[0], c[1], c[2], c[3], c[4],
            c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12])

        br = np.zeros((nr,nz))
        bz = np.zeros((nr,nz))
        bphi = ( kwargs["r0"] / rg ) * kwargs["bphi0"]

        # search for magnetic axis if not given
        if not "psi0" in kwargs:
            x = psifun.find_axis(kwargs["r0"], kwargs["z0"], c[0], c[1], c[2],
                                    c[3], c[4], c[5],
                                    c[6], c[7], c[8], c[9], c[10], c[11], c[12])
            psi0 = kwargs["psimult"]*psifun.psi0(
                x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12])
            r0 = x[0]*kwargs["r0"]
            z0 = x[1]*kwargs["r0"]
            psi1 = 0
            # Padding
            psi0 = psi0 - 1e-4 if psi0 < psi1 else psi0 + 1e-4

        return {"rmin" : rmin, "rmax" : rmax, "nr" : nr, "zmin" : zmin,
                "zmax" : zmax, "nz" : nz, "axisr" : r0, "axisz" : z0,
                "psi" : psirz, "psi0" : psi0, "psi1" : psi1, "br" : br,
                "bphi" : bphi, "bz" : bz}


    def convert_B_GS(rmin, rmax, nr, zmin, zmax, nz, phimin, phimax, nphi,
                        **kwargs):
        """Convert :class:`B_GS` input to `B_3DS` input.

        The resulting field is far from being divergence-free as the
        perturbation BR and Bz components are omitted.

        Parameters
        ----------
        rmin : float
            Magnetic field data R grid min edge [m].
        rmax : float
            Magnetic field data R grid max edge [m].
        nr : int
            Number of R grid points in magnetic field data.
        zmin : float
            Magnetic field data z grid min edge [m].
        zmax : float
            Magnetic field data z grid max edge [m].
        nz : int
            Number of z grid points in magnetic field data.
        phimin : float
            Beginning of the toroidal period [deg].
        phimax : float
            End of the toroidal period [deg].
        nphi : int
            Number of phi grid points in magnetic field data.
        **kwargs
            Arguments passed to :meth:`B_GS.write_hdf5` excluding ``fn`` and
            ``desc``.

        Returns
        -------
        out : dict
            :class:`B_GS` converted as an input for :meth:`write_hdf5`.
        """
        rgrid = np.linspace(rmin, rmax, nr)
        zgrid = np.linspace(zmin, zmax, nz)

        zg, rg = np.meshgrid(zgrid, rgrid);

        c = kwargs["coefficients"] # For shorter notation.
        psirz = kwargs["psimult"]*psifun.psi0(
            rg/kwargs["r0"],zg/kwargs["r0"],
            c[0],c[1],c[2],c[3],c[4],c[5],
            c[6],c[7],c[8],c[9],c[10],c[11],c[12])

        br = np.zeros((nr, nphi, nz))
        bz = np.zeros((nr, nphi, nz))
        bphi = np.zeros((nr, nphi, nz))

        # Ripple
        radius = np.sqrt( ( rg - kwargs["r0"] )**2 + ( zg - kwargs["z0"] )**2 )
        theta = np.arctan2( zg - kwargs["z0"], rg - kwargs["r0"] )
        delta = kwargs["delta0"] * np.exp(-0.5*theta**2) \
            * np.power( radius / kwargs["a0"], kwargs["alpha0"] )

        phigrid = np.linspace(phimin, phimax, nphi+1)[:-1]
        for i in range(nphi):
            cos = np.cos( kwargs["nripple"] * phigrid[i] * np.pi / 180 )
            bphi[:,i,:] = ( kwargs["r0"] / rg ) * kwargs["bphi0"] \
                * ( 1 + delta * cos )

        # search for magnetic axis if not given
        if not "psi0" in kwargs:
            x = psifun.find_axis(kwargs["r0"], kwargs["z0"], c[0], c[1], c[2],
                                    c[3], c[4], c[5],
                                    c[6], c[7], c[8], c[9], c[10], c[11], c[12])
            psi0 = kwargs["psimult"]*psifun.psi0(
                x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
            raxis = x[0]*kwargs["r0"]
            zaxis = x[1]*kwargs["r0"]

            psi1 = 0

        return {"b_rmin":rmin, "b_rmax":rmax, "b_nr":nr, "b_zmin":zmin,
                "b_zmax":zmax, "b_nz":nz, "b_phimin":phimin, "b_phimax":phimax,
                "b_nphi":nphi, "axisr":raxis, "axisz":zaxis, "psi":psirz,
                "psi0":psi0, "psi1":psi1, "br":br, "bphi":bphi, "bz":bz}

    def convert_B_GS(R0, z0, B_phi0, psi_mult, psi_coeff,
                        Rmin, Rmax, nR, zmin, zmax, nz, phimin, phimax, nphi,
                        psi0=None, raxis=None, zaxis=None, desc=None,
                        **kwargs):
        """Convert :class:`B_GS` input to `B_STS` input.

        The resulting field is far from being divergence-free as the
        perturbation BR and Bz components are omitted.

        Parameters
        ----------
        R0, z0 : real
            Major axis Rz-coordinates.
        B_phi0 : real
            Toroidal field at axis.
        psi_mult : real
            Scaling factor for psi.
        psi_coeff : real 13 x 1 numpy array
            Coefficients defining psi.
        Rmin, Rmax, zmin, zmax, phimin, phimax : real
            Edges of the uniform Rphiz-grid.
        nR, nz, nphi : int
            Number of Rphiz-grid points.
        psi0 : float, optional <br>
            Poloidal flux at magnetic axis.
        raxis : float, optional <br>
            Magnetic axis R coordinate [m].
        zaxis : float, optional <br>
            Magnetic axis z coordinate [m].
        **kwargs
            Arguments passed to :meth:`B_GS.write_hdf5` excluding ``fn`` and
            ``desc``.

        Returns
        -------
        out : dict
            :class:`B_GS` converted as an input for :meth:`write_hdf5`.
        """


        rgrid = np.linspace(Rmin, Rmax, nR)
        zgrid = np.linspace(zmin, zmax, nz)

        zg, Rg = np.meshgrid(zgrid, rgrid)

        c = psi_coeff # For shorter notation.
        psiRz = psi_mult*psifun.psi0(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                        c[7],c[8],c[9],c[10],c[11],c[12])
        psiRz_R = psi_mult*psifun.psiX(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                        c[7],c[8],c[9],c[10],c[11],c[12])
        psiRz_z = psi_mult*psifun.psiY(Rg/R0,zg/R0,c[0],c[1],c[2],c[3],c[4],c[5],c[6],
                                        c[7],c[8],c[9],c[10],c[11],c[12])

        Br = np.zeros((nR, nphi, nz))
        Bz = np.zeros((nR, nphi, nz))
        Bphi = np.zeros((nR, nphi, nz))
        psi  = np.zeros((nR, nphi, nz))
        axisr= np.zeros((nphi,))
        axisz= np.zeros((nphi,))


        phigrid = np.linspace(phimin,phimax,nphi+1)
        phigrid = phigrid[0:-1]

        # search for magnetic axis if not given
        if psi0 == None:
            x = psifun.find_axis(R0, z0, c[0], c[1], c[2], c[3], c[4], c[5], c[6],
                c[7], c[8], c[9], c[10], c[11], c[12])
            psi0 = psi_mult*psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12]) # At axis.
            raxis = x[0]*R0
            zaxis = x[1]*R0

        psi1=0 #??

        #Normalize psi (It tries to impersonate the VMEC psi)
        #print("Initial psi0: {}".format(psi0))
        #print("Initial psi1: {}".format(psi1))
        psiRz -= psi0 # move axis to 0
        psi1  -= psi0 # ..
        psiRz /= psi1 # scale LCFS to 1
        psi0 = 0.0
        psi1 = 1.0


        # The magnetic field one gets from psi as:
        #  B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
        #  B_z &=  \frac{1}{R}\frac{\partial\psi}{\partial R}
        #  f
        #  B_\phi = \frac{B_0R_0}{R}
        #  f

        for i in range(0,nphi):
            Bphi[:,i,:] = ((R0/Rg)*B_phi0 )
            psi[ :,i,:] = psiRz
            Br[  :,i,:] = ( -1.0 / Rg / R0 ) * psiRz_z
            Bz[  :,i,:] = (  1.0 / Rg / R0 ) * psiRz_R
            axisr[ i  ] = raxis
            axisz[ i  ] = zaxis


        return { "b_rmin"  : Rmin,   "b_rmax"  : Rmax,  "b_nr"  : nR,
                    "b_zmin"  : zmin,   "b_zmax"  : zmax,  "b_nz"  : nz,
                    "b_phimin": phimin, "b_phimax": phimax,"b_nphi": nphi,
                    "axisr" : axisr,  "axisz"   : axisz,
                    "psi" : psi, "psi0" : psi0, "psi1"  : psi1,
                    "br"      : Br,     "bphi"    : Bphi,  "bz"    : Bz,
                    "axis_phimin":phimin, "axis_phimax":phimax, "axis_nphi":nphi}


class ToroidallyAveragedBfield(InputTemplate):

    def convert_B_3DS(**kwargs):
        """Convert :class:`B_3DS` input to `B_2DS` input.

        This function takes a toroidal average of Bphi and sets BR and Bz to
        zero.

        Parameters
        ----------
        **kwargs
            Arguments passed to :meth:`B_3DS.write_hdf5` excluding ``fn`` and
            ``desc``.

        Returns
        -------
        out : dict
            :class:`B_3DS` converted as an input for :meth:`write_hdf5`.
        """
        bphi = np.mean(kwargs["bphi"], axis=1)
        br   = np.mean(kwargs["bphi"], axis=1) * 0
        bz   = np.mean(kwargs["bphi"], axis=1) * 0
        return {
            "rmin":kwargs["b_rmin"], "rmax":kwargs["b_rmax"],
            "nr":kwargs["b_nr"], "zmin":kwargs["b_zmin"],
            "zmax":kwargs["b_zmax"], "nz":kwargs["b_nz"],
            "axisr":kwargs["axisr"], "axisz":kwargs["axisz"],
            "psi0":kwargs["psi0"], "psi1":kwargs["psi1"],
            "psi":kwargs["psi"], "br":br, "bz":bz,
            "bphi":bphi}
