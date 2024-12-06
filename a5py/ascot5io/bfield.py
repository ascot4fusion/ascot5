"""Magnetic field input.

Magnetic field is present in every simulation and plays for an important role
in accuracy of the results. High-fidelity magnetic field is essential. Most of
the time magnetic field is the most memory-intensive input (and also most of
the CPU time is spent interpolating it).

The poloidal flux psi is also used as a radial coordinate in simulations e.g. to
interpolate 1D data. More precisely,
rho = sqrt( (psi - psi0) / ( psi1 - psi0 ) ), where psi0 is psi on axis and psi1
is psi on separatrix, is used. Note that psi0 and psi1 must be provided
explicitly in all magnetic field inputs. It is a good idea to include some
padding to psi0, as otherwise rho might have a complex value near separatrix
(which leads to markers being aborted there) if the given psi0 differs from the
actual extrema of the interpolation scheme.
"""
import h5py
import numpy as np
import desc.io as dscio
import desc.grid as dscg
import netCDF4 as nc
import scipy.interpolate as si

from .coreio.fileapi import add_group
from .coreio.treedata import DataGroup

import a5py.physlib.analyticequilibrium as psifun

class B_TC(DataGroup):
    """Magnetic field in Cartesian basis for testing purposes.

    This input defines the magnetic field vector on the Cartesian basis, and
    uses constant Jacobian to make the field non-uniform if needed. This field
    is only used to test the validity of the orbit-integrators.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]
        return out

    @staticmethod
    def write_hdf5(fn, bxyz, jacobian, rhoval, psival=None, axisr=1, axisz=0,
                   desc=None):
        """Write input data to the HDF5 file.

        Parameters
        ----------
        fn : str
            Full path to the HDF5 file.
        bxyz : array_like (3,1)
            Magnetic field in cartesian coordinates at origo.
        jacobian : array_like (3,3)
            Magnetic field Jacobian, jacobian[i,j] = dB_i/dx_j.
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
        if bxyz.shape != (3,1) and bxyz.shape != (3,):
            raise ValueError("Invalid shape for Bxyz.")

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
        bxyz = np.array([1,0,3])
        jac  = np.array([ [0,0,0.01], [0,0,0], [0,0,0] ])
        return {"bxyz":bxyz, "jacobian":jac, "rhoval":0.5,
                "psival":1.5, "axisr":6, "axisz":0}

class B_GS(DataGroup):
    """Analytical tokamak field.

    This field can either be axisymmetric or it can include a ripple-like
    perturbation. The field is very fast to interpolate and it is not memory
    intensive. However, it is still intended mostly for testing purposes.

    The axisymmetric field is divergence-free whereas the 3D field is not
    as it included only the toroidal component of the perturbation.

    See:
    A.J. Cerfon, J.P. Freidberg. "One size fits all" analytic solutions to
    the Grad-Shafranov equation. Physics of Plasmas 17 (3) (2010) 032502.
    http://scitation.aip.org/content/aip/journal/pop/17/3/10.1063/1.3328818
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]

        return out

    @staticmethod
    def write_hdf5(fn, r0, z0, bphi0, psimult, coefficients, psi0=None, psi1=0,
                   raxis=None, zaxis=None, nripple=0, a0=2, alpha0=2,
                   delta0=0.05, desc=None):
        """Write input data to the HDF5 file.

        If psi0 is None, magnetic axis location is found numerically and psi0
        value is evaluated at that point. The analytical field is defined by
        assuming midplane is at z=0, but it is moved to z=z0 here.

        Setting number of toroidal field coils to non-zero makes the field 3D.

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
            x = psifun.find_axis(r0, z0, c[0], c[1], c[2], c[3], c[4], c[5],
                                 c[6], c[7], c[8], c[9], c[10], c[11], c[12])
            psi0 = psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                               c[5], c[6], c[7], c[8], c[9], c[10], c[11],
                               c[12]) * psimult
            psi1  = 0
            raxis = x[0]*r0
            zaxis = x[1]*r0
            psi0 = psi0 - 1e-8 if psi0 < psi1 else psi0 + 1e-8

        with h5py.File(fn, "a") as f:
            g = add_group(f, parent, group, desc=desc)
            gname = g.name.split("/")[-1]

            g.create_dataset("r0",           (1,),   data=r0,      dtype='f8')
            g.create_dataset("z0",           (1,),   data=z0,      dtype='f8')
            g.create_dataset("raxis",        (1,),   data=raxis,   dtype='f8')
            g.create_dataset("zaxis",        (1,),   data=zaxis,   dtype='f8')
            g.create_dataset("bphi0",        (1,),   data=bphi0,   dtype='f8')
            g.create_dataset("psi0",         (1,),   data=psi0,    dtype='f8')
            g.create_dataset("psi1",         (1,),   data=psi1,    dtype='f8')
            g.create_dataset("psimult",      (1,),   data=psimult, dtype='f8')
            g.create_dataset("nripple",      (1,),   data=nripple, dtype='i8')
            g.create_dataset("a0",           (1,),   data=a0,      dtype='f8')
            g.create_dataset("alpha0",       (1,),   data=alpha0,  dtype='f8')
            g.create_dataset("delta0",       (1,),   data=delta0,  dtype='f8')
            g.create_dataset("coefficients", (13,1), data=coefficients,
                             dtype='f8')

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
        coefficients = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                                  6.200e-03, -1.205e-03, -3.701e-05,  0,
                                  0,          0,          0,          0,
                                  -0.155])
        return {"r0":6.2, "z0":0, "bphi0":5.3, "psimult":200,
                "coefficients":coefficients, "nripple":1, "a0":2, "alpha0":2,
                "delta0":0.05}


class B_2DS(DataGroup):
    """Axisymmetric tokamak field interpolated with splines.

    This field is suitable for simulations where toroidal ripple or other
    3D perturbations are not relevant (for MHD one can use the dedicated input).
    Using axisymmetric field is **much** faster than 3D fields, so using this
    input when applicable is strongly recommended.

    The input consists of BR, Bphi, Bz, and psi tabulated on an uniform
    (R,z) grid which are interpolated with cubic splines in simulation.
    During the simulation, the total BR (and Bz) are computed as a sum of input
    BR (Bz) and the component calculated from the gradient of psi. In other
    words, the equilibrium BR and Bz are already contained in psi so in most
    cases the input BR and Bz should be set to zero. If for some reason psi has
    a poor quality, it can be scaled to an insignificant value and the field can
    be provided explicitly in BR and Bz. However, in such case the field is no
    longer divergence-free.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]

        out["psi"]  = np.transpose(out["psi"])
        out["br"]   = np.transpose(out["br"])
        out["bphi"] = np.transpose(out["bphi"])
        out["bz"]   = np.transpose(out["bz"])
        return out

    @staticmethod
    def write_hdf5(fn, rmin, rmax, nr, zmin, zmax, nz,
                   axisr, axisz, psi, psi0, psi1,
                   br, bphi, bz, desc=None):
        """Write input data to the HDF5 file.

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

        if psi.shape  != (nr,nz):
            raise ValueError("Inconsistent shape for psi.")
        if br.shape   != (nr,nz):
            raise ValueError("Inconsistent shape for br.")
        if bphi.shape != (nr,nz):
            raise ValueError("Inconsistent shape for bphi.")
        if bz.shape   != (nr,nz):
            raise ValueError("Inconsistent shape for bz.")

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
    def create_dummy():
        """Create dummy data that has correct format and is valid, but can be
        non-sensical.

        This method is intended for testing purposes or to provide data whose
        presence is needed but which is not actually used in simulation.

        The dummy output is an ITER-like but circular equilibrium.

        Returns
        -------
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        # ITER-like but circular equilibrium
        coefficients = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                                  6.200e-03, -1.205e-03, -3.701e-05,  0,
                                  0,          0,          0,          0,
                                  -0.155])
        gs = {"rmin":4, "rmax":8, "nr":50, "zmin":-4, "zmax":4, "nz":100,
              "r0":6.2, "z0":0, "bphi0":5.3, "psimult":200,
              "coefficients":coefficients}
        return B_2DS.convert_B_GS(**gs)

    @staticmethod
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
        zg, rg = np.meshgrid(zgrid, rgrid);

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

    @staticmethod
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

class B_3DS(DataGroup):
    """Non-axisymmetric tokamak field.

    The input consists of BR, Bphi, and Bz tabulated on an uniform cylindrical
    grid and psi tabulated on (R,z) grid. Data is interpolated with cubic
    splines during the simulation, and the total BR (and Bz) are computed as
    a sum of input BR (Bz) and the component calculated from the gradient of
    psi. In other words, the equilibrium BR and Bz are already contained in psi
    so BR and Bz should only contain the components from 3D perturbation.

    This input is suitable for artificial 3D perturbations. For MHD, consider
    the dedicated input.

    Since this input interpolates BR and Bz with splines, the resulting B is
    not divergence free. This can be amended by using a dense grid for magnetic
    field evaluation. Note however that it can become memory-intensive.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]

        out["psi"]  = np.transpose(out["psi"])
        out["br"]   = np.transpose(out["br"],   (2,1,0))
        out["bphi"] = np.transpose(out["bphi"], (2,1,0))
        out["bz"]   = np.transpose(out["bz"],   (2,1,0))
        return out

    @staticmethod
    def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                   b_phimin, b_phimax, b_nphi,
                   axisr, axisz, psi, psi0, psi1, br, bphi, bz,
                   psi_rmin=None, psi_rmax=None, psi_nr=None,
                   psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
        """Write input data to the HDF5 file.

        It is possible to use different (R,z) grids for psi and magnetic field
        components by giving the (R,z)-grid for psi separately. 3D data can be
        memory intensive which necessitates sparser grid for B components, but
        psi can still be evaluated on a dense grid.

        The toroidal angle phi is treated as a periodic coordinate, meaning
        ``A(phi=phimin) == A(phi=phimax)``. However, the phi grid, where input
        arrays are tabulated, is ``linspace(phimin, phimax, nphi+1)[:-1]``
        to avoid storing duplicate data.

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
            Beginning of the toroidal period [deg].
        b_phimax : float
            End of the toroidal period [deg].
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
            raise ValueError("Inconsistent shape for psi.")
        if br.shape   != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape for br.")
        if bphi.shape != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape for bphi.")
        if bz.shape   != (b_nr,b_nphi,b_nz):
            raise ValueError("Inconsistent shape for bz.")

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

            g.create_dataset("psi",  (psi_nz,psi_nr),    data=psi,  dtype="f8")
            g.create_dataset("br",   (b_nz,b_nphi,b_nr), data=br,   dtype="f8")
            g.create_dataset("bphi", (b_nz,b_nphi,b_nr), data=bphi, dtype="f8")
            g.create_dataset("bz",   (b_nz,b_nphi,b_nr), data=bz,   dtype="f8")

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
        coefficients = np.array([ 2.218e-02, -1.288e-01, -4.177e-02, -6.227e-02,
                                  6.200e-03, -1.205e-03, -3.701e-05,  0,
                                  0,          0,          0,          0,
                                  -0.155])
        gs = {"rmin":4, "rmax":8, "nr":50, "zmin":-4, "zmax":4, "nz":100,
              "phimin":0, "phimax":360, "nphi":100,
              "r0":6.2, "z0":0, "bphi0":5.3, "psimult":200,
              "coefficients":coefficients, "nripple":18, "a0":2, "alpha0":2,
              "delta0":0.05}
        return B_3DS.convert_B_GS(**gs)

    @staticmethod
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

class B_3DST(DataGroup):
    """Time-dependent 3D tokamak field interpolated with splines.

    This input is equivalent to :class:`B_3DS` except that the data includes
    fourth dimension for time.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]

        out["psi"]  = np.transpose(out["psi"])
        out["br"]   = np.transpose(out["br"],   (3,2,1,0))
        out["bphi"] = np.transpose(out["bphi"], (3,2,1,0))
        out["bz"]   = np.transpose(out["bz"],   (3,2,1,0))
        return out

    @staticmethod
    def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                   b_phimin, b_phimax, b_nphi, b_tmin, b_tmax, b_nt,
                   axisr, axisz, psi, psi0, psi1, br, bphi, bz,
                   psi_rmin=None, psi_rmax=None, psi_nr=None,
                   psi_zmin=None, psi_zmax=None, psi_nz=None, desc=None):
        """Write input data to the HDF5 file.

        It is possible to use different (R,z) grids for psi and magnetic field
        components by giving the (R,z)-grid for psi separately. 3D data can be
        memory intensive which necessitates sparser grid for B components, but
        psi can still be evaluated on a dense grid.

        The toroidal angle phi is treated as a periodic coordinate, meaning
        ``A(phi=phimin) == A(phi=phimax)``. However, the phi grid, where input
        arrays are tabulated, is ``linspace(phimin, phimax, nphi+1)[:-1]``
        to avoid storing duplicate data.

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
            Beginning of the toroidal period [deg].
        b_phimax : float
            End of the toroidal period [deg].
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

            g.create_dataset("psi",  (psi_nz, psi_nr), data=psi,  dtype="f8")
            g.create_dataset("br",   (b_nt, b_nz, b_nphi, b_nr), data=br,
                             dtype="f8")
            g.create_dataset("bphi", (b_nt, b_nz, b_nphi, b_nr), data=bphi,
                             dtype="f8")
            g.create_dataset("bz",   (b_nt, b_nz, b_nphi, b_nr), data=bz,
                             dtype="f8")
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
        return {"b_rmin":4, "b_rmax":8, "b_nr":3, "b_zmin":-2, "b_zmax":2,
                "b_nz":3, "b_phimin":0, "b_phimax":360, "b_nphi":3, "b_tmin":1,
                "b_tmax":1, "b_nt":3, "axisr":4, "axisz":0, "psi0":0,
                "psi1":1, "br":np.zeros((3,3,3,3)), "bphi":np.ones((3,3,3,3)),
                "bz":np.zeros((3,3,3,3)), "psi":0.5*np.ones((3,3)),
                "psi_rmin":4, "psi_rmax":8, "psi_nr":3,
                "psi_zmin":-2, "psi_zmax":2, "psi_nz":3}

class B_STS(DataGroup):
    """Stellarator magnetic field.
    """

    def read(self):
        """Read data from HDF5 file.

        Returns
        -------
        data : dict
            Data read from HDF5 stored in the same format as is passed to
            :meth:`write_hdf5`.
        """
        out = {}
        with self as f:
            for key in f:
                out[key] = f[key][:]

        out["psi"]  = np.transpose(out["psi"],  (2,1,0))
        out["br"]   = np.transpose(out["br"],   (2,1,0))
        out["bphi"] = np.transpose(out["bphi"], (2,1,0))
        out["bz"]   = np.transpose(out["bz"],   (2,1,0))
        return out

    @staticmethod
    def vmec_field(ncfile,phimin=0,phimax=361,nphi=361,ntheta=120,
                   nr=100,nz=100,psipad=0.0):
    """Load magnetic field data from a VMEC equilibrium.

    Notes
    -----
    VMEC coordinates:
        s = psi/psi1 (normalized toroidal flux)
        u = theta (poloidal angle)
        v = phi (toroidal angle)

    VMEC outputs quantities on two different radial meshes:
        `s_full = np.linspace(0, 1, ns) # full mesh`
        `s_half = np.insert(s_full[0:-1] + 0.5 / (ns - 1),0,np.nan) # half mesh`
    The interpolation/extrapolation scheme to convert from the half mesh 
           to the full mesh
    is adapted from Hirshman et al. 1990: 
           https://doi.org/10.1016/0021-9991(90)90259-4.

    The toroidal magnetic flux is saved in the VMEC output as the variable `phi`
    
    In B_STS B_STS_eval_rho defines psi as:
        `rho[0] = sqrt( (psi - Bdata->psi0) / delta );`
    So psi is the toroidal magnetic flux.

    Values outside the LCFS are interpolated to the nearest values. Simulations
        should not extend past rho>1.

    Parameters
    ----------
    ncfile : str
        File path to VMEC NetCDF output.
    phimin : float, optional
        Minimum toroidal angle phi (deg). Default = 0.
    phimax : float, optional
        Maximum toroidal angle phi (deg). Default = 360.
    nphi : int, optional
        Number of toroidal angle phi grid points. Default = 360.
    ntheta : int, optional
        Number of poloidal angle theta grid points. Default = 120.
    nr : int, optional
        Number of radial coordinate R grid points. Default = 100.
    nz : int, optional
        Number of vertical coordinate Z grid points. Default = 100.
    psipad : float, optional
        Padding to slightly alter flux on axis. 

    Returns
    -------
    out : dict
        Dictionary with the following items:
        - `'axis_nphi'`, `'b_nphi'`, `'psi_nphi'`: nphi phi bins
        - `'b_nr'`, `'psi_nr'`: nr radial bins
        - `'b_nz'`, `'psi_nz'`: nz z-bins
        - `'axis_phimin'`, `'b_phimin'`, `'psi_phimin'`: phimin (deg)
        - `'axis_phimax'`, `'b_phimax'`, `'psi_phimax'`: (phimax-phimin)*(nphi-1)/nphi (deg)
        - `'b_rmin'`, `'psi_rmin'`: minimum radial coordinate R of output grids (m)
        - `'b_rmax'`, `'psi_rmax'`: maximum radial coordinate R of output grids (m)
        - `'b_zmin'`, `'psi_zmin'`: minimum vertical coordinate Z of output grids (m)
        - `'b_zmax'`, `'psi_zmax'`: maximum vertical coordinate Z of output grids (m)
        - `'axis_r'`: R(phi) on the magnetic axis (m)
        - `'axis_z'`: Z(phi) on the magnetic axis (m)
        - `'rlcfs'`: R(phi,len(u)) for the LCFS (m)
        - `'zlcfs'`: Z(phi,len(u)) for the LCFS (m)
        - `'psi0'`: toroidal magnetic flux on the magnetic axis (Wb)
        - `'psi1'`: toroidal magnetic flux through the last closed flux surface (Wb)
        - `'psi'`: toroidal magnetic flux psi(R,phi,Z) (Wb)
        - `'br'`: radial magnetic field B_R(R,phi,Z) (T)
        - `'bphi'`: toroidal magnetic field B_phi(R,phi,Z) (T)
        - `'bz'`: vertical magnetic field B_Z(R,phi,Z) (T)
    """
        # read VMEC NetCDF file
        data = nc.Dataset(ncfile)
        xm = np.array(data.variables["xm"])  # poloidal mode numbers
        xn = np.array(data.variables["xn"])  # toroidal mode numbers
        xn_nyq = np.array(data.variables["xn_nyq"])  # tor mode numbers (Nyquist)
        xm_nyq = np.array(data.variables["xm_nyq"])  # pol mode numbers (Nyquist)
        psi_1d = np.array(data.variables["phi"])  # toroidal flux, full mesh
        rmnc = np.array(data.variables["rmnc"])  # cos(mn) comp of cyl R, full mesh
        zmns = np.array(data.variables["zmns"])  # sin(mn) comp of cyl Z, full mesh
        bsupumnc = np.array(data.variables["bsupumnc"])  #cos(mn) comp B^u,half mesh
        bsupvmnc = np.array(data.variables["bsupvmnc"])  #cos(mn) comp B^v,half mesh

        # poloidal angle array
        theta = np.linspace(0, 2.0 * np.pi, ntheta)  # rad

        # toroidal angle array
        # note phi should start at 0 and end on 360, inclusive
        phi = np.deg2rad(np.linspace(phimin, phimax, nphi, endpoint=False))  # rad
   
        # derivatives
        rumns = rmnc * (-1 * xm)  # drmn*cos(m*u-n*v)/du = -m*rmn*sin(m*u-n*v)
        zumnc = zmns * (xm)  # dzmn*sin(m*u-n*v)/du = m*zmn*cos(m*u-n*v)
        rvmns = rmnc * (xn)  # drmn*cos(m*u-n*v)/dv = n*rmn*sin(m*u-n*v)
        zvmnc = zmns * (-1 * xn)  # dzmn*sin(m*u-n*v)/dv = -n*zmn*cos(m*u-n*v)

        # convert from half mesh to full mesh
        ns = len(psi_1d)  # number of flux surfaces
        midx_e = np.nonzero(np.mod(xm_nyq, 2) == 0)[0]  # indices of even m modes
        midx_o = np.nonzero(np.mod(xm_nyq, 2) != 0)[0]  # indices of odd m modes
        w1 = np.atleast_2d(  # interpolation weights
            0.5 * np.sqrt((np.arange(ns-2) + 2 - 1)/(np.arange(ns-2) + 2 - 1.5))).T
        w2 = np.atleast_2d(  # interpolation weights
            0.5 * np.sqrt((np.arange(ns-2) + 2 - 1)/(np.arange(ns-2) + 2 - 0.5))).T
        w3 = 2 * np.sqrt((ns - 2) / (ns - 1.5))  # extrapolation weight

        # interpolation on intermediate surfaces
        bsupumnc[1:-1, midx_e] = 0.5*bsupumnc[1:-1,midx_e] + 0.5*bsupumnc[2:,midx_e]
        bsupvmnc[1:-1, midx_e] = 0.5*bsupvmnc[1:-1,midx_e] + 0.5*bsupvmnc[2:,midx_e]
        bsupumnc[1:-1, midx_o] = w1*bsupumnc[1:-1, midx_o] + w2*bsupumnc[2:, midx_o]
        bsupvmnc[1:-1, midx_o] = w1*bsupvmnc[1:-1, midx_o] + w2*bsupvmnc[2:, midx_o]

        # extrapolation to magnetic axis
        bsupumnc[0, midx_e] = 2 * bsupumnc[1, midx_e] - bsupumnc[2, midx_e]
        bsupvmnc[0, midx_e] = 2 * bsupvmnc[1, midx_e] - bsupvmnc[2, midx_e]
        bsupumnc[0, midx_o] = 0.0
        bsupvmnc[0, midx_o] = 0.0

        # extrapolation to boundary surface
        bsupumnc[-1, midx_e] = 2 * bsupumnc[-1, midx_e] - bsupumnc[-2, midx_e]
        bsupvmnc[-1, midx_e] = 2 * bsupvmnc[-1, midx_e] - bsupvmnc[-2, midx_e]
        bsupumnc[-1, midx_o] = w3 * bsupumnc[-1, midx_o] - bsupumnc[-2, midx_o]
        bsupvmnc[-1, midx_o] = w3 * bsupvmnc[-1, midx_o] - bsupvmnc[-2, midx_o]

        # inverse Fourier transform to (s,u,v) coordinates
        r_grid = costransform(theta, phi, rmnc, xm, xn)  # R (m)
        z_grid = sintransform(theta, phi, zmns, xm, xn)  # Z (m)
        ru_grid = sintransform(theta, phi, rumns, xm, xn)  # dR/du (m/rad)
        zu_grid = costransform(theta, phi, zumnc, xm, xn)  # dZ/du (m/rad)
        rv_grid = sintransform(theta, phi, rvmns, xm, xn)  # dR/dv (m/rad)
        zv_grid = costransform(theta, phi, zvmnc, xm, xn)  # dZ/dv (m/rad)
        bu_grid = costransform(theta, phi, bsupumnc, xm_nyq, xn_nyq)  # B^u (T)
        bv_grid = costransform(theta, phi, bsupvmnc, xm_nyq, xn_nyq)  # B^v (T)

        # magnetic axis
        axis_r = r_grid[0, 0, :]  # m
        axis_z = z_grid[0, 0, :]  # m

        #get lcfs
        lcfs_r = r_grid[-1, :, :] # m
        lcfs_z = z_grid[-1, :, :] # m

        # calculate B_R, B_phi, B_z
        br_grid = bu_grid * ru_grid + bv_grid * rv_grid
        bphi_grid = bv_grid * r_grid
        bz_grid = bu_grid * zu_grid + bv_grid * zv_grid

        # range in R and Z, make interpolation arrays
        rmin = np.amin(r_grid)
        rmax = np.amax(r_grid)
        zmin = np.amin(z_grid)
        zmax = np.amax(z_grid)
        r_1d = np.linspace(rmin, rmax, nr)  # m
        z_1d = np.linspace(zmin, zmax, nz)  # m
        z_2d, r_2d = np.meshgrid(z_1d, r_1d)
    
        # toroidal magnetic flux
        psi0 = psi_1d[0]  # axis
        psi1 = psi_1d[-1]  # LCFS
        psi_2d = np.tile(psi_1d, (ntheta, 1)).T  # Wb

        # interpolate psi, B_R, B_phi, B_Z to cylindircal coordinates
        psi = np.zeros([nr, nz, nphi])
        br = np.zeros([nr, nz, nphi])
        bphi = np.zeros([nr, nz, nphi])
        bz = np.zeros([nr, nz, nphi])

        # interpolate to cylindrical grid, iterate through toroidal angle
        for i in range(nphi):
            # interpolate data inside VMEC domain
            print(i)
            psi[:, :, i] = si.griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                psi_2d.flatten(),
                (r_2d, z_2d),
                fill_value=psi1,
            )
            br[:, :, i] = si.griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                br_grid[:, :, i].flatten(),
                (r_2d, z_2d),
            )
            bphi[:, :, i] = si.griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                bphi_grid[:, :, i].flatten(),
                (r_2d, z_2d),
            )
            bz[:, :, i] = si.griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                bz_grid[:, :, i].flatten(),
                (r_2d, z_2d),
            )

            #Replace br, bphi, bz NaN values outside LCFS with closest values
            data = br[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            br[:,:,i] = filled_data
       
            data = bz[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bz[:,:,i] = filled_data

            data = bphi[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bphi[:,:,i] = filled_data

        # change order from [R,Z,phi] to [R,phi,Z]
        psi = np.transpose(psi, (0, 2, 1))
        br = np.transpose(br, (0, 2, 1))
        bphi = np.transpose(bphi, (0, 2, 1))
        bz = np.transpose(bz, (0, 2, 1))

        #pad psi0 if needed
        if psipad != 0.0:
            print('Warning: Padding psi0 with',psipad)
            psi0 += psipad

        out = {
            "axis_phimin": phimin,  # deg
            "axis_phimax": np.rad2deg(phi[-1]),  # deg
            "axis_nphi": nphi,
            "axisr": axis_r,  # m
            "axisz": axis_z,  # m
            "rlcfs": lcfs_r,  # m
            "zlcfs": lcfs_z,  # m
            "b_rmin": rmin,  # m
            "b_rmax": rmax,  # m
            "b_nr": nr,
            "b_zmin": zmin,  # m
            "b_zmax": zmax,  # m
            "b_nz": nz,
            "b_phimin": phimin,  # deg
            "b_phimax": np.rad2deg(phi[-1]),  # deg
            "b_nphi": nphi,
            "br": br,  # T
            "bphi": bphi,  # T
            "bz": bz,  # T
            "psi": psi,  # Wb
            "psi0": psi0,  # Wb
            "psi1": psi1,  # Wb
            "psi_rmin": rmin,  # m
            "psi_rmax": rmax,  # m
            "psi_nr": nr,
            "psi_zmin": zmin,  # m
            "psi_zmax": zmax,  # m
            "psi_nz": nz,
            "psi_phimin": phimin,  # deg
            "psi_phimax": np.rad2deg(phi[-1]),  # deg
            "psi_nphi": nphi,
        }

        return out

    @staticmethod
    def desc_field(h5file,phimin=0,phimax=361,nphi=361,ntheta=120,
                   nr=100,nz=100,psipad=0.0):
    """Load magnetic field data from a DESC equilibrium.

    Notes
    -----
    The toroidal magnetic flux is saved in the DESC output as the variable `Psi`
    
    In B_STS B_STS_eval_rho defines psi as:
        `rho[0] = sqrt( (psi - Bdata->psi0) / delta );`
    So psi is the toroidal magnetic flux.

    Values outside the LCFS are interpolated to the nearest values. Simulations
        should not extend past rho>1.

    Parameters
    ----------
    h5file : str
        File path to DESC HDF5 output.
    phimin : float, optional
        Minimum toroidal angle phi (deg). Default = 0.
    phimax : float, optional
        Maximum toroidal angle phi (deg). Default = 360.
    nphi : int, optional
        Number of toroidal angle phi grid points. Default = 360.
    ntheta : int, optional
        Number of poloidal angle theta grid points. Default = 120.
    nr : int, optional
        Number of radial coordinate R grid points. Default = 100.
    nz : int, optional
        Number of vertical coordinate Z grid points. Default = 100.
    psipad : float, optional
        Padding to slightly alter flux on axis.

    Returns
    -------
    out : dict
        Dictionary with the following items:
        - `'axis_nphi'`, `'b_nphi'`, `'psi_nphi'`: nphi phi bins
        - `'b_nr'`, `'psi_nr'`: nr radial bins
        - `'b_nz'`, `'psi_nz'`: nz z-bins
        - `'axis_phimin'`, `'b_phimin'`, `'psi_phimin'`: phimin (deg)
        - `'axis_phimax'`, `'b_phimax'`, `'psi_phimax'`: (phimax-phimin)*(nphi-1)/nphi (deg)
        - `'b_rmin'`, `'psi_rmin'`: minimum radial coordinate R of output grids (m)
        - `'b_rmax'`, `'psi_rmax'`: maximum radial coordinate R of output grids (m)
        - `'b_zmin'`, `'psi_zmin'`: minimum vertical coordinate Z of output grids (m)
        - `'b_zmax'`, `'psi_zmax'`: maximum vertical coordinate Z of output grids (m)
        - `'axis_r'`: R(phi) on the magnetic axis (m)
        - `'axis_z'`: Z(phi) on the magnetic axis (m)
        - `'rlcfs'`: R(phi,len(u)) for the LCFS (m)
        - `'zlcfs'`: Z(phi,len(u)) for the LCFS (m)
        - `'psi0'`: toroidal magnetic flux on the magnetic axis (Wb)
        - `'psi1'`: toroidal magnetic flux through the last closed flux surface (Wb)
        - `'psi'`: toroidal magnetic flux psi(R,phi,Z) (Wb)
        - `'br'`: radial magnetic field B_R(R,phi,Z) (T)
        - `'bphi'`: toroidal magnetic field B_phi(R,phi,Z) (T)
        - `'bz'`: vertical magnetic field B_Z(R,phi,Z) (T)
    """
        fam = dscio.load(h5file)
        try:  # if file is an EquilibriaFamily, use final Equilibrium
            eq = fam[-1]
        except:  # file is already an Equilibrium
            eq = fam
        eq.resolution_summary()

        # poloidal angle array
        theta = np.linspace(0, 2.0 * np.pi, ntheta)  # rad

        # toroidal angle array
        # note: phi should start at 0 and end on 360, inclusive
        phi = np.deg2rad(np.linspace(phimin, phimax, nphi, endpoint=False))  # rad

        # magnetic axis
        grid_axis = dscg.LinearGrid(rho=0.0, zeta=phi, NFP=eq.NFP)
        data_axis = eq.compute(["R", "Z"], grid=grid_axis)
        axis_r = data_axis["R"]  # m
        axis_z = data_axis["Z"]  # m
        psi0 = 0  # Wb

        # boundary
        grid_bdry = dscg.LinearGrid(rho=1.0, theta=theta, zeta=phi, NFP=eq.NFP)
        data_bdry = eq.compute(["R", "Z"], grid=grid_bdry)
        bdry_r = data_bdry["R"]
        bdry_z = data_bdry["Z"]

        # boundary
        grid_bdry = dscg.LinearGrid(rho=1.0, theta=theta, zeta=phi, NFP=eq.NFP)
        data_bdry = eq.compute(["R", "Z"], grid=grid_bdry)
        rmin = np.min(data_bdry["R"])  # m
        rmax = np.max(data_bdry["R"])  # m
        zmin = np.min(data_bdry["Z"])  # m
        zmax = np.max(data_bdry["Z"])  # m
        psi1 = eq.Psi  # Wb

        # output domain
        R_1d = np.linspace(rmin, rmax, nr)  # m
        Z_1d = np.linspace(zmin, zmax, nz)  # m
        Z_2d, R_2d = np.meshgrid(Z_1d, R_1d)

        # interpolate psi, B_R, B_phi, B_Z to cylindircal coordinates
        psi = np.zeros([nr, nz, nphi])
        br = np.zeros([nr, nz, nphi])
        bphi = np.zeros([nr, nz, nphi])
        bz = np.zeros([nr, nz, nphi])

        # interpolate to cylindrical grid, iterate through toroidal angle
        for k in range(nphi):
            print(k)
            grid = dscg.ConcentricGrid(
                L=eq.L_grid, M=eq.M_grid, N=0, NFP=eq.NFP, node_pattern="linear")
            grid._nodes[:, 2] = phi[k]
            data = eq.compute(["R", "Z", "psi", "B_R", "B_phi", "B_Z"], grid=grid)
            
            # interpolate data inside DESC domain
            psi[:, :, k] = si.griddata(
                (data["R"], data["Z"]),
                data["psi"] * 2 * np.pi,  # DESC `psi` is normalized by 2 pi
                (R_2d, Z_2d),
                fill_value=psi1,
            )
            br[:,:,k] = si.griddata((data["R"],data["Z"]),data["B_R"],(R_2d,Z_2d))
            bphi[:,:,k]=si.griddata((data["R"],data["Z"]),data["B_phi"],(R_2d,Z_2d))
            bz[:,:,k] = si.griddata((data["R"],data["Z"]),data["B_Z"], (R_2d,Z_2d))

            #Replace br, bphi, bz NaN values outside LCFS with closest values
            data = br[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            br[:,:,i] = filled_data

            data = bz[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bz[:,:,i] = filled_data

            data = bphi[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bphi[:,:,i] = filled_data

        # change order from [R,Z,phiang] to [R,phiang,Z]
        psi = np.transpose(psi, (0, 2, 1))
        br = np.transpose(br, (0, 2, 1))
        bphi = np.transpose(bphi, (0, 2, 1))
        bz = np.transpose(bz, (0, 2, 1))

        #pad psi0 if needed
        if psipad != 0.0:
            print('Warning: Padding psi0 with',psipad)
            psi0 += psipad

        out = {
            "axis_phimin": phimin,  # deg
            "axis_phimax": np.rad2deg(phi[-1]),  # deg
            "axis_nphi": nphi,
            "axisr": axis_r,  # m
            "axisz": axis_z,  # m
            "rlcfs": bdry_r,  # m
            "zlcfs": bdry_z,  # m
            "b_rmin": rmin,  # m
            "b_rmax": rmax,  # m
            "b_nr": nr,
            "b_zmin": zmin,  # m
            "b_zmax": zmax,  # m
            "b_nz": nz,
            "b_phimin": phimin,  # deg
            "b_phimax": np.rad2deg(phi[-1]),  # deg
            "b_nphi": nphi,
            "br": br,  # T
            "bphi": bphi,  # T
            "bz": bz,  # T
            "psi": psi,  # Wb
            "psi0": psi0,  # Wb
            "psi1": psi1,  # Wb
            "psi_rmin": rmin,  # m
            "psi_rmax": rmax,  # m
            "psi_nr": nr,
            "psi_zmin": zmin,  # m
            "psi_zmax": zmax,  # m
            "psi_nz": nz,
            "psi_phimin": phimin,  # deg
            "psi_phimax": np.rad2deg(phi[-1]),  # deg
            "psi_nphi": nphi,
        }

        return out

    @staticmethod
    def vmec_vacuum(ncfile,extfile,ntheta=120,psipad=0.0):
    """Load magnetic field data from VMEC and EXTENDER outputs into ASCOT input format.

     Notes
    -----
    The toroidal magnetic flux is saved in the VMEC output as the variable `phi`
    
    In B_STS B_STS_eval_rho defines psi as:
        `rho[0] = sqrt( (psi - Bdata->psi0) / delta );`
    So psi is the toroidal magnetic flux.

    Psi values outside the LCFS are defined as psi1, so while markers may exist in the
        vacuum region (beyond LCFS) the mapping to the rho variable will not be valid
        as the rho=1 everywhere in the vacuum region. 

    The magnetic field (Br,Bphi,Bz) is defined via the EXTENDER input file while the
        toroidal flux (psi) is supplied via the VMEC input. Care should be taken to 
        make sure that the VMEC psi calculations agree with EXTENDER for rho<1. 

    EXTENDER input is assumed to cover one field period!!

    Parameters
    ----------
    ncfile : str
        File path to VMEC NetCDF output (free-boundary solution).
    extfile : str
        File path to EXTENDER NetCDF output that corresponds to the VMEC solution.
    ntheta : int, optional
        Number of poloidal angle theta grid points. Default = 120.
    psipad : float, optional
        Padding to slightly alter flux on axis.

    Returns
    -------
    out : dict
        Dictionary with the following items:
        - `'axis_nphi'`, `'b_nphi'`, `'psi_nphi'`: nphi phi bins
        - `'b_nr'`, `'psi_nr'`: nr radial bins
        - `'b_nz'`, `'psi_nz'`: nz z-bins
        - `'axis_phimin'`, `'b_phimin'`, `'psi_phimin'`: phimin (deg)
        - `'axis_phimax'`, `'b_phimax'`, `'psi_phimax'`: (phimax-phimin)*(nphi-1)/nphi (deg)
        - `'b_rmin'`, `'psi_rmin'`: minimum radial coordinate R of output grids (m)
        - `'b_rmax'`, `'psi_rmax'`: maximum radial coordinate R of output grids (m)
        - `'b_zmin'`, `'psi_zmin'`: minimum vertical coordinate Z of output grids (m)
        - `'b_zmax'`, `'psi_zmax'`: maximum vertical coordinate Z of output grids (m)
        - `'axis_r'`: R(phi) on the magnetic axis (m)
        - `'axis_z'`: Z(phi) on the magnetic axis (m)
        - `'rlcfs'`: R(phi,len(u)) for the LCFS (m)
        - `'zlcfs'`: Z(phi,len(u)) for the LCFS (m)
        - `'psi0'`: toroidal magnetic flux on the magnetic axis (Wb)
        - `'psi1'`: toroidal magnetic flux through the last closed flux surface (Wb)
        - `'psi'`: toroidal magnetic flux psi(R,phi,Z) (Wb)
        - `'br'`: radial magnetic field B_R(R,phi,Z) (T)
        - `'bphi'`: toroidal magnetic field B_phi(R,phi,Z) (T)
        - `'bz'`: vertical magnetic field B_Z(R,phi,Z) (T)
    """
        # load NetCDF data
        vmec = nc.Dataset(ncfile)
        extender = nc.Dataset(extfile)
    
        # array dimenstions
        nr = int(extender.variables["ir"].getValue())
        nz = int(extender.variables["jz"].getValue())
        nphi = int(extender.variables["kp"].getValue())
        nfp = int(extender.variables["nfp"].getValue())
   
        # coordinate bounds (m)
        rmin = float(extender.variables["rmin"].getValue())
        rmax = float(extender.variables["rmax"].getValue())
        zmin = float(extender.variables["zmin"].getValue())
        zmax = float(extender.variables["zmax"].getValue())
    
        # magnetic field (T)
        br = np.array(extender.variables["br_001"])
        bphi = np.array(extender.variables["bp_001"])
        bz = np.array(extender.variables["bz_001"])

        # poloidal and toroidal angles (rad)
        theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=True)
        phi = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=False)
    
        # VMEC spectral coefficients
        xm = np.array(vmec.variables["xm"])  # poloidal mode numbers
        xn = np.array(vmec.variables["xn"])  # toroidal mode numbers
        rmnc = np.array(vmec.variables["rmnc"])  # cos(mn) comp of cyl R, full mesh
        zmns = np.array(vmec.variables["zmns"])  # sin(mn) comp of cyl Z, full mesh
        psi_1d = np.array(vmec.variables["phi"])  # toroidal flux, full mesh

        # inverse Fourier transform to (psi,theta,phi) coordinates (m)
        r_grid = costransform(theta, phi, rmnc, xm, xn)
        z_grid = sintransform(theta, phi, zmns, xm, xn)

        # magnetic axis (m)
        axisr = r_grid[0, 0, :]
        axisz = z_grid[0, 0, :]

        #get lcfs (m)
        lcfs_r = r_grid[-1, :, :] 
        lcfs_z = z_grid[-1, :, :]

        # toroidal magnetic flux (Wb)
        psi0 = psi_1d[0] #axis
        psi1 = psi_1d[-1] #LCFS
        psi_2d = np.tile(psi_1d, (ntheta, 1)).T
    
        # R,Z interpolation points (m)
        r_1d = np.linspace(rmin, rmax, nr)
        z_1d = np.linspace(zmin, zmax, nz)
        z_2d, r_2d = np.meshgrid(z_1d, r_1d)

        # interpolate psi inside VMEC domain to cylindircal coordinates
        psi = np.zeros([nr, nz, nphi])
        for i in range(nphi):
            print(i)
            psi[:, :, i] = si.griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                psi_2d.flatten(),
                (r_2d, z_2d),
                fill_value=psi1,
            )
    
        # close NetCDF data
        vmec.close()
        extender.close()
    
        # change order from [R,Z,phi] to [phi,Z,R] for below
        psi = np.transpose(psi, (2, 1, 0))
    
        # repeat for each field period
        phidum = phi
        for i in range(0,nfp):
            phi = np.append(phi,phidum+i*2*np.pi/nfp)
        axisr = np.tile(axisr, nfp)
        axisz = np.tile(axisz, nfp)
        lcfs_r = np.tile(lcfs_r,(1,nfp))
        lcfs_z = np.tile(lcfs_z,(1,nfp))
        br = np.tile(br, (nfp, 1, 1))
        bphi = np.tile(bphi, (nfp, 1, 1))
        bz = np.tile(bz, (nfp, 1, 1))
        psi = np.tile(psi, (nfp, 1, 1))
   
        # repeat endpoint phi=0 == phi=360
        phi = np.append(phi,2*np.pi)
        nphi = len(phi)
        axisr = np.append(axisr, axisr[0])
        axisz = np.append(axisz, axisz[0])
        br = np.append(br, br[0, :, :][None, :], axis=0)
        bphi = np.append(bphi, bphi[0, :, :][None, :], axis=0)
        bz = np.append(bz, bz[0, :, :][None, :], axis=0)
        psi = np.append(psi, psi[0, :, :][None, :], axis=0)
    
        #dumb way to append to lcfs endpoint
        dumx,dumy = lcfs_r.shape
        new_lcfsr = np.zeros([dumx,dumy+1])
        new_lcfsr[:,0:dumy] = lcfs_r
        new_lcfsr[:,-1] = lcfs_r[:,0]
        new_lcfsz = np.zeros([dumx,dumy+1])
        new_lcfsz[:,0:dumy] = lcfs_z
        new_lcfsz[:,-1] = lcfs_z[:,0]
        lcfs_r = new_lcfsr
        lcfs_z = new_lcfsz

        # change order from [phi,Z,R] to [R,phi,Z]
        br = np.transpose(br, (2, 0, 1))
        bphi = np.transpose(bphi, (2, 0, 1))
        bz = np.transpose(bz, (2, 0, 1))
        psi = np.transpose(psi, (2, 0, 1))

        #pad psi0 if needed
        if psipad != 0.0:
            print('Warning: Padding psi0 with',psipad)
            psi0 += psipad
    
        out = {
            "axis_phimin": 0,  # deg
            "axis_phimax": np.rad2deg(phi[-1]),  # deg
            "axis_nphi": nphi,
            "axisr": axisr,  # m
            "axisz": axisz,  # m
            "rlcfs": lcfs_r, # m
            "zlcfs": lcfs_z, #m
            "b_rmin": rmin,  # m
            "b_rmax": rmax,  # m
            "b_nr": nr,
            "b_zmin": zmin,  # m
            "b_zmax": zmax,  # m
            "b_nz": nz,
            "b_phimin": 0,  # deg
            "b_phimax": np.rad2deg(phi[-1]),  # deg
            "b_nphi": nphi,
            "br": br,  # T
            "bphi": bphi,  # T
            "bz": bz,  # T
            "psi": psi,  # Wb
            "psi0": psi0,  # Wb
            "psi1": psi1,  # Wb
            "psi_rmin": rmin,  # m
            "psi_rmax": rmax,  # m
            "psi_nr": nr,
            "psi_zmin": zmin,  # m
            "psi_zmax": zmax,  # m
            "psi_nz": nz,
            "psi_phimin": 0,  # deg
            "psi_phimax": np.rad2deg(phi[-1]),  # deg
            "psi_nphi": nphi,
        }

        return out

    @staticmethod
    def write_hdf5(fn, b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                   b_phimin, b_phimax, b_nphi, psi0, psi1,
                   br, bphi, bz, psi,
                   axis_phimin, axis_phimax, axis_nphi, axisr, axisz,
                   psi_rmin=None, psi_rmax=None, psi_nr=None,
                   psi_zmin=None, psi_zmax=None, psi_nz=None,
                   psi_phimin=None, psi_phimax=None, psi_nphi=None,
                   desc=None):
        """Write input data to the HDF5 file.

        It is possible to use different (R,z) grids for psi and magnetic field
        components by giving the (R,z)-grid for psi separately. 3D data can be
        memory intensive which necessitates sparser grid for B components, but
        psi can still be evaluated on a dense grid.

        The toroidal angle phi is treated as a periodic coordinate, meaning
        ``A(phi=phimin) == A(phi=phimax)``. However, the phi grid, where input
        arrays are tabulated, is ``linspace(phimin, phimax, nphi+1)[:-1]``
        to avoid storing duplicate data.

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
            Beginning of the toroidal period [deg].
        b_phimax : float
            End of the toroidal period [deg].
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
        psi_phimin : float, optional
            Psi data beginning of the toroidal period.
        psi_phimax : float, optional
            Psi data end of the toroidal period.
        psi_nphi : int, optional
            Number of phi grid points in psi data.
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

            g.create_dataset("b_rmin",      (1,), data=b_rmin,      dtype="f8")
            g.create_dataset("b_rmax",      (1,), data=b_rmax,      dtype="f8")
            g.create_dataset("b_nr",        (1,), data=b_nr,        dtype="i4")
            g.create_dataset("b_phimin",    (1,), data=b_phimin,    dtype="f8")
            g.create_dataset("b_phimax",    (1,), data=b_phimax,    dtype="f8")
            g.create_dataset("b_nphi",      (1,), data=b_nphi,      dtype="i4")
            g.create_dataset("b_zmin",      (1,), data=b_zmin,      dtype="f8")
            g.create_dataset("b_zmax",      (1,), data=b_zmax,      dtype="f8")
            g.create_dataset("b_nz",        (1,), data=b_nz,        dtype="i4")
            g.create_dataset("psi_rmin",    (1,), data=psi_rmin,    dtype="f8")
            g.create_dataset("psi_rmax",    (1,), data=psi_rmax,    dtype="f8")
            g.create_dataset("psi_nr",      (1,), data=psi_nr,      dtype="i4")
            g.create_dataset("psi_phimin",  (1,), data=psi_phimin,  dtype="f8")
            g.create_dataset("psi_phimax",  (1,), data=psi_phimax,  dtype="f8")
            g.create_dataset("psi_nphi",    (1,), data=psi_nphi,    dtype="i4")
            g.create_dataset("psi_zmin",    (1,), data=psi_zmin,    dtype="f8")
            g.create_dataset("psi_zmax",    (1,), data=psi_zmax,    dtype="f8")
            g.create_dataset("psi_nz",      (1,), data=psi_nz,      dtype="i4")
            g.create_dataset("axis_phimin", (1,), data=axis_phimin, dtype="f8")
            g.create_dataset("axis_phimax", (1,), data=axis_phimax, dtype="f8")
            g.create_dataset("axis_nphi",   (1,), data=axis_nphi,   dtype="i4")
            g.create_dataset("psi0",        (1,), data=psi0,        dtype="f8")
            g.create_dataset("psi1",        (1,), data=psi1,        dtype="f8")

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
        return {"b_rmin":4, "b_rmax":8, "b_nr":3, "b_zmin":-2, "b_zmax":2,
                "b_nz":3, "b_phimin":0, "b_phimax":360, "b_nphi":3, "psi0":0.5,
                "psi1":1, "br":np.zeros((3,3,3)), "bphi":np.ones((3,3,3)),
                "bz":np.zeros((3,3,3)), "psi":0.5*np.ones((3,3,3)),
                "axis_phimin":0, "axis_phimax":360, "axis_nphi":3,
                "axisr":np.array([6,6,6]), "axisz":np.array([0,0,0]),
                "psi_rmin":4, "psi_rmax":8, "psi_nr":3,
                "psi_zmin":-2, "psi_zmax":2, "psi_nz":3,
                "psi_phimin":0, "psi_phimax":360, "psi_nphi":3}

    @staticmethod
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
