"""Biot-Savart solver for computing magnetic field from a coil geometry.
"""
import numpy as np
import ctypes
from scipy import interpolate
from copy import deepcopy

from a5py.exceptions import AscotInitException
from a5py.ascotpy.libascot import _LIBASCOT, PTR_REAL

class BioSaw():
    """Biot-Savart solver for computing magnetic field from a coil geometry.

    Attributes
    ----------
    _ascot : :class:`.Ascot`
        Ascot object used to run BioSaw.
    """

    def __init__(self, ascot):
        self._ascot = ascot

    def calculate(self, r, phi, z, coilxyz, current=1.0, revolve=None):
        """Calculate magnetic field at given points due to given coil(s).

        Parameters
        ----------
        r : array_like (nr,)
            R grid for the coordinates where the field is evaluated.
        phi : array_like (nphi,)
            phi grid for the coordinates where the field is evaluated.
        z : array_like (nz,)
            z grid for the coordinates where the field is evaluated.
        coilxyz : list or array_like, (ncoil,3)
            xyz coordinates defining the coil segments.

            This argument can also be a list of coordinate arrays, in which case
            each list element is assumed to be an individual coil and
            the resulting field is summed together.
        current : list or float
            Coil current in Ampères either same for all coils or given
            individually.
        revolve : int, optional
            If set, the field is rotated repeatedly in toroidal direction by
            a given amount of phi indices and summed, as if the coils were
            revolved.

            In other words, the total field is given by:

            B_tot = B_coil(phi0 at 0) + B_coil(phi0 at phi[revolve])
            + B_coil(phi0 at phi[2*revolve]) + ...,

            where the field is rotated through the whole ``phi`` array.

        Returns
        -------
        br : array_like, (nr,nphi,nz)
            Magnetic field R component.
        bphi : array_like, (nr,nphi,nz)
            Magnetic field phi component.
        bz : array_like, (nr,nphi,nz)
            Magnetic field z component.
        """
        if _LIBASCOT is None:
            raise AscotInitException(
                "Failed to load libascot.so which is required by BioSaw")

        if not isinstance(coilxyz, list): coilxyz = [coilxyz]
        if not isinstance(current, list): current = [current] * len(coilxyz)

        [R,P,Z] = np.meshgrid(r, phi, z, indexing="ij")
        X = R * np.cos(P)
        Y = R * np.sin(P)
        n = R.size

        bx = np.zeros(n, dtype="f8")
        by = np.zeros(n, dtype="f8")
        bz = np.zeros(n, dtype="f8")

        fun = _LIBASCOT.biosaw_calc_B
        fun.restype  = None
        fun.argtypes = [ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                        ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL,
                                      PTR_REAL, PTR_REAL, PTR_REAL]

        for xyz, I in zip(coilxyz, current):
            coiln = xyz.shape[0]
            bx0 = np.zeros(n, dtype="f8")
            by0 = np.zeros(n, dtype="f8")
            bz0 = np.zeros(n, dtype="f8")
            if any(  np.diff(xyz[:,0])**2 + np.diff(xyz[:,1])**2
                   + np.diff(xyz[:,2])**2 == 0):
                raise ValueError("A coil has zero-length elements")
            _LIBASCOT.biosaw_calc_B(n, X.ravel(), Y.ravel(), Z.ravel(),
                                    coiln, xyz[:,0], xyz[:,1], xyz[:,2],
                                    bx0, by0, bz0)
            bx += bx0 * I
            by += by0 * I
            bz += bz0 * I

        shape = (r.size,phi.size,z.size)
        bx = bx.reshape(shape)
        by = by.reshape(shape)
        bz = bz.reshape(shape)

        br   =  bx*np.cos(P) + by*np.sin(P)
        bphi = -bx*np.sin(P) + by*np.cos(P)
        bz   =  bz

        if revolve is not None:
            br0   = np.copy(br)
            bphi0 = np.copy(bphi)
            bz0   = np.copy(bz)
            i = 1
            while i*revolve < phi.size:
                br   += np.roll(br0,   i*revolve, axis=1)
                bphi += np.roll(bphi0, i*revolve, axis=1)
                bz   += np.roll(bz0,   i*revolve, axis=1)
                i += 1

        return br, bphi, bz

    def addto2d(self, coilxyz, phimin=0, phimax=360, nphi=360,
                rmin=None, rmax=None, nr=None, zmin=None, zmax=None, nz=None,
                current=1.0, revolve=None, b0=None, b2d=None):
        """Generate 3D tokamak magnetic field input from coil geometry and
           2D input.

        This function uses the active :class:`.B_2DS` input to generate
        the :class:`.B_3DS` input.

        Note that the toroidal grid is defined as:

        ``np.linspace(phimin, phimax, nphi+1)[:-1]``

        Parameters
        ----------
        coilxyz : list or array_like, (ncoil,3)
            xyz coordinates defining the coil segments.

            This argument can also be a list of coordinate arrays, in which case
            each list element is assumed to be an individual coil and
            the resulting field is summed together.
        phimin : float
            Toroidal grid minimum value.
        phimax : float
            Toroidal grid maximum value.
        nphi : int
            Number of toroidal grid points.
        rmin : float, optional
            Minimum value in R grid unless taken from the input data.
        rmax : float, optional
            Maximum value in R grid unless taken from the input data.
        nr : int, optional
            Number of R grid points unless taken from the input data.
        zmin : float, optional
            Minimum value in z grid unless taken from the input data.
        zmax : float, optional
            Maximum value in z grid unless taken from the input data.
        nz : int, optional
            Number of z grid points unless taken from the input data.
        current : list or float
            Coil current in Ampères either same for all coils or given
            individually.

            Note that the absolute value does not matter here (only relative
            values between the coils) if ``b0`` is provided.
        revolve : int, optional
            If set, the field is rotated repeatedly in toroidal direction by
            a given amount of phi indices and summed, as if the coils were
            revolved.

            In other words, the total field is given by:

            B_tot = B_coil(phi0 at 0) + B_coil(phi0 at phi[revolve])
            + B_coil(phi0 at phi[2*revolve]) + ...,

            where the field is rotated through the whole ``phi`` array.
        b0 : float, optional
            If given, the current in the coils is scaled (after revolving them
            if requested) so that the toroidal field on axis has this value on
            average.
        b2d : dict, optional
            :class:`.B_2DS` data to be used or otherwise read from active bfield
            input.

        Returns
        -------
        out : dict
            :class:`.B_3DS` data.
        """
        b3d = {}
        if b2d is None:
            b2d = self._ascot.data.bfield.active
            if b2d.get_type() != "B_2DS":
                raise ValueError("Active bfield input is not B_2DS")
            b2d = b2d.read()

        b3d["b_rmin"]   = b2d["rmin"]   if rmin   is None else rmin
        b3d["b_rmax"]   = b2d["rmax"]   if rmax   is None else rmax
        b3d["b_nr"]     = b2d["nr"]     if nr     is None else nr
        b3d["b_zmin"]   = b2d["zmin"]   if zmin   is None else zmin
        b3d["b_zmax"]   = b2d["zmax"]   if zmax   is None else zmax
        b3d["b_nz"]     = b2d["nz"]     if nz     is None else nz
        b3d["b_phimin"] = phimin
        b3d["b_phimax"] = phimax
        b3d["b_nphi"]   = nphi

        r = np.linspace(b3d["b_rmin"], b3d["b_rmax"], int(b3d["b_nr"]))
        z = np.linspace(b3d["b_zmin"], b3d["b_zmax"], int(b3d["b_nz"]))
        p = np.linspace(phimin, phimax, nphi+1)[:-1]
        br, bphi, bz = self.calculate(r, p*np.pi/180, z, coilxyz,
                                      current=current, revolve=revolve)

        # Scale axis Bphi to given value if requested
        scaling = 1.0
        if b0 is not None:
            f = interpolate.interp2d(r, z, np.mean(bphi, axis=1).T,
                                     kind='cubic')
            scaling = b0 / f(b2d["axisr"], b2d["axisz"])

        b3d["br"]   = scaling * br
        b3d["bphi"] = scaling * bphi
        b3d["bz"]   = scaling * bz
        b3d.update(
            {"axisr":b2d["axisr"], "axisz":b2d["axisz"], "psi":b2d["psi"],
             "psi0":b2d["psi0"], "psi1":b2d["psi1"], "psi_rmin":b2d["rmin"],
             "psi_rmax":b2d["rmax"], "psi_nr":b2d["nr"], "psi_zmin":b2d["zmin"],
             "psi_zmax":b2d["zmax"], "psi_nz":b2d["nz"]})
        return b3d

    def addto3d(self, coilxyz, current=1.0, revolve=None, b3d=None):
        """Generate new 3D tokamak magnetic field input where contribution from
           given coils are included.

        Parameters
        ----------
        coilxyz : list or array_like, (ncoil,3)
            xyz coordinates defining the coil segments.

            This argument can also be a list of coordinate arrays, in which case
            each list element is assumed to be an individual coil and
            the resulting field is summed together.
        current : list or float
            Coil current in Ampères either same for all coils or given
            individually.
        revolve : int, optional
            If set, the field is rotated repeatedly in toroidal direction by
            a given amount of phi indices and summed, as if the coils were
            revolved.

            In other words, the total field is given by:

            B_tot = B_coil(phi0 at 0) + B_coil(phi0 at phi[revolve])
            + B_coil(phi0 at phi[2*revolve]) + ...,

            where the field is rotated through the whole ``phi`` array.
        b3d : dict, optional
            :class:`.B_3DS` data to be used or otherwise read from active bfield
            input.

        Returns
        -------
        out : dict
            :class:`.B_3DS` data.
        """
        if b3d is None:
            b3d = self._ascot.data.bfield.active
            if b3d.get_type() != "B_3DS":
                raise ValueError("Active bfield input is not B_3DS")
            b3d = b3d.read()
        else:
            b3d = deepcopy(b3d)

        r = np.linspace(b3d["b_rmin"], b3d["b_rmax"], int(b3d["b_nr"]))
        z = np.linspace(b3d["b_zmin"], b3d["b_zmax"], int(b3d["b_nz"]))
        p = np.linspace(b3d["b_phimin"], b3d["b_phimax"], int(b3d["b_nphi"])+1)
        p = p[:-1]
        br, bphi, bz = self.calculate(r, p*np.pi/180, z, coilxyz,
                                      current=current, revolve=revolve)
        b3d["br"]   += br
        b3d["bphi"] += bphi
        b3d["bz"]   += bz

        return b3d
