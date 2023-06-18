"""Analytical inputs.
"""
import numpy as np
import unyt

import a5py.ascot5io.marker as marker
import a5py.ascot5io.options as options

from a5py.ascot5io.bfield import B_GS, B_2DS, B_3DS
from a5py.ascot5io.wall import wall_2D, wall_3D

class AnalyticalInputs():
    """Inputs that can be constructed analytically or otherwise from scratch.
    """

    def bfield_analytical_iter_circular(self, axisymmetric=True, splines=False):
        """Create ITER like but circular magnetic field.

        Parameters
        ----------
        axisymmetric : bool, optional
            If True, the field contains toroidal ripple.

            Note that the 3D field is not divergence free.
        splines : bool, optional
            If True, the analytical magnetic field format class:`B_GS` is
            converted to :class:`B_2DS` or :class:`B_3DS` depending on whether
            the field is axisymmetric or not.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """

        out = {}
        out["coefficients"] = np.array([
             2.21808016e-02,  -1.28841781e-01,  -4.17718173e-02,
            -6.22680280e-02,   6.20083978e-03,  -1.20524711e-03,
            -3.70147050e-05,   0.00000000e+00,   0.00000000e+00,
             0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
            -0.155])
        out.update({"r0" : 6.2, "z0" : 0, "bphi0" : 5.3, "psimult" : 200})

        if axisymmetric:
            if splines:
                out.update({"rmin" : 4, "rmax" : 8.5, "nr" : 120, "zmin" : -4,
                            "zmax" : 4, "nz" : 200})
                out = B_GS.write_hdf5_B_2DS(fn="", **out)
                return ("B_2DS", B_2DS.convert_B_GS(**out))
            return ("B_GS", out)

        out.update({"nripple" : 18, "a0" : 2, "alpha0" : 0.2, "delta0" : 0.5})
        if splines:
            out.update({"rmin" : 4, "rmax" : 8.5, "nr" : 120, "zmin" : -4,
                        "zmax" : 4, "nz" : 200, "pmin" : 0, "pmax" : 360,
                        "nphi" : 360})
            return ("B_3DS", B_3DS.convert_B_GS(**out))

        return ("B_GS", out)

    def plasma_flat(self, density=10e20, temperature=10e3):
        """Create uniform hydrogen plasma that decays is flat inside the
        separatrix but inexistent outside.

        Parameters
        ----------
        density : float, optional
            Plasma density.
        temperature : float, optional
            Plasma temperature.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        nrho   = 100
        rho    = np.transpose( np.linspace(0, 2, nrho) )
        prof   = np.ones((nrho, 1))
        edens  = density     * prof
        etemp  = temperature * prof
        idens  = density     * prof
        itemp  = temperature * prof
        edens[rho>1] = 1
        idens[rho>1] = 1

        out = {"nrho" : nrho, "nion" : 1, "anum" : np.array([1]),
               "znum" : np.array([1]), "mass" : np.array([1]),
               "charge" : np.array([1]), "rho" : rho,
               "edensity" : edens, "etemperature" : etemp, "idensity" : idens,
               "itemperature" : itemp}
        return ("plasma_1D", out)

    def wall_rectangular(self, nphi=1):
        """Create wall with a rectangular cross section.

        Parameters
        ----------
        nphi : int, optional
            Number of sectors in 3D wall.

            If larger than one, the wall input will be :class:`wall_3D` instead
            of :class:`wall_2D`.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        out = {"nelements" : 20,
               "r" :  np.concatenate( (np.linspace(4.1, 8.4, 6)[1:],
                                       np.linspace(8.4, 8.4, 6)[1:],
                                       np.linspace(8.4, 4.1, 6)[1:],
                                       np.linspace(4.1, 4.1, 6)[1:]) ),
               "z" :  np.concatenate( (np.linspace(-3.9, -3.9, 6)[1:],
                                       np.linspace(-3.9, 3.9, 6)[1:],
                                       np.linspace(3.9, 3.9, 6)[1:],
                                       np.linspace(3.9, -3.9, 6)[1:]) )
        }
        if nphi > 1:
            out.update({"nphi" : nphi})
            return ("wall_3D", wall_3D.convert_wall_2D(**out))
        return ("wall_2D", out)
