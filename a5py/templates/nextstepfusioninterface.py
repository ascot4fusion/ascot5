"""Interface to the Next-Step Fusion simulation ecosystem.
"""
import numpy as np
import unyt

class ImportNSF():

    def nsf_wall(self, x_limiter, y_limiter):
        """Initialize 2D wall from given coordinates.

        Parameters
        ----------
        x_limiter : array_like
            x coordinates of the wall limiter [cm].
        y_limiter : array_like
            y coordinates of the wall limiter [cm].

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        out = {
            "nelements": x_limiter.size,
            "r": x_limiter / 100,
            "z": y_limiter / 100,
            }
        return ("wall_2D", out)


    def nsf_plasma(
        self,
        normalized_rho,
        electron_density_profile,
        electron_temperature_profile,
        ion_temperature_profile,
        zeff,
        maxrho=3.0,
        ):
        """Initialize plasma.

        For now it is assumed that there's a single impurity species and that is
        C12.

        Parameters
        ----------
        normalized_rho : array_like
            Normalized rho values where profiles are defined.
        electron_density_profile : array_like
            Electron density profile.
        electron_temperature_profile : array_like
            Electron temperature profile.
        ion_temperature_profile : array_like
            Ion temperature profile.
        zeff : float
            Effective charge of the plasma.
        maxrho : float, optional
            Maximum rho value to which the profiles are extrapolated.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """

        rho = np.append(normalized_rho, maxrho)

        etemperature = electron_temperature_profile
        etemperature = np.append(etemperature, etemperature[-1])
        itemperature = ion_temperature_profile
        itemperature = np.append(itemperature, itemperature[-1])

        edensity = electron_density_profile * 1e19
        edensity = np.append(edensity, edensity[-1] / 1e10)

        # Calculate impurity fraction assuming C12
        if zeff <= 1.0:
            zeff = 1.000001
        nc = (zeff - 1) * edensity / (6**2 - 6)
        nh = edensity - 6*nc
        idensity = np.zeros((rho.size, 2))
        idensity[:,0] = nh
        idensity[:,1] = nc
        out = {
            "nrho": int(rho.size),
            "nion": 2,
            "anum": np.array([2, 12]),
            "znum": np.array([1, 6]),
            "mass": np.array([1.008, 12.000]),
            "charge": np.array([1.0, 6.0]),
            "rho": rho,
            "vtor": rho * 0,
            "edensity": edensity,
            "etemperature": etemperature,
            "idensity": idensity,
            "itemperature": itemperature,
            }
        return ("plasma_1D", out)


    def nsf_bfield(
        self,
        x_psi,
        y_psi,
        magnetic_center,
        psi,
        psi_axis,
        psi_separatrix,
        b0,
        ):
        """Iniitialize magnetic field.

        Parameters
        ----------
        x_psi : array_like
            Psi grid in the R direction.
        y_psi : array_like
            Psi grid in the z direction.
        magnetic_center : array_like
            Magnetic center coordinates.
        psi : array_like
            Psi values on the grid.
        psi_axis : float
            Psi axis value.
        psi_separatrix : float
            Psi separatrix value.
        b0 : float
            Toroidal field strength at the magnetic axis.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        r = x_psi.ravel()
        bphi = np.tile(b0 * magnetic_center[0] / r, (y_psi.size,1)).T
        scale = 1e5
        out = {
            "rmin": x_psi[0],
            "rmax": x_psi[-1],
            "nr": x_psi.size,
            "zmin": y_psi[0],
            "zmax": y_psi[-1],
            "nz": y_psi.size,
            "axisr": magnetic_center[0],
            "axisz": magnetic_center[1],
            "psi": psi / scale,
            "psi0": psi_axis / scale,
            "psi1": psi_separatrix / scale,
            "br": psi * 0,
            "bphi": bphi,
            "bz": psi * 0,
            }
        return ("B_2DS", out)
