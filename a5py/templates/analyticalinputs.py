"""Analytical inputs.
"""
import numpy as np
import unyt

from scipy.interpolate import interpn
from scipy.integrate   import quad

from a5py.ascot5io.bfield import B_2DS, B_3DS
from a5py.ascot5io.wall import wall_3D

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
                return ("B_2DS", B_2DS.convert_B_GS(**out))
            return ("B_GS", out)

        out.update({"nripple" : 18, "a0" : 2, "alpha0" : 0.2, "delta0" : 0.5})
        if splines:
            out.update({"rmin" : 4, "rmax" : 8.5, "nr" : 120, "zmin" : -4,
                        "zmax" : 4, "nz" : 200, "phimin" : 0, "phimax" : 360,
                        "nphi" : 360})
            return ("B_3DS", B_3DS.convert_B_GS(**out))

        return ("B_GS", out)

    def bfield_analytical_step_ripple(self, b2d=None, ncoil=12, rcoil=8.0,
                                      nphi=180, nperiod=1):
        """Create magnetic field ripple assuming rectangular TF-coils as in STEP.

        Parameters
        ----------
        b2d : dict
            Dictionary containing B2DS data.
        ncoil : int, optional
            Number of TF-coils.
        rcoil : float, optional
            Coil width [m].
        nphi : int, optional
            Number of toroidal grid points in the output field.
        nperiod : int, optional
            Assume that the field has 360 deg / nperiod periodicity.

            This value is used to generate output data that covers only a single
            toroidal period and which has `nphi` grid points.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        out = {}
        if b2d is None:
            raise ValueError("Provide dictionary for 'b2d' containg B2DS data")

        r0  = b2d['axisr']
        z0  = b2d['axisz']
        r   = np.linspace(b2d['rmin'], b2d['rmax'], int(b2d['nr'][0])).ravel()
        z   = np.linspace(b2d['zmin'], b2d['zmax'], int(b2d['nz'][0])).ravel()
        phi = np.linspace(0, 2*np.pi/nperiod, nphi+1)[:-1]
        b0  = interpn((r,z), b2d['bphi'], (r0,z0))

        R, P, Z = np.meshgrid(r, phi, z, indexing='ij')

        dbr   = ( b0 * r0 / R ) * np.power(R/rcoil, ncoil) * np.sin(ncoil * P)
        dbphi = ( b0 * r0 / R ) * np.power(R/rcoil, ncoil) * np.cos(ncoil * P)
        dbz   = 0 * P

        # Add the 2D components
        def to3d(comp):
            return np.transpose(
                np.multiply.outer(comp, np.ones(phi.shape)), (0,2,1))

        out['br']   = to3d(b2d['br'])   + dbr
        out['bphi'] = to3d(b2d['bphi']) + dbphi
        out['bz']   = to3d(b2d['bz'])   + dbz

        # Populate remaining data
        out['axisr']    = r0
        out['axisz']    = z0
        out['psi0']     = b2d['psi0']
        out['psi1']     = b2d['psi1']
        out['psi']      = b2d['psi']
        out['b_rmin']   = r[0]
        out['b_rmax']   = r[-1]
        out['b_nr']     = r.size
        out['b_phimin'] = 0.0
        out['b_phimax'] = 360.0 / nperiod
        out['b_nphi']   = phi.size
        out['b_zmin']   = z[0]
        out['b_zmax']   = z[-1]
        out['b_nz']     = z.size

        return ("B_3DS", out)

    def plasma_flat(self, density=10e20, temperature=10e3, anum=1, znum=1,
                    charge=1, mass=1):
        """Create uniform single-species plasma that is flat inside the
        separatrix but inexistent outside.

        Parameters
        ----------
        density : float, optional
            Plasma density.
        temperature : float, optional
            Plasma temperature.
        anum : int, optional
            Ion atomic number.
        znum : int, optional
            Ion charge number.
        charge : int, optional
            Ion charge.
        mass : int, optional
            Ion mass.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        nrho   = 100
        rho    = np.transpose( np.linspace(0, 10, nrho) )
        prof   = np.ones((nrho, 1))
        edens  = density     * prof
        etemp  = temperature * prof
        idens  = density     * prof
        itemp  = temperature * prof
        edens[rho>1] = 1
        idens[rho>1] = 1

        out = {"nrho" : nrho, "nion" : 1, "anum" : np.array([anum]),
               "znum" : np.array([znum]), "mass" : np.array([mass]),
               "charge" : np.array([charge]), "rho" : rho,
               "edensity" : edens, "etemperature" : etemp, "idensity" : idens,
               "itemperature" : itemp}
        return ("plasma_1D", out)

    def neutral_flat(self, density=10e20, temperature=10e3, anum=1, znum=1):
        """Create uniform single-species constant neutral data.
        """
        density     = np.ones((100,1)) * density
        temperature = np.ones((100,1)) * temperature
        out = {"rhomin" : 0, "rhomax" : 10, "nrho" : 100, "nspecies" : 1,
               "anum" : np.array([anum]), "znum" : np.array([znum]),
               "density" : density, "temperature" : temperature,
               "maxwellian" : 1}
        return ("N0_1D", out)

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

    def asigma_chebyshev_cx_hh0(
            self, z1=1, a1=2, m1=3.34358377e-27*unyt.kg, z2=1, a2=2,
            m2=3.34449439e-27*unyt.kg, emin=5e1*unyt.eV, emax=1.5e5*unyt.eV,
            nekin=200, tmin=1e2*unyt.eV, tmax=1e4*unyt.eV, ntemp=20):
        """Create CX data based on cross-sections from Chebyshev fit for H-H CX
        reaction.

        Calculate a set of cross-section data based on an analytical Chebyshev
        fit to data for the hydrogen-hydrogen charge-exchange reaction as
        described in the report ORNL-6086 by Barnett. In the report, this
        reaction is called, Electron Capture Cross Sections for
        H+ + H --> H + H+, and is found on page A-22.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        # Check that wanted energy span is possible and construct abscissae
        E_min = 1.2E-01 # eV/amu
        E_max = 6.3E+05
        if emin / a1 < E_min or emax / a1 > E_max:
            raise ValueError(
                f"Energy span must be between {E_min} and {E_max} eV/amu")

        # Prepare for integration
        def sigma(ekin):
            """Calculates sigma(eperamu) [m^2]"""
            if ekin < E_min or ekin > E_max:
                return 0.0
            X = ( 2*np.log(ekin)  - np.log(E_min) - np.log(E_max) ) \
                / ( np.log(E_max) - np.log(E_min) )
            T = np.array([
                -5.49142   * (     X ),
                -3.42948   * (   2*X**2 -   1 ),
                -1.98377   * (   4*X**3 -   3*X ),
                -0.878009  * (   8*X**4 -   8*X**2 +   1 ),
                -0.198932  * (  16*X**5 -  20*X**3 +   5*X ),
                0.0837431 * (  32*X**6 -  48*X**4 +  18*X**2 -  1 ),
                0.121252  * (  64*X**7 - 112*X**5 +  56*X**3 -  7*X ),
                0.0827182 * ( 128*X**8 - 256*X**6 + 160*X**4 - 32*X**2 + 1 )
                ])
            # Convert cm^2 to m^2
            return 1e-4*np.exp(-72.6656 / 2 + np.sum(T)) * unyt.m**2

        def integrand(u, ekin, temperature):
            """Calculates the integrand as a function of u.

            Before using this, ekin and temperature should be parameterized.
            """
            u *= unyt.m / unyt.s

            v1 = np.sqrt(2*ekin/m1).to('m/s')
            energycollperamu = (0.5 * m1 * u**2 / a1).to('eV')
            sigmaval = sigma(energycollperamu)
            if sigmaval == 0:
                return 0.0

            # Calculate the value of the integrand
            e1 = -m2 * v1**2 / ( 2*temperature )
            e2 = -m2 * u**2 / ( 2*temperature )
            hyperbolicsinarg =  m2 * u * v1 / temperature
            val = sigmaval * u**2 * np.sqrt( m2 / ( 2 * np.pi * temperature ) )\
                 / v1
            if(hyperbolicsinarg >= 10.0):
                # For x >= 10, the relative difference between sinh(x) and 1/2*exp(x)
                # is less than 2e-9, allowing us to use the approximation
                return np.exp(hyperbolicsinarg + e1 + e2) * val.to('m**2')
            else:
                return (2 * np.sinh(hyperbolicsinarg)* np.exp(e1 + e2) * val).to('m**2')

        egrid = np.linspace(emin, emax, nekin)
        tgrid = np.linspace(tmin,tmax,ntemp)
        sigmav = np.zeros((ntemp*nekin,1))

        # Loop through all (energy,temperature) grid points and calculate
        # the corresponding rate coefficient value
        for itemp in range(ntemp):
            for iekin in range(nekin):
                sigmav[itemp*nekin+iekin] = \
                    quad(lambda u: integrand(u, egrid[iekin], tgrid[itemp]),
                         0, 2e7, points=(4800,1e7))[0]

        data = {
            "nreac":1, "z1":np.array([z1]), "a1":np.array([a1]),
            "z2":np.array([z2]), "a2":np.array([a2]), "reactype":np.array([6]),
            "nenergy":np.array([nekin]), "energymin":np.array([emin]),
            "energymax":np.array([emax]), "ndensity":np.array([1]),
            "densitymin":np.array([1]), "densitymax":np.array([1e30]),
            "ntemperature":np.array([ntemp]), "temperaturemin":np.array([tmin]),
            "temperaturemax":np.array([tmax]), "sigma":sigmav.T}
        return ("asigma_loc", data)
