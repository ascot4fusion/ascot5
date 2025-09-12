"""Templates that use analytical models to create inputs without any
actual data."""
import numpy as np
import unyt

from scipy.interpolate import interpn
from scipy.integrate   import quad

from a5py.physlib import aeq, species2properties

from .template import InputTemplate


class PremadeMagneticField(InputTemplate):

    def __init__(self, ascot, field, axisymmetric=True, splines=False):
        """Create one of the premade magnetic fields.

        Parameters
        ----------
        ascot : :class:`Ascot`
            The Ascot object where the input will be created.
        field : str
            Name of the premade magnetic field to create.

            The options are:
            - `"iter-circular"`: Circular ITER-like field.
            - `"iter-baseline"`: more realistic ITER-like field with elongation
              and similar parameters as in the baseline scenario.
            - `"nstx-doublenull"`: Double-null configuration of the NSTX field.
        axisymmetric : bool, optional
            If True, the field contains toroidal ripple.

            Unfortunately, the model for the TF ripple is not divergence free.
        splines : bool, optional
            If True, the field is interpolated with splines in the simulation.

            By default, the field is analytical.
        """
        data = {}
        if field == "iter-circular":
            c = aeq.parameters2coefficients(
                A=1.0, epsilon=0.32, kappa=1.0, delta=0.0,
                xpoint=None, symmetric=True,
                )
            data.update({
                "rmajor":6.2 * unyt.m, "axisz":0.0 * unyt.m,
                "axisb":5.3 * unyt.T, "psiscaling":200.0,
                "coefficients":np.append(c, 1.0),
                })
            if splines:
                data.update({
                    "rgrid":np.linspace(4.0, 8.5, 120) * unyt.m,
                    "zgrid":np.linspace(-4.0, 4.0, 200) * unyt.m,
                    })
            if not axisymmetric:
                data.update({
                    "nripple":18, "rminor":2.0 * unyt.m,
                    "ripplepenetration":0.2, "ripplescaling":0.5,
                    })
                if splines:
                    data["phigrid"] = np.linspace(0.0, 360, 180) * unyt.deg
        elif field == "iter-baseline":
            c = aeq.parameters2coefficients(
                A=-0.155, epsilon=0.32, kappa=1.7, delta=0.33,
                xpoint=(0.88, -0.6), symmetric=False,
                )
            data.update({
                "rmajor":6.2 * unyt.m, "axisz":0.0 * unyt.m,
                "axisb":5.3 * unyt.T, "psiscaling":200.0,
                "coefficients":np.append(c, -0.155),
                })
            if splines:
                data.update({
                    "rgrid":np.linspace(4.0, 8.5, 120) * unyt.m,
                    "zgrid":np.linspace(-4.0, 4.0, 200) * unyt.m,
                    })
            if not axisymmetric:
                data.update({
                    "nripple":18, "rminor":2.0 * unyt.m,
                    "ripplepenetration":0.2, "ripplescaling":0.5,
                    })
                if splines:
                    data["phigrid"] = np.linspace(0.0, 360, 180) * unyt.deg
        elif field == "nstx-doublenull":
            # The X-point is from the "One-size-fits-all" paper:
            # (1-1.1*delta*eps, 1.1*kappa*eps)
            c = aeq.parameters2coefficients(
                A=0.0, epsilon=0.78, kappa=2.0, delta=0.35,
                xpoint=(1.0 - 1.1*0.35*0.78, 1.1*2.0*0.78), symmetric=True,
                )
            data.update({
                "rmajor":0.85 * unyt.m, "axisz":0.0 * unyt.m,
                "axisb":0.3 * unyt.T, "psiscaling":200.0,
                "coefficients":np.append(c, 0.0),
                })
            if splines:
                data.update({
                    "rgrid":np.linspace(0.1, 1.7, 120) * unyt.m,
                    "zgrid":np.linspace(-2.0, 2.0, 200) * unyt.m,
                    })
            if not axisymmetric:
                data.update({
                    "nripple":12, "rminor":0.68 * unyt.m,
                    "ripplepenetration":0.2, "ripplescaling":0.5,
                    })
                if splines:
                    data["phigrid"] = np.linspace(0.0, 360, 180) * unyt.deg
        else:
            raise ValueError(f"No premade field found: {field}")

        if splines and axisymmetric:
            input_type = "bfield2d"
        elif splines and not axisymmetric:
            input_type = "bfield3d"
        else:
            input_type = "bfieldanalytical"
        super().__init__(ascot, input_type, data)


def add_step_like_ripple(self, b2d=None, ncoil=12, rcoil=8.0,
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


class FlatPlasma(InputTemplate):

    def __init__(
            self,
            ascot,
            density=10e20*unyt.m**-3,
            temperature=10e3*unyt.eV,
            rotation=0.0*unyt.m/unyt.s,
            species="H"
            ):
        """Create flat plasma density that is constant inside the separatrix
        and (almost) zero elsewhere.

        Parameters
        ----------
        ascot : :class:`Ascot`
            The Ascot object where the input will be created.
        density : float, optional
            Plasma density inside the separatrix.
        temperature : float, optional
            Plasma temperature (constant everywhere).
        rotation : float, optional
            Plasma rotation (constant everywhere).
        species : str or [str], optional
            Plasma ion species.

            Density of each species is chosen so that the charge density for
            each species is equal.
        """
        if not isinstance(species, list) or not isinstance(species, tuple):
            species = [species]

        qtot = np.sum([species2properties(s).charge.v for s in species])
        iondensity = np.ones((4,len(species))) * density
        iondensity[2:,:] = 1.0
        iondensity /= qtot
        data = {
            "rhogrid":np.array([0.0, 1.0, 1.05, 10.0]),
            "species":species,
            "iondensity":iondensity,
            "iontemperature":np.array([temperature]*4),
            "electrondensity":np.array([density, density, 1.0, 1.0]),
            "electrontemperature":np.array([temperature]*4),
            "rotation":np.array([rotation]*4),
            }
        super().__init__(ascot, "plasma1d", data)


def create_flat_neutral_profile(self, density=10e20, temperature=10e3, anum=1, znum=1):
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
            -5.4914200 * ( X ),
            -3.4294800 * ( 2*X**2 - 1 ),
            -1.9837700 * ( 4*X**3 - 3*X ),
            -0.8780090 * ( 8*X**4 - 8*X**2 + 1 ),
            -0.1989320 * ( 16*X**5 - 20*X**3 + 5*X ),
            +0.0837431 * ( 32*X**6 - 48*X**4 + 18*X**2 - 1 ),
            +0.1212520 * ( 64*X**7 - 112*X**5 + 56*X**3 - 7*X ),
            +0.0827182 * ( 128*X**8 - 256*X**6 + 160*X**4 - 32*X**2 + 1 )
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


def generic_nbi(r, phi, z, tanrad, focallen, dgrid, nbeamlet,
                anum, znum, mass, energy, efrac, power, div, tilt=None,
                tanz=None, divhalofrac=None, divhalo=None, ids=1):
    """Generate a generic injector.

    The injector grid can be either circular or rectangular, and will be
    randomly populated with beamlets with a common focus point.

    The contents of this function can be used as a template for implementing
    machine-specific injectors with an actual geometry.

    Parameters
    ----------
    r : float
        Injector center point R-coordinate [m].
    phi : float
        Injector center point phi-coordinate [deg].
    z : float
        Injector center point z-coordinate [m].
    tanrad : float
        Tangency radius of the injector centerline [m].
    focallen : float
        Distance to the focus point from the center of the grid
    dgrid : float or array_like (2,)
        Diameter of circular grid or width and height of rectangular grid
    nbeamlet : int
        Number of beamlets in this injector
    anum : int
        Mass number of injected species
    znum : int
        Nuclear charge number of injected species
    mass : float
        Mass of the injected species [amu]
    energy : float
        Full injection energy [eV]
    efrac : array_like (3,)
        Particle fractions for full, 1/2 and 1/3 energies
    power : float
        Injected power [W]
    div : array_like (2,)
        Horizontal and vertical divergences [rad]
    tilt : float, optional
        Vertical tilt angle of the beam centerline [deg]
    tanz : float, optional
        Vertical shift of the tangency radius point [m]
    divhalofrac : float, optional
        Fraction of particles with halo divergence
    divhalo : float, optional
        Horizontal and vertical divergences [rad]
    ids : int, optional
        Unique numerical identifier for this injector.

    Returns
    -------
    inj : :class:`Injector`
        A single injector.
    """
    phi    = phi * np.pi/180.0
    phidir = np.pi + phi - np.arcsin(tanrad/r);
    center = np.array([r * np.cos(phi), r * np.sin(phi), z])

    unitd = np.array([np.cos(phidir), np.sin(phidir), 0.0])
    unitz = np.array([0.0, 0.0, 1.0])
    unith = np.cross(unitd, unitz)

    if tilt is not None:
        tilt = tilt * np.pi/180.0
    elif tanz is not None:
        tilt = tanz / r * np.cos(phidir)
    else:
        tilt = 0

    focus = center + focallen * unitd + focallen * np.tan(tilt) * unitz

    n = {}
    n["ids"]       = ids
    n["nbeamlet"]  = nbeamlet
    n["beamletx"]  = np.zeros(nbeamlet)
    n["beamlety"]  = np.zeros(nbeamlet)
    n["beamletz"]  = np.zeros(nbeamlet)
    n["beamletdx"] = np.zeros(nbeamlet)
    n["beamletdy"] = np.zeros(nbeamlet)
    n["beamletdz"] = np.zeros(nbeamlet)
    n["divh"]      = div[0]
    n["divv"]      = div[1]

    if divhalofrac is not None:
        n["divhalofrac"] = divhalofrac
        n["divhaloh"]    = divhalo[0]
        n["divhalov"]    = divhalo[1]
    else:
        n["divhalofrac"] = 0.0
        n["divhaloh"]    = 1e-10
        n["divhalov"]    = 1e-10

    # Injected species
    n["anum"]   = anum
    n["znum"]   = znum
    n["mass"]   = mass
    n["energy"] = energy
    n["efrac"]  = np.array(efrac)
    n["power"]  = power

    for i in range(nbeamlet):
        if np.array(dgrid).size == 1:
            # circular grid
            d = np.inf
            while d > 0.5 * dgrid:
                dx = 0.5 * dgrid * ( 2 * np.random.rand(1) - 1 )
                dy = 0.5 * dgrid * ( 2 * np.random.rand(1) - 1 )
                d  = np.sqrt(dx**2 + dy**2)
        else:
            # rectangular grid
            dx = 0.5 * dgrid[0] * ( 2 * np.random.rand(1) - 1 )
            dy = 0.5 * dgrid[1] * ( 2 * np.random.rand(1) - 1 )

        beamletxyz = center + dx*unith + dy*unitz
        dirxyz = focus - beamletxyz
        dirxyz = dirxyz / np.linalg.norm(dirxyz)

        n["beamletx"][i]  = beamletxyz[0]
        n["beamlety"][i]  = beamletxyz[1]
        n["beamletz"][i]  = beamletxyz[2]
        n["beamletdx"][i] = dirxyz[0]
        n["beamletdy"][i] = dirxyz[1]
        n["beamletdz"][i] = dirxyz[2]

    return Injector(**n)
