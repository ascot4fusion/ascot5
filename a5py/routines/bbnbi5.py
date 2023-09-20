"""
Generate a generic neutral beam injector geometry for BBNBI5.

Jari Varje <jari.varje@tokamakenergy.co.uk>
"""
import numpy as np
import scipy.constants as const
import a5py.ascot5io.nbi as nbi

def generate(fn, Rgrid, phigrid, zgrid, tanrad, focal_length, dgrid, nbeamlet, anum, znum, mass, energy, efrac, power, div, tilt=None, tanz=None, div_halo_frac=None, div_halo=None, desc=None):
    """
    Generate generic injector. The injector grid can be either circular or
    rectangular, and will be randomly populated with beamlets with a common
    focus point.

    Args:
        fn : str
            Full path to the HDF5 file.
        Rgrid, phigrid, zgrid : float
            Position of the center point of the injector grid [m], [deg], [m]
        tanrad : float
            Tangency radius of the injector centerline
        focal_length : float
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
            Full injection energy [keV]
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
        div_halo_frac : float, optional
            Fraction of particles with halo divergence
        div_halo : float, optional
            Horizontal and vertical divergences [rad]
        desc : str, optional
            Input description
    """

    phigrid = phigrid * np.pi/180.0
    phidir = np.pi + phigrid - np.arcsin(tanrad/Rgrid);

    grid = np.array([Rgrid * np.cos(phigrid), Rgrid * np.sin(phigrid), zgrid])

    unitd = np.array([np.cos(phidir), np.sin(phidir), 0])
    unitz = np.array([0, 0, 1])
    unith = np.cross(unitd, unitz)

    if tilt is not None:
        tilt = tilt * np.pi/180.0
    elif tanz is not None:
        tilt = tanz / Rgrid*np.cos(phidir)
    else:
        tilt = 0

    focus = grid + focal_length*unitd + focal_length*np.tan(tilt)*unitz

    n=dict()
    n["id"] = 1
    n["nbeamlet"] = nbeamlet
    n["beamletx"] = np.zeros(nbeamlet)
    n["beamlety"] = np.zeros(nbeamlet)
    n["beamletz"] = np.zeros(nbeamlet)
    n["beamletdx"] = np.zeros(nbeamlet)
    n["beamletdy"] = np.zeros(nbeamlet)
    n["beamletdz"] = np.zeros(nbeamlet)
    n["div_h"] = div[0]
    n["div_v"] = div[1]

    if div_halo_frac is not None:
        n["div_halo_frac"] = div_halo_frac
        n["div_halo_h"] = div_halo[0]
        n["div_halo_v"] = div_halo[1]
    else:
        n["div_halo_frac"] = 0.0
        n["div_halo_h"] = 1e-10
        n["div_halo_v"] = 1e-10

    n["anum"] = anum
    n["znum"] = znum
    n["mass"] = mass*const.u
    n["energy"] = energy*const.e
    n["efrac"] = np.array(efrac)
    n["power"] = power

    for i in range(nbeamlet):
        if np.array(dgrid).size == 1:    # circular grid
            d = np.inf
            while d > dgrid/2:
                dx = dgrid/2 * (2*np.random.rand(1)-1)
                dy = dgrid/2 * (2*np.random.rand(1)-1)
                d = np.sqrt(dx**2+dy**2)
        else:    # rectangular grid
            dx = dgrid[0]/2 * (2*np.random.rand(1)-1)
            dy = dgrid[1]/2 * (2*np.random.rand(1)-1)

        beamletxyz = grid + dx*unith + dy*unitz
        dirxyz = focus - beamletxyz
        dirxyz = dirxyz / np.linalg.norm(dirxyz)

        n["beamletx"][i] = beamletxyz[0]
        n["beamlety"][i] = beamletxyz[1]
        n["beamletz"][i] = beamletxyz[2]
        n["beamletdx"][i] = dirxyz[0]
        n["beamletdy"][i] = dirxyz[1]
        n["beamletdz"][i] = dirxyz[2]

    nbi.write_hdf5(fn, [n], desc)
