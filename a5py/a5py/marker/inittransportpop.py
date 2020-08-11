"""
Initialize markers for transport coefficient evaluations and checks.

Markers are initialized at the outer mid-plane. The point is to divide the
phase space in to a grid, and initialize a number of markers at each grid point.
This way one obtains the transport coefficients as a function phase-space at
these grid nodes.
"""

import numpy as np

import a5py.ascot5io.mrk_gc as mrk_gc
import a5py.ascot5io.mrk_fl as mrk_fl

from a5py.ascotpy import Ascotpy


def init_rho(fn, n, rhogrid, rgrid=None, randomize_rho=False, bfield=True,
             desc=None):
    """
    Initialize field line markers in radius.

    Markers are initialized so that at every rho coordinate there are n markers.
    The coordinates are given in rhogrid array at any order. If randomize_rho is
    True, n markers are randomly and uniformly distributed between min and max
    rhogrid values.

    Use rgrid to give R values explicitly in which case rhogrid is ignored.
    """
    rhogrid = np.atleast_1d(rhogrid)

    def getvals(grid):
        if randomize_rho:
            vals = np.amin(grid) + ( np.amax(grid)- np.amin(grid) ) \
                      * np.random.rand(n)
        else:
            vals = np.zeros((grid.size*n,))
            for i in range(grid.size):
                vals[i*n:(i+1)*n] = grid[i]

        return vals

    if rgrid is None:
        rhovals = getvals(rhogrid)
    else:
        rhovals = getvals(rgrid)

    # Find OMP R,z values
    a5 = Ascotpy(fn)
    a5.init(bfield=bfield)
    if rgrid is None:
        rz_omp = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    else:
        rz_omp = a5.get_rhotheta_rz( 0.5*np.ones(rhovals.shape), 0, 0, 0 )
        rz_omp[0] = rhovals
    a5.free(bfield=bfield)

    mrk = {}
    mrk["n"]      = rz_omp[0].size
    mrk["ids"]    = np.arange(1, mrk["n"]+1)
    mrk["r"]      = rz_omp[0]
    mrk["z"]      = rz_omp[1]
    mrk["phi"]    = 360*np.random.rand(mrk["n"])
    mrk["pitch"]  = 1 * np.ones((mrk["n"],))
    mrk["time"]   = 0 * np.ones((mrk["n"],))
    mrk["weight"] = 1 * np.ones((mrk["n"],))

    mix = np.random.permutation(mrk["n"])
    mrk["ids"] = mrk["ids"][mix]
    mrk["r"]   = mrk["r"][mix]

    mrk_fl.write_hdf5(fn, **mrk, desc=desc)


def init_rhoenergypitch(fn, n, mass, charge, anum, znum, rhogrid, energygrid,
                        pitchgrid, rgrid=None,
                        randomize_rho=False, randomize_energy=False,
                        randomize_pitch=False, time=0, desc=None):
    """
    Initialize guiding center markers in radius, energy, and pitch.

    See init_rho for how markers are distributed. This function distributes
    markers also in energy and pitch. At each grid point n markers are
    initialized.

    Use rgrid to give R values explicitly in which case rhogrid is ignored.
    """

    rhogrid    = np.atleast_1d(rhogrid)
    energygrid = np.atleast_1d(energygrid)
    pitchgrid  = np.atleast_1d(pitchgrid)

    ntotal = n
    if not randomize_rho:
        ntotal *= np.unique(rhogrid).size
    if not randomize_energy:
        ntotal *= np.unique(energygrid).size
    if not randomize_pitch:
        ntotal *= np.unique(pitchgrid).size

    rhovals    = np.zeros((ntotal,))
    energyvals = np.zeros((ntotal,))
    pitchvals  = np.zeros((ntotal,))
    for irho in range(rhogrid.size):
        for ienergy in range(energygrid.size):
            for ipitch in range(pitchgrid.size):
                idx =   ipitch * (energygrid.size * rhogrid.size) \
                      + ienergy * rhogrid.size + irho
                if randomize_rho and rgrid is None:
                    rhovals[idx*n:(idx+1)*n] = \
                    np.amin(rhogrid) + ( np.amax(rhogrid)- np.amin(rhogrid) ) \
                        * np.random.rand(n)
                elif rgrid is None:
                    rhovals[idx*n:(idx+1)*n] = rhogrid[irho]

                if randomize_rho and rgrid is not None:
                    rhovals[idx*n:(idx+1)*n] = \
                    np.amin(rgrid) + ( np.amax(rgrid)- np.amin(rgrid) ) \
                        * np.random.rand(n)
                elif rgrid is not None:
                    rhovals[idx*n:(idx+1)*n] = rgrid[irho]

                if randomize_energy:
                    energyvals[idx*n:(idx+1)*n] = \
                      np.amin(energygrid) \
                    + ( np.amax(energygrid)- np.amin(energygrid) ) \
                    * np.random.rand(n)
                else:
                    energyvals[idx*n:(idx+1)*n] = energygrid[ienergy]

                if randomize_pitch:
                    pitchvals[idx*n:(idx+1)*n] = \
                      np.amin(pitchgrid)
                    + ( np.amax(pitchgrid)- np.amin(pitchgrid) ) \
                    * np.random.rand(n)
                else:
                    pitchvals[idx*n:(idx+1)*n] = pitchgrid[ipitch]

    # Find OMP R,z values
    a5 = Ascotpy(fn)
    a5.init(bfield=True)
    if rgrid is None:
        rz_omp = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    else:
        rz_omp = a5.get_rhotheta_rz( 0.5*np.ones(rhovals.shape), 0, 0, 0 )
        rz_omp[0] = rhovals
    a5.free(bfield=True)

    mrk = {}
    mrk["n"]      = rz_omp[0].size
    mrk["ids"]    = np.arange(1, mrk["n"]+1)
    mrk["r"]      = rz_omp[0]
    mrk["z"]      = rz_omp[1]
    mrk["phi"]    = 360*np.random.rand(mrk["n"])
    mrk["zeta"]   = 2*np.pi*np.random.rand(mrk["n"])
    mrk["pitch"]  = pitchvals
    mrk["energy"] = energyvals
    mrk["time"]   = 0 * np.ones((mrk["n"],))
    mrk["weight"] = 1 * np.ones((mrk["n"],))

    mix = np.random.permutation(mrk["n"])
    mrk["ids"]    = mrk["ids"][mix]
    mrk["r"]      = mrk["r"][mix]
    mrk["pitch"]  = mrk["pitch"][mix]
    mrk["energy"] = mrk["energy"][mix]

    mrk["mass"]   = mass.to("amu").d   * np.ones((mrk["n"],))
    mrk["charge"] = charge * np.ones((mrk["n"],))
    mrk["anum"]   = anum   * np.ones((mrk["n"],))
    mrk["znum"]   = znum   * np.ones((mrk["n"],))

    mrk_gc.write_hdf5(fn=fn, **mrk, desc=desc)


def init_rhopparapperp(fn, n, mass, charge, anum, znum, rhogrid, pparagrid,
                       pperpgrid, rgrid=None,
                       randomize_rho=False, randomize_ppara=False,
                       randomize_pperp=False, time=0, desc=None):
    """
    Initialize guiding center markers in radius, ppara, and pperp.

    See init_rho for how markers are distributed. This function distributes
    markers also in ppara and pperp. At each grid point n markers are
    initialized.

    Use rgrid to give R values explicitly in which case rhogrid is ignored.
    """

    rhogrid   = np.atleast_1d(rhogrid)
    pparagrid = np.atleast_1d(pparagrid)
    pperpgrid = np.atleast_1d(pperpgrid)

    ntotal = n
    if not randomize_rho:
        ntotal *= np.unique(rhogrid).size
    if not randomize_ppara:
        ntotal *= np.unique(pparagrid).size
    if not randomize_pperp:
        ntotal *= np.unique(pperpgrid).size

    rhovals   = np.zeros((ntotal,))
    pparavals = np.zeros((ntotal,))
    pperpvals = np.zeros((ntotal,))
    for irho in range(rhogrid.size):
        for ippara in range(pparagrid.size):
            for ipperp in range(pperpgrid.size):
                idx =   ipperp * (pparagrid.size + rhogrid.size) \
                      + ippara * rhogrid.size + irho
                if randomize_rho:
                    rhovals[idx*n:(idx+1)*n] = \
                    np.amin(rhogrid) + ( np.amax(rhogrid)- np.amin(rhogrid) ) \
                        * np.random.rand(n)
                else:
                    rhovals[idx*n:(idx+1)*n] = rhogrid[irho]

                if randomize_ppara:
                    pparavals[idx*n:(idx+1)*n] = \
                      np.amin(pparagrid) \
                    + ( np.amax(pparagrid)- np.amin(pparagrid) ) \
                    * np.random.rand(n)
                else:
                    pparavals[idx*n:(idx+1)*n] = pparagrid[ippara]

                if randomize_pperp:
                    pitchvals[idx*n:(idx+1)*n] = \
                      np.amin(pperpgrid)
                    + ( np.amax(pperpgrid)- np.amin(pperpgrid) ) \
                    * np.random.rand(n)
                else:
                    pperpvals[idx*n:(idx+1)*n] = pperpgrid[ipperp]

    # Find OMP R,z values
    a5 = Ascotpy(fn)
    a5.init(bfield=True)
    if rgrid is None:
        rz_omp = a5.get_rhotheta_rz( rhovals, 0, 0, 0 )
    else:
        rz_omp = a5.get_rhotheta_rz( 0.5*np.ones(rhovals.shape), 0, 0, 0 )
        rz_omp[0] = rhovals
    a5.free(bfield=True)

    mrk = {}
    mrk["n"]      = rz_omp[0].size
    mrk["ids"]    = np.arange(1, mrk["n"]+1)
    mrk["r"]      = rz_omp[0]
    mrk["z"]      = rz_omp[1]
    mrk["phi"]    = 360*np.random.rand(mrk["n"])
    mrk["zeta"]   = 2*np.pi*np.random.rand(mrk["n"])
    mrk["pitch"]  = pparavals / np.sqrt(pparavals**2 + pperpvals**2)
    mrk["time"]   = 0 * np.ones((mrk["n"],))
    mrk["weight"] = 1 * np.ones((mrk["n"],))

    pnorm2 =  pparavals**2 + pperpvals**2
    mrk["energy"] = np.sqrt(restmass**2 + pnorm2 * c**2) - restmass
    mrk["mass"]   = mass   * np.ones((mrk["n"],))
    mrk["charge"] = charge * np.ones((mrk["n"],))
    mrk["anum"]   = anum   * np.ones((mrk["n"],))
    mrk["znum"]   = znum   * np.ones((mrk["n"],))

    mix = np.random.permutation(mrk["n"])
    mrk["ids"]    = mrk["ids"][mix]
    mrk["r"]      = mrk["r"][mix]
    mrk["pitch"]  = mrk["pitch"][mix]
    mrk["energy"] = mrk["energy"][mix]

    mrk_gc.write_hdf5(fn=fn, **mrk, desc=desc)
