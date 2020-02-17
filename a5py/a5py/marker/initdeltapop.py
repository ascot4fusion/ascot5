"""
Initializes delta population of markers at given intervals.
"""

import numpy as np

import scipy.constants as constants
import a5py.ascot5io.mrk_gc as mrkmod
from a5py.ascotpy import Ascotpy


def initmarkers(fn, n, mass, charge, r, energy, pitch, weight=1, time=0,
                desc=None, separateruns=False):

    # Find the OMP z value
    a5 = Ascotpy(fn)
    a5.init(bfield=True)
    z0 = a5.evaluate(100, 0, 0, 0, "axis")["axisz"]
    a5.free(bfield=True)

    # Generate radial coordinates
    def gencoords(x):
        if isinstance(x, list):
            if x[-1] == 0:
                x  = x[0] + (x[1]-x[0]) * np.random.rand(n,)
                nx = 1
            else:
                nx = x[2]
                x  = np.linspace(x[0], x[1], x[2])
        else:
            x  = x*np.ones((n,))
            nx = 1

        return (x,nx)

    r, nr      = gencoords(r)
    energy, ne = gencoords(energy)
    pitch, nx  = gencoords(pitch)

    rcoords = np.array([])
    ecoords = np.array([])
    pcoords = np.array([])
    for ir in range(nr):
        for ie in range(ne):
            for ip in range(nx):

                if nr == 1:
                    rcoords = np.append(rcoords, r)
                else:
                    rcoords = np.append(rcoords, r[ir]*np.ones((n,)))

                if ne == 1:
                    ecoords = np.append(ecoords, energy)
                else:
                    ecoords = np.append(ecoords, energy[ie]*np.ones((n,)))

                if nx == 1:
                    pcoords = np.append(pcoords, pitch)
                else:
                    pcoords = np.append(pcoords, pitch[ip]*np.ones((n,)))

    r      = rcoords
    energy = ecoords
    pitch  = pcoords

    N = rcoords.size
    ids = np.arange(N)+1

    if separateruns:
        pass

    else:
        mix = np.random.permutation(N)
        ids    = ids[mix]
        r      = r[mix]
        energy = energy[mix]
        pitch  = pitch[mix]

        mrkmod.write_hdf5(
            fn     = fn,
            n      = N,
            ids    = ids,
            mass   = mass*np.ones((N,1)),
            charge = charge*np.ones((N,1)),
            r      = r,
            phi    = 360*np.random.rand(N, 1),
            z      = z0*np.ones((N,1)),
            energy = energy,
            pitch  = pitch,
            zeta   = 2*np.pi*np.random.rand(N, 1),
            anum   = np.ones((N,1)),
            znum   = np.ones((N,1)),
            weight = weight*np.ones((N,1)),
            time   = time*np.ones((N,1)),
            desc   = desc)
