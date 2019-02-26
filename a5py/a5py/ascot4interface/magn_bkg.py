import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator as rgi

def read_magn_bkg(fn,hdrfn):
    str = dict()
    with open(fn) as fh:
        tmp = list(map(float,fh.readline().split()))
        str['phi0']       = tmp[0]
        str['nSector']    = int(tmp[1])
        str['nPhi']       = int(tmp[2])
        str['nCoil']      = int(tmp[3])
        str['zeroAtCoil'] = int(tmp[4])

        r1,r2,nr = [float(number) for number in fh.readline().split()]
        nr = int(nr)
        z1,z2,nz = [float(number) for number in fh.readline().split()]
        nz = int(nz)

        str['r'] = np.linspace(r1,r2,nr)
        str['z'] = np.linspace(z1,z2,nz)
        if(str['nPhi'] > 1):
            dphi = 360.0 / ( str['nSector'] * str['nPhi'] )
            str['phi'] = ( np.linspace( str['phi0'] + dphi * 0.5,
                                     str['phi0'] + 360.0 / str['nSector'] - dphi * 0.5,
                                     str['nPhi'] ) )

        str['phimap_tor'] = np.array([int(number) for number in fh.readline().split()])
        str['phimap_pol'] = np.array([int(number) for number in fh.readline().split()])

        data = np.array(fh.read().split(), dtype=float).flatten()

        sz2d = nr*nz
        sz3d = nr*nz*str['nPhi']
        str['psi']  = data[:sz2d].reshape(nz,nr)
        print(np.shape(data),nz,str['nPhi'],nr,sz2d,sz3d)
        str['br']   = data[sz2d+0*sz3d:sz2d+1*sz3d].reshape(nz,str['nPhi'],nr).squeeze()
        str['bphi'] = data[sz2d+1*sz3d:sz2d+2*sz3d].reshape(nz,str['nPhi'],nr).squeeze()
        str['bz']   = data[sz2d+2*sz3d:sz2d+3*sz3d].reshape(nz,str['nPhi'],nr).squeeze()

    read_magn_header(hdrfn,str)
    return str

def read_magn_header(fn,str):
    with open(fn) as fh:
        # first four lines contain nothing interesting
        fh.readline()
        fh.readline()
        fh.readline()
        fh.readline()

        # Next three lines contain axis psi, R, and z values
        tmp = [float(number) for number in fh.readline().split()]
        str['psi0']       = tmp[0]
        str['psi1']       = tmp[1]

        tmp = [float(number) for number in fh.readline().split()]
        str['axis_r']       = tmp[0]

        tmp = [float(number) for number in fh.readline().split()]
        str['axis_z']       = tmp[0]

    return str

def read_magn_bkg_stellarator(fn):

    with h5py.File(fn, 'r') as f: # Open for reading
        str = dict()

        str['r'] = f['bfield/stellarator/r'][:]
        str['phi'] = f['bfield/stellarator/phi'][:]
        str['z'] = f['bfield/stellarator/z'][:]

        str['br'] = f['bfield/stellarator/br'][:]
        str['bphi'] = f['bfield/stellarator/bphi'][:]
        str['bz'] = f['bfield/stellarator/bz'][:]
        str['s'] = f['bfield/stellarator/s'][:]

        str['axis_r'] = f['bfield/stellarator/axis_R'][:]
        str['axis_phi'] = f['bfield/stellarator/axis_phi'][:]
        str['axis_z'] = f['bfield/stellarator/axis_z'][:]

        str['n_periods'] = f['bfield/stellarator/toroidalPeriods'][:]
        try:
            str['symmetrymode'] = f['bfield/stellarator/symmetrymode'][:]
        except KeyError:
            print("Warning! No symmetry mode specified in input.h5/bfield")
            print("Defaulting to stellarator symmetry")
            str['symmetrymode'] = 0

    return str

def stellarator_bfield_sector2full(data):
    out = data
    # Data is in the format (r, phi, z)
    out['r'] = data['r']
    out['z'] = data['z']
    out['phi'] = np.concatenate((data['phi'][:-1], data['phi'] + data['phi'][-1]))
    # br
    br = data['br'][:, :-1, :]
    br_sym = -data['br'][:, -1::-1, -1::-1]
    out['br'] = np.concatenate((br, br_sym),1)
    # bphi
    bphi = data['bphi'][:, :-1, :]
    bphi_sym = data['bphi'][:, -1::-1, -1::-1]
    out['bphi'] = np.concatenate((bphi, bphi_sym),1)
    # bz
    bz = data['bz'][:, :-1, :]
    bz_sym = data['bz'][:, -1::-1, -1::-1]
    out['bz'] = np.concatenate((bz, bz_sym),1)
    # s
    s = data['s'][:, :-1, :]
    s_sym = data['s'][:, -1::-1, -1::-1]
    out['s'] = np.concatenate((s, s_sym),1)
    out['symmetrymode'] == 1

    return out

def stellarator_psi_lims(data):
    ndense = 3;
    rvec   = np.linspace(data['r'][0], data['r'][-1],
                         ndense * data['r'].size)
    phivec = np.linspace(data['phi'][0], data['phi'][-1],
                         ndense * data['phi'].size)
    zvec   = np.linspace(data['z'][0], data['z'][-1],
                         ndense * data['z'].size)
    psi_int = rgi((data['r'].flatten(), data['phi'].flatten(),
                   data['z'].flatten()),
                  data['s'])

    min_psi = np.amin(data['s'])
    max_psi = np.amax(data['s'])
    rgrid, phigrid, zgrid = np.meshgrid(rvec, phivec, zvec)
    for r in rvec:
        for phi in phivec:
            psi = psi_int(np.array([r * np.ones(zvec.shape),
                                    phi * np.ones(zvec.shape),
                                    zvec]).T)
            min_psi = np.amin([min_psi, np.amin(psi)])
            max_psi = np.amax([max_psi, np.amax(psi)])
    return min_psi, max_psi
