import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator as rgi
from a5py.ascotpy import Ascotpy
from scipy.optimize import fmin

def read_magn_bkg(fn,hdrfn):
    data = dict()
    with open(fn) as fh:
        tmp = list(map(float,fh.readline().split()))
        data['phi0']       = tmp[0]
        data['nSector']    = int(tmp[1])
        data['nPhi']       = int(tmp[2])
        data['nCoil']      = int(tmp[3])
        data['zeroAtCoil'] = int(tmp[4])

        r1,r2,nr = [float(number) for number in fh.readline().split()]
        nr = int(nr)
        z1,z2,nz = [float(number) for number in fh.readline().split()]
        nz = int(nz)

        data['r'] = np.linspace(r1,r2,nr)
        data['z'] = np.linspace(z1,z2,nz)
        if(data['nPhi'] > 1):
            dphi = 360.0 / ( data['nSector'] * data['nPhi'] )
            data['phi'] = ( np.linspace( data['phi0'] + dphi * 0.5,
                                     data['phi0'] + 360.0 / data['nSector'] - dphi * 0.5,
                                     data['nPhi'] ) )

        data['phimap_tor'] = np.array([int(number) for number in fh.readline().split()])
        data['phimap_pol'] = np.array([int(number) for number in fh.readline().split()])

        data = np.array(fh.read().split(), dtype=float).flatten()

        sz2d = nr*nz
        sz3d = nr*nz*data['nPhi']
        data['psi']  = data[:sz2d].reshape(nz,nr)
        print(np.shape(data),nz,data['nPhi'],nr,sz2d,sz3d)
        data['br']   = data[sz2d+0*sz3d:sz2d+1*sz3d].reshape(nz,data['nPhi'],nr).squeeze()
        data['bphi'] = data[sz2d+1*sz3d:sz2d+2*sz3d].reshape(nz,data['nPhi'],nr).squeeze()
        data['bz']   = data[sz2d+2*sz3d:sz2d+3*sz3d].reshape(nz,data['nPhi'],nr).squeeze()

        data['psi']  = np.transpose(data['psi'])
        data['br']   = np.transpose(data['br'])
        data['bphi'] = np.transpose(data['bphi'])
        data['bz']   = np.transpose(data['bz'])

    read_magn_header(hdrfn,data)
    return data

def read_magn_header(fn,data):
    with open(fn) as fh:
        # first four lines contain nothing interesting
        fh.readline()
        fh.readline()
        fh.readline()
        fh.readline()

        # Next three lines contain axis psi, R, and z values
        tmp = [float(number) for number in fh.readline().split()]
        data['psi0']       = tmp[0]
        data['psi1']       = tmp[1]

        tmp = [float(number) for number in fh.readline().split()]
        data['axis_r']       = tmp[0]

        tmp = [float(number) for number in fh.readline().split()]
        data['axis_z']       = tmp[0]

    return data

def read_magn_bkg_stellarator(fn):

    with h5py.File(fn, 'r') as f: # Open for reading
        data = dict()

        data['r'] = f['bfield/stellarator/r'][:]
        data['phi'] = f['bfield/stellarator/phi'][:]
        data['z'] = f['bfield/stellarator/z'][:]

        data['br'] = f['bfield/stellarator/br'][:]
        data['bphi'] = f['bfield/stellarator/bphi'][:]
        data['bz'] = f['bfield/stellarator/bz'][:]
        data['s'] = f['bfield/stellarator/s'][:]

        data['axis_r'] = f['bfield/stellarator/axis_R'][:]
        data['axis_phi'] = f['bfield/stellarator/axis_phi'][:]
        data['axis_z'] = f['bfield/stellarator/axis_z'][:]

        data['n_periods'] = f['bfield/stellarator/toroidalPeriods'][:]
        try:
            data['symmetrymode'] = f['bfield/stellarator/symmetrymode'][:]
        except KeyError:
            print("Warning! No symmetry mode specified in input.h5/bfield")
            print("Defaulting to stellarator symmetry")
            data['symmetrymode'] = 0
        if (data['phi'][0] == np.mod(data['phi'][-1],360/data['n_periods'])):
            print("Warning! Removing duplicate bfield data point.")
            data = bfield_remove_duplicate_phi(data)
        if(data['symmetrymode'] == 0):
            print("Converting stellarator symmetric input to periodic.")
            data = stellarator_bfield_sector2full(data)
        if (data['axis_phi'][0] == np.mod(data['axis_phi'][-1],360)):
            print("Warning! Removing duplicated axis datapoint.")
            data['axis_r'] = data['axis_r'][0:-1]
            data['axis_phi'] = data['axis_phi'][0:-1]
            data['axis_z'] = data['axis_z'][0:-1]
        # Transpose to ascot5io format
        data['br']   = np.transpose(data['br'],   (2,1,0))
        data['bphi'] = np.transpose(data['bphi'], (2,1,0))
        data['bz']   = np.transpose(data['bz'],   (2,1,0))
        data['s']    = np.transpose(data['s'],    (2,1,0))

    return data

def stellarator_bfield_sector2full(data):
    out = data
    # Data is in the format f(z, phi, r)
    out['r'] = data['r']
    out['z'] = data['z']
    out['phi'] = np.concatenate( (data['phi'][:-1],
                                  data['phi'][:] + data['phi'][-1]) )
    # br
    br = data['br'][:, :-1, :]
    br_sym = - data['br'][-1::-1, -1:0:-1, :]
    out['br'] = np.concatenate((br, br_sym),1)
    # bphi
    bphi = data['bphi'][:, :-1, :]
    bphi_sym = data['bphi'][-1::-1, -1:0:-1, :]
    out['bphi'] = np.concatenate((bphi, bphi_sym),1)
    # bz
    bz = data['bz'][:, :-1, :]
    bz_sym = data['bz'][-1::-1, -1:0:-1, :]
    out['bz'] = np.concatenate((bz, bz_sym),1)
    # s
    s = data['s'][:, :-1, :]
    s_sym = data['s'][-1::-1, -1:0:-1, :]
    out['s'] = np.concatenate((s, s_sym),1)
    out['symmetrymode'] == 1
    return out

def bfield_remove_duplicate_phi(data):
    out = data
    out['br'] = data['br'][:, :-1, :]
    out['bphi'] = data['bphi'][:, :-1, :]
    out['bz'] = data['bz'][:, :-1, :]
    out['s'] = data['s'][:, :-1, :]
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

def bfield_psi_lims(data, h5fn):
    try:
        a5 = Ascotpy(h5fn=h5fn)
        a5.init(bfield=1)
    except OSError as err:
        raise err
    # Initial guess for psi0
    psi0 = np.amin(data['s'])
    for i in range(data['axis_phi'].size):
        psii = fmin(
            lambda x: a5.eval_bfield(x[0],x[1],x[2],0,evalpsi=1)['psi'],
            [data['axis_r'][i], data['axis_phi'][i], data['axis_z'][i]],
            disp=0, full_output=1)[1]
        psi0 = np.amin([psi0, psii])
    # Initial guess for psi1
    psi1 = np.amax(data['s'])
    for i in range(data['axis_phi'].size):
        psii = fmin(
            lambda x: -a5.eval_bfield(x[0],x[1],x[2],0,evalpsi=1)['psi'],
            [data['axis_r'][i], data['axis_phi'][i], data['axis_z'][i]],
            disp=0, full_output=1)[1]
        psii = -psii            # Back to a positive value
        psi1 = np.amax([psi1, psii])
    return psi0, psi1
