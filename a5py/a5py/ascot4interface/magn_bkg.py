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

        data['phimap_tor'] = np.array([int(float(number)) for number in fh.readline().split()])
        data['phimap_pol'] = np.array([int(float(number)) for number in fh.readline().split()])

        h5data = np.array(fh.read().split(), dtype=float).flatten()

        sz2d = nr*nz
        sz3d = nr*nz*data['nPhi']
        data['psi']  = h5data[:sz2d].reshape(nz,nr)
        print(np.shape(h5data),nz,data['nPhi'],nr,sz2d,sz3d)
        data['br']   = h5data[sz2d+0*sz3d:sz2d+1*sz3d].reshape(nz,data['nPhi'],nr).squeeze()
        data['bphi'] = h5data[sz2d+1*sz3d:sz2d+2*sz3d].reshape(nz,data['nPhi'],nr).squeeze()
        data['bz']   = h5data[sz2d+2*sz3d:sz2d+3*sz3d].reshape(nz,data['nPhi'],nr).squeeze()

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
        bstr = 'bfield/stellarator/'

        data['r']    = f[bstr+'r'][:]
        data['phi']  = f[bstr+'phi'][:]
        data['z']    = f[bstr+'z'][:]

        data['br']   = f[bstr+'br'][:]
        data['bphi'] = f[bstr+'bphi'][:]
        data['bz']   = f[bstr+'bz'][:]
        data['s']    = f[bstr+'s'][:]

        data['axis_r']   = f[bstr+'axis_R'][:]
        data['axis_phi'] = f[bstr+'axis_phi'][:]
        data['axis_z']   = f[bstr+'axis_z'][:]

        data['n_periods'] = f[bstr+'toroidalPeriods'][:]
        try:
            data['symmetrymode'] = f[bstr+'symmetrymode'][:]
        except KeyError:
            print("Warning! No symmetry mode specified in input.h5/bfield")
            print("Defaulting to stellarator symmetry")
            data['symmetrymode'] = 0
    philim = 360/data['n_periods']/(2 if data['symmetrymode'] == 0 else 1)
    if(data['symmetrymode'] == 0):
        print("Converting stellarator symmetric input to periodic.")
        data = stellarator_bfield_sector2full(data)
    elif (data['phi'][0] == np.mod(data['phi'][-1],philim)):
        print("Warning! Removing duplicate bfield data point.")
        data = bfield_remove_duplicate_phi(data)
    if (data['axis_phi'][0] == np.mod(data['axis_phi'][-1],360)):
        print("Warning! Removing duplicated axis datapoint.")
        data['axis_r'] = data['axis_r'][0:-1]
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
    out['symmetrymode'] = 1
    return out

def bfield_remove_duplicate_phi(data):
    out = data
    out['br'] = data['br'][:, :-1, :]
    out['bphi'] = data['bphi'][:, :-1, :]
    out['bz'] = data['bz'][:, :-1, :]
    out['s'] = data['s'][:, :-1, :]
    return out
