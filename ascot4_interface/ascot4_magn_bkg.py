import sys
sys.path.append('ui')

from pylab import *
import ui_B_2D

def read_magn_bkg(fn):
    str = dict()
    fh = open(fn)
    tmp = map(float,fh.readline().split())
    str['phi0']       = tmp[0]
    str['nSector']    = int(tmp[1])
    str['nPhi']       = int(tmp[2])
    str['nCoil']      = int(tmp[3])
    str['zeroAtCoil'] = int(tmp[4])    

    r1,r2,nr = map(float,fh.readline().split())
    nr = int(nr)
    z1,z2,nz = map(float,fh.readline().split())
    nz = int(nz)

    str['r'] = linspace(r1,r2,nr)
    str['z'] = linspace(z1,z2,nz)
    if(str['nPhi'] > 1):
        dphi = 360.0 / ( str['nSector'] * str['nPhi'] )
        str['phi'] = ( linspace( str['phi0'] + dphi * 0.5,
                                 str['phi0'] + 360.0 / str['nSector'] - dphi * 0.5,
                                 str['nPhi'] ) )

    str['phimap_tor'] = array(map(int,fh.readline().split()))
    str['phimap_pol'] = array(map(int,fh.readline().split()))

    #data = loadtxt(fh).flatten()
    data = array(fh.read().split(), dtype=float).flatten()

    sz2d = nr*nz
    sz3d = nr*nz*str['nPhi']
    str['psi']  = data[:sz2d].reshape(nz,nr)
    print shape(data),nz,str['nPhi'],nr,sz2d,sz3d
    str['br']   = data[sz2d+0*sz3d:sz2d+1*sz3d].reshape(nz,str['nPhi'],nr).squeeze()
    str['bphi'] = data[sz2d+1*sz3d:sz2d+2*sz3d].reshape(nz,str['nPhi'],nr).squeeze()
    str['bz']   = data[sz2d+2*sz3d:sz2d+3*sz3d].reshape(nz,str['nPhi'],nr).squeeze()

    fh.close()
    
    return str

def read_magn_header(fn):
    str = dict()
    return str

def write_magn_bkg(f, m):
    ui_B_2D.write_hdf5('ascot.h5',
                       np.array([m['r'][0], m['r'][-1]]),
                       np.array([m['z'][0], m['z'][-1]]),
                       m['psi']/(2*np.pi), m['br']*0, m['bphi'], m['bz']*0,
                       np.array([6.2, 0.6]), np.array([-12, 0]))
    return
