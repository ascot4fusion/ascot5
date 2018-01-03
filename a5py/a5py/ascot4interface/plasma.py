from scipy import *
import numpy as np


def read_plasma(fn):
    f = open(fn,'r')
    c1 = f.readline()
    c2 = f.readline()
    c3 = f.readline()
    l = f.readline().split()
    if(fn[-2:] == '1d'):
        f.seek(0)
        return read_1d(f)
    elif(fn[-2:] == '2d'):
        f.seek(0)
        return read_2d(f)
    else:
        print('unrecognized first line after comments')
        print(l)

def read_1d(fh):
    pls = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nrho,nion = list(map(int,fh.readline().split()[:2]))
    pls['znum'] = np.array(list(map(int,fh.readline().split()[:nion])))
    pls['anum'] = np.array(list(map(int,fh.readline().split()[:nion])))
    pls['coll'] = np.array(list(map(int,fh.readline().split()[:nion+1])))
    pls['nrho'] = nrho
    pls['nion'] = nion
    fh.readline() # ignore headers
    data = loadtxt(fh)
    pls['rho'] = np.array(data[:,0])
    pls['te'] = np.array(data[:,1])
    pls['ne'] = np.array(data[:,2])
    pls['vtor'] = np.array(data[:,3])
    pls['ti1'] = np.array(data[:,4])
    for i in range(0,nion):
        pls['ni'+str(i+1)] = np.array(data[:,5+i])
    fh.close()
    return pls

def read_2d(fh):
    str = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nr,nz = map(int,fh.readline().split()[:2])
    rmin,rmax,zmin,zmax = list(map(float,fh.readline().split()[:4]))
    nion,rho2d = list(map(float,fh.readline().split()[:2]))
    nion = int(nion)
    str['znum'] = list(map(int,fh.readline().split()[:nion]))
    str['anum'] = list(map(int,fh.readline().split()[:nion]))
    str['coll'] = list(map(int,fh.readline().split()[:nion+1]))
    fieldnames = fh.readline().split()[0:-1:2]
    data = loadtxt(fh)
    str['r'] = linspace(rmin,rmax,nr)
    str['z'] = linspace(zmin,zmax,nz)
    for i,name in enumerate(fieldnames):
        str[name.lower()] = data[:,i].reshape(nr,nz).T

    if(i+1 < shape(data)[1]):
        print('Not all data fields assigned to struct')
        print('Something is probably wrong with the data file')
    
    fh.close()
    return str
    
def write_plasma_1d(fn, p):
    dens_i = np.array([p['ni'+str(i)] for i in range(1,p['nion']+1)])
    p1D.write_hdf5(fn, p['nrho'], p['nion'], p['znum'], p['znum'], p['rho'], np.zeros(p['rho'].shape), np.zeros(p['rho'].shape), p['ne'], p['te'], dens_i, p['ti1'])


