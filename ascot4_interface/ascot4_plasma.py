from pylab import *
import numpy as np

def read_plasma(fn):
    f = open(fn,'r')
    c1 = f.readline()
    c2 = f.readline()
    c3 = f.readline()
    l = f.readline().split()
    if(len(l) == 4 and l[3] == 'Nrad,Nion'):
        f.seek(0)
        return read_1d(f)
    elif( len(l) == 8 and l[3] == 'size'):
        f.seek(0)
        return read_2d(f)
    else:
        print 'unrecognized first line after comments'
        print l

def read_1d(fh):
    str = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nrho,nion = map(int,fh.readline().split()[:2])
    str['znum'] = map(int,fh.readline().split()[:nion])
    str['anum'] = map(int,fh.readline().split()[:nion])
    str['coll'] = map(int,fh.readline().split()[:nion+1])
    str['nrho'] = nrho
    str['nion'] = nion
    fieldnames = fh.readline().split()[0:-1:2]
    data = loadtxt(fh)
    for i,name in enumerate(fieldnames):
        str[name.lower()] = data[:,i]

    fh.close()
    return str

def read_2d(fh):
    str = {'comm1' : fh.readline(),'comm2' : fh.readline(),'comm3' : fh.readline()}
    nr,nz = map(int,fh.readline().split()[:2])
    rmin,rmax,zmin,zmax = map(float,fh.readline().split()[:4])
    nion,rho2d = map(float,fh.readline().split()[:2])
    nion = int(nion)
    str['znum'] = map(int,fh.readline().split()[:nion])
    str['anum'] = map(int,fh.readline().split()[:nion])
    str['coll'] = map(int,fh.readline().split()[:nion+1])
    fieldnames = fh.readline().split()[0:-1:2]
    data = loadtxt(fh)
    str['r'] = linspace(rmin,rmax,nr)
    str['z'] = linspace(zmin,zmax,nz)
    for i,name in enumerate(fieldnames):
        str[name.lower()] = data[:,i].reshape(nr,nz).T

    if(i+1 < shape(data)[1]):
        print 'Not all data fields assigned to struct'
        print 'Something is probably wrong with the data file'
    
    fh.close()
    return str
    
def write_plasma_1d(f, p):

    f.create_group('plasma')
    f['plasma'].attrs["type"] = np.string_("p1d")
    f.create_dataset('plasma/Z_num', data=p['znum'])
    f.create_dataset('plasma/A_mass', data=p['anum'])
    f['plasma'].attrs['n_ions'] = p['nion']
    f['plasma'].attrs['n_neutrals'] = 0 # No neutral data in ascot4 input

    # 1D plasma properties
    f.create_dataset('plasma/1D/rho', data=p['rho'])
    f.create_dataset('plasma/1D/temp_0', data=np.zeros(p['rho'].shape))
    f.create_dataset('plasma/1D/dens_0', data=np.zeros(p['rho'].shape))
    f.create_dataset('plasma/1D/temp_e', data=p['te']*8.6173e-05)
    f.create_dataset('plasma/1D/dens_e', data=p['ne'])
    f.create_dataset('plasma/1D/temp_i', data=p['ti1']*8.6173e-05)
    dens_i = np.array([p['ni'+str(i)] for i in range(1,p['nion']+1)])
    f.create_dataset('plasma/1D/dens_i', data=np.transpose(dens_i))
    f['plasma/1D'].attrs['n_rho'] = p['nrho']

