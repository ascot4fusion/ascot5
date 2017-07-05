from pylab import *
import numpy as np

def read_erad(fn):

    f = open(fn,'r')
    str = {'comm1' : f.readline()} #Skip comment line
    str['n_rho'] = f.readline()

    data = loadtxt(f)

    str['rho'] = data[:,0]
    str['dV_drho'] = data[:,1]
    # For data in format dV/rho, we can ignore effective minor radius
    
    return str
    
def write_erad(f, e):

    f.create_group('efield/erad')
    f['efield/erad'].attrs['n_rho'] = e['n_rho']
    f['efield/erad'].attrs['r_eff'] = 1.0
    f['efield/erad'].attrs['rho_min'] = np.amin(e['rho'])
    f['efield/erad'].attrs['rho_max'] = np.amax(e['rho'])
    f.create_dataset('efield/erad/rho', data=e['rho'])
    f.create_dataset('efield/erad/dV_drho', data=e['dV_drho'])

