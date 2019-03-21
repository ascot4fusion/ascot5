import numpy as np

def read_erad(fn):

    with open(fn,'r') as f:
        data = {'comm1' : f.readline()} #Skip comment line
        data['n_rho'] = f.readline().split()[0]

        h5data = np.loadtxt(f)

        data['rho'] = h5data[:,0]
        data['dV_drho'] = h5data[:,1]
        # For data in format dV/rho, we can ignore effective minor radius

    return data
