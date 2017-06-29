import h5py
from ascot4_particles import *
from ascot4_magn_bkg import *
from ascot4_plasma import *
from ascot4_wall_2d import *
from ascot4_wall_3d import *

def main():
    overwrite_fields = False
    h5file = 'ascot.h5'
    f = h5py.File(h5file, 'a') # Open for reading or writing    

    # Particle input to h5
    if overwrite_fields or (not 'inistate' in f):
        if 'inistate' in f:
            del f['inistate']
        fname = 'input.particles'
        data = read_particles(fname)
        write_particles(f, data)
    # Magnetic field to h5

    #TODO inputs for all kinds of different ASCOT4 magn_bkg
    # at least 2D and 3D
    fname = 'input.magn_bkg'
    data = read_magn_bkg(fname)
    write_magn_bkg(f, data)

    # plasma profile
    if overwrite_fields or (not 'plasma' in f):
        if 'plasma' in f:
            del f['plasma']
        fname = 'input.plasma_1d'
        data = read_plasma(fname)
        write_plasma_1d(f, data)
        f['plasma'].attrs['type'] = '1D'

    # 2D wall
    if overwrite_fields or (not 'wall/2D' in f):
        if 'wall/2D' in f:
            del f['wall/2D']
        fname = 'input.wall_2d'
        data = read_wall_2d(fname)
        write_wall_2d(f, data)
        f['wall'].attrs['type'] = '2D'

    # 3D wall
    if overwrite_fields or (not 'wall/3D' in f):
        if 'wall/3D' in f:
            del f['wall/3D']
        fname = 'input.wall_3d'
        data = read_wall_3d(fname)
        write_wall_3d(f, data)
        f['wall'].attrs['type'] = '3D'
        
if __name__ == '__main__':
    main()
