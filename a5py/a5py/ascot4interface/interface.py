"""
Converting ASCOT4 input files to ASCOT5 input HDF5.
"""
from ascot4_particles import *
from ascot4_magn_bkg import *
from ascot4_plasma import *
from ascot4_erad import *
from ascot4_wall_2d import *
from ascot4_wall_3d import *


def run(a4folder, h5fn, overwrite=True):
    """
    Convert ASCOT4 input files to ASCOT5 input HDF5 file.

    Parameters
    ----------
    
    a4folder : str
        Path to folder where ASCOT4 input files are located. Input
        files should be named "input.magn_bkg", "input.magn_header",  
        "input.plasma_1d",   "input.wall_2d",   "input.wall_3d", 
        "input.particles" and "input.h5".

    h5fn : str
        Full path to output HDF5 file.
    overwrite : bool, optional
        Indicates whether any existing fields within HDF5 file
        will bee overwritten. Default is true.

    Notes
    -----

    Standard ASCOT4 input does not contain electric field so a dummy
    field is written instead (if radial electric field is not given).

    Avoid having multiple files for same input, e.g. 3D wall both in
    input.h5 and input.wall_3d as this is not supported. Having both
    input.wall_2d and input.wall_3d is okay as these are different 
    inputs.
    """

    f = h5py.File(h5fn, 'r')

    if a4folder[-1] != "/":
        a4folder += "/"

    # Markers.
    if overwrite or (not "markers" in f):
        fname = a4folder + "input.particles"
        if (os.path.isfile(fname)):
            data = read_particles(fname)
            markers.write_hdf5(h5fn, data)

    # Magnetic field.
    if overwrite or (not "bfield" in f):
        fnamebkg = a4folder + "input.magn_bkg"
        fnamehdr = a4folder + "input.magn_header"
        fnameh5  = a4folder + "input.h5"
        if (os.path.isfile(fnamebkg)) and (os.path.isfile(fnamehdr)):
            data = read_magn_bkg(fnamebkg, fnamehdr)
            
            if data["n_phi"] > 1:
                B_3D.write_hdf5(h5fn, data)
            else:
                B_2D.write_hdf5(h5fn, data)

        elif os.path.isfile(fnameh5):
            data = read_magn_bkg_stellarator(fnameh5)
            B_ST.write_hdf5(h5fn, data)

    # Plasma profiles.
    if overwrite or (not "plasma" in f):
        fname = a4folder + "input.plasma_1d"
        if (os.path.isfile(fname)):
            data = read_plasma(fname)
            plasma_1D.write_hdf5(h5fn, data)

    # Electric field.
    if overwrite or (not "efield" in f):
        fname = a4folder + "input.erad"
        if (os.path.isfile(fname)):
            data = read_erad(fname)
            E_1D.write_hdf5(h5fn, data)
        else:
            E = np.array([0, 0, 0])
            E_TC.write_hdf5(h5fn, E)

    # Wall.
    if overwrite or (not "wall" in f):
        fname = a4folder + "input.wall_2d"
        if (os.path.isfile(fname)):
            data = read_wall_2d(fname)
            wall_2D.write_hdf5(h5fn, data)

        fname   = a4folder + "input.wall_3d"
        fnameh5 = a4folder + "input.h5"
        if (os.path.isfile(fname)):
            data = read_wall_3d(fname)
            wall_3D.write_hdf5(h5fn, data)
        elif (os.path.isfile(fnameh5)):
            data = read_wall_3d_hdf5(fname)
            wall_3D.write_hdf5(h5fn, data)

    f.close()

