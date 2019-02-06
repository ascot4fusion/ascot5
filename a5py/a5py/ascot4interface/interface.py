"""
Converting ASCOT4 input files to ASCOT5 input HDF5.
"""
import os.path

import numpy as np
import h5py

import a5py.ascot4interface.markers as a4markers
import a5py.ascot4interface.magn_bkg as a4magn_bkg
import a5py.ascot4interface.plasma as a4plasma
import a5py.ascot4interface.erad as a4erad
import a5py.ascot4interface.wall_2d as a4wall_2d
import a5py.ascot4interface.wall_3d as a4wall_3d

import a5py.ascot5io.B_2DS as B_2DS
import a5py.ascot5io.B_3DS as B_3DS
import a5py.ascot5io.B_STS as B_STS
import a5py.ascot5io.N0_3D as N0_3D
import a5py.ascot5io.plasma_1D as plasma_1D
import a5py.ascot5io.mrk_prt as mrk_prt
import a5py.ascot5io.mrk_gc as mrk_gc
import a5py.ascot5io.E_TC as E_TC
import a5py.ascot5io.E_1DS as E_1DS
import a5py.ascot5io.wall_2D as wall_2D
import a5py.ascot5io.wall_3D as wall_3D

from a5py.postprocessing.physicslib import guessMass

def read_markers(a4folder, h5fn):
    fname = a4folder + "input.particles"
    if (os.path.isfile(fname)):
        data = a4markers.read_particles(fname)
        if 'charge' not in data['fieldNames']:
            print("Converting Znum to charge.")
            data["fields"]['charge'] = data["fields"]['Znum'].astype('float')
        if 'mass' not in data['fieldNames']:
            print("Converting Anum to mass.")
            data["fields"]['mass'] = np.array(
                list(map(guessMass,
                         data["fields"]['Anum'],
                         data["fields"]['Znum'],
                         data["fields"]['charge'])))
        if 'id' not in data['fieldNames']:
            print("Generating unique ids.")
            data["fields"]['id'] = np.array(
                range(1,data["fields"]['charge'].size + 1))
        if (min(data["fields"]["id"]) <= 0):
            zero_ind = np.where(data["fields"]["id"] == 0)[0];
            data["fields"]["id"][zero_ind] = max(data["fields"]["id"] ) + 1
            print("Converting id 0 to new unique id: "
                  + str(int(max(data["fields"]["id"]))))
        if 'vphi' in data['fieldNames']:
            # We have particles
            print("Warning! Forcing time to zero for all markers.")
            mrk_prt.write_hdf5(
                h5fn, data["fields"]["id"].size, data["fields"]["id"],
                data["fields"]["mass"], data["fields"]["charge"],
                data["fields"]["Rprt"], data["fields"]["phiprt"],
                data["fields"]["zprt"],
                data["fields"]["vR"], data["fields"]["vphi"],
                data["fields"]["vz"],
                data["fields"]["weight"], data["fields"]["weight"]*0)
        elif 'energy' in data['fieldNames']:
            # We have guiding centers (theta is random)
            print("Warning! Forcing time to zero and "
                  "randomizing theta for all markers.")
            theta = 2*np.pi*np.random.rand(data["fields"]["id"].size)
            mrk_gc.write_hdf5(
                h5fn, data["fields"]["id"].size, data["fields"]["id"],
                data["fields"]["mass"], data["fields"]["charge"],
                data["fields"]["R"], data["fields"]["phi"],
                data["fields"]["z"],
                data["fields"]["energy"], data["fields"]["pitch"], theta,
                data["fields"]["weight"], data["fields"]["weight"]*0 )

def read_bfield(a4folder, h5fn):
    fnamebkg = a4folder + "input.magn_bkg"
    fnamehdr = a4folder + "input.magn_header"
    fnameh5  = a4folder + "input.h5"
    if (os.path.isfile(fnamebkg)) and (os.path.isfile(fnamehdr)):
        data = a4magn_bkg.read_magn_bkg(fnamebkg, fnamehdr)
        if data["nPhi"] > 1:
            B_3DS.write_hdf5(
                h5fn,
                data['r'][0], data['r'][-1], data['r'].size,
                data['z'][0], data['z'][-1], data['z'].size,
                0, 360, data['nPhi'],
                data['axis_r'], data['axis_z'], data['psi']/(2*np.pi),
                data['psi0']/(2*np.pi), data['psi1']/(2*np.pi),
                data['br'], data['bphi'], data['bz'])
        else:
            B_2DS.write_hdf5(
                h5fn,
                data['r'][0], data['r'][-1], data['r'].size,
                data['z'][0], data['z'][-1], data['z'].size,
                data['axis_r'], data['axis_z'], data['psi']/(2*np.pi),
                data['psi0']/(2*np.pi), data['psi1']/(2*np.pi),
                data['br']*0, data['bphi'], data['bz']*0)

    elif os.path.isfile(fnameh5):
        with h5py.File(fnameh5, 'r') as f:
            if (not "/bfield" in f):
                return
        data = a4magn_bkg.read_magn_bkg_stellarator(fnameh5)
        if (data['axis_phi'][0] == np.mod(data['axis_phi'][-1],360)):
            # Remove duplicated datapoint
            print("Warning! Removing duplicated axis datapoint.")
            data['axis_r'] = data['axis_r'][0:-1]
            data['axis_phi'] = data['axis_phi'][0:-1]
            data['axis_z'] = data['axis_z'][0:-1]
        B_STS.write_hdf5(
            h5fn,
            data['r'][0], data['r'][-1], data['r'].size,
            data['z'][0], data['z'][-1], data['z'].size,
            data['phi'][0], data['phi'][-1], data['phi'].size,
            data['br'], data['bphi'], data['bz'], data['s'],
            data['n_periods'],
            data['axis_phi'][0], data['axis_phi'][-1], data['axis_phi'].size,
            data['axis_r'], data['axis_z'],
            sym_mode=data['symmetrymode'],
            psiaxis=0, psisepx=1)

def read_plasma(a4folder, h5fn):
    fname1d = a4folder + "input.plasma_1d"
    fname2d = a4folder + "input.plasma_2d"
    if (os.path.isfile(fname1d)):
        data = a4plasma.read_plasma(fname1d)
        # Make sure the input is linearly spaced. If not, interpolate
        tol = 1.0001
        diff = np.diff(data['rho'])
        if ( max(diff)/min(diff) > tol):
            print("Warning! Interpolating plasma data to uniform grid")
            new_rho = np.linspace(np.amin(data['rho']),
                                  np.amax(data['rho']),
                                  data['nrho'])
            data['ne'] = np.interp(new_rho, data['rho'], data['ne'])
            data['te'] = np.interp(new_rho, data['rho'], data['te'])
            for i in range(1, data['nion']+1):
                data['ni'+str(i)] = np.interp(new_rho, data['rho'],
                                              data['ni'+str(i)])
            data['ti1'] = np.interp(new_rho, data['rho'], data['ti1'])
            data['rho'] = new_rho
        dens_i = np.array([data['ni'+str(i)] for i in range(1,data['nion']+1)])
        plasma_1D.write_hdf5(
            h5fn, data['nrho'], data['nion'], data['znum'], data['anum'],
            data['rho'], data['ne'], data['te'], dens_i, data['ti1'])
    if (os.path.isfile(fname2d)):
        data = a4plasma.read_plasma(fname2d)
        dens_i = np.array([data['ni'+str(i)] for i in range(1,data['nion']+1)])
        print("2D plasma data not yet implemented for ASCOT4. "
              "Skipping 2D plasma input.")
        # plasma_2D.write_hdf5(h5fn, data['nrho'], data['nion'], data['znum'],
        #                      data['anum'],
        #                      data['rho'], np.zeros(data['rho'].shape),
        #                      np.zeros(data['rho'].shape),
        #                      data['ne'], data['te'], dens_i, data['ti1'])

def read_erad(a4folder, h5fn):
    fname = a4folder + "input.erad"
    if (os.path.isfile(fname)):
        data = a4erad.read_erad(fname)
        # Make sure the input is linearly spaced. If not, interpolate
        tol = 1.0001
        diff = np.diff(data['rho'])
        if ( max(diff)/min(diff) > tol):
            print("Warning! Interpolating dV_drho to uniform grid")
            new_rho = np.linspace(np.amin(data['rho']),
                                  np.amax(data['rho']),
                                  data['n_rho'])
            data['dV_drho'] = np.interp(new_rho, data['rho'],
                                        data['dV_drho'])
            data['rho'] = new_rho
        E_1DS.write_hdf5(
            h5fn, int(data['n_rho']), np.amin(data['rho']),
            np.amax(data['rho']), data['dV_drho'], 1.0)
    else:
        E = np.array([0, 0, 0])
        E_TC.write_hdf5(h5fn, E)

def read_wall(a4folder, h5fn):
    fname = a4folder + "input.wall_2d"
    if (os.path.isfile(fname)):
        data = a4wall_2d.read_wall_2d(fname)
        wall_2D.write_hdf5(h5fn, data['r'].size, data['r'], data['z'])
    fname   = a4folder + "input.wall_3d"
    fnameh5 = a4folder + "input.h5"
    if (os.path.isfile(fname)):
        data = a4wall_3d.read_wall_3d(fname)
        wall_3D.write_hdf5(
            h5fn, data['id'].size, data['x1x2x3'],
            data['y1y2y3'], data['z1z2z3'], data['id'])
    elif (os.path.isfile(fnameh5)):
        with h5py.File(fnameh5, 'r') as f:
            if (not "/wall" in f):
                return
        data = a4wall_3d.read_wall_3d_hdf5(fnameh5)
        wall_3D.write_hdf5(
            h5fn, data['id'].size, data['x1x2x3'],
            data['y1y2y3'], data['z1z2z3'], data['id'])

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

    if not os.path.isfile(h5fn):
        with h5py.File(h5fn, 'a') as f:
            pass

    if a4folder[-1] != "/":
        a4folder += "/"

    with h5py.File(h5fn, 'r') as f:
        groups = f.keys()

    if overwrite or (not "markers" in groups):
        read_markers(a4folder, h5fn)

    # Magnetic field.
    if overwrite or (not "bfield" in groups):
        read_bfield(a4folder, h5fn)

    # Plasma profiles.
    if overwrite or (not "plasma" in groups):
        read_plasma(a4folder, h5fn)

    # Electric field.
    if overwrite or (not "efield" in groups):
        read_erad(a4folder, h5fn)

    # Neutral density
    if overwrite or (not "neutral" in groups):
        # No ASCOT4 neutral density
        N0 = np.array([ [ [0,0] , [0,0] ], [ [0,0] , [0,0] ] ])
        N0_3D.write_hdf5(h5fn, -1, 1, 2, -1, 1, 2, 0, 2*np.pi, 2, N0)

    # Wall.
    if overwrite or (not "wall" in groups):
        read_wall(a4folder, h5fn)
