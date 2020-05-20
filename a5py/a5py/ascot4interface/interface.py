"""
Converting ASCOT4 input files to ASCOT5 input HDF5.
"""
import os.path

import numpy as np
import h5py

import a5py.ascot4interface.markers  as a4markers
import a5py.ascot4interface.magn_bkg as a4magn_bkg
import a5py.ascot4interface.plasma   as a4plasma
import a5py.ascot4interface.erad     as a4erad
import a5py.ascot4interface.wall_2d  as a4wall_2d
import a5py.ascot4interface.wall_3d  as a4wall_3d
import a5py.ascot4interface.mhdinput as a4mhdinput

import a5py.ascot5io.B_2DS       as B_2DS
import a5py.ascot5io.B_3DS       as B_3DS
import a5py.ascot5io.B_STS       as B_STS
import a5py.ascot5io.N0_3D       as N0_3D
import a5py.ascot5io.options     as options
import a5py.ascot5io.plasma_1D   as plasma_1D
import a5py.ascot5io.mrk_prt     as mrk_prt
import a5py.ascot5io.mrk_gc      as mrk_gc
import a5py.ascot5io.E_TC        as E_TC
import a5py.ascot5io.E_1DS       as E_1DS
import a5py.ascot5io.wall_2D     as wall_2D
import a5py.ascot5io.wall_3D     as wall_3D
import a5py.ascot5io.boozer      as boozer
import a5py.ascot5io.mhd         as mhd
import a5py.ascot5io.ascot5tools as a5tools
import a5py.testascot.helpers    as helpers

from a5py.preprocessing.boozermaps import Boozermaps
from a5py.postprocessing.physicslib import guessMass


import a5py.preprocessing.psilims as psilims

def read_markers(a4folder, h5fn):
    fname = a4folder + "input.particles"
    if (os.path.isfile(fname)):
        data = a4markers.read_particles(fname)
        if 'vphi' in data['fieldNames']:
            # We have particles
            print("Warning! Forcing time to zero for all markers.")
            mrk_prt.write_hdf5(
                fn=h5fn, n=data["fields"]["id"].size, ids=data["fields"]["id"],
                mass=data["fields"]["mass"], charge=data["fields"]["charge"],
                r=data["fields"]["Rprt"], phi=data["fields"]["phiprt"],
                z=data["fields"]["zprt"],
                vr=data["fields"]["vR"], vphi=data["fields"]["vphi"],
                vz=data["fields"]["vz"],
                anum=data["fields"]['Anum'], znum=data["fields"]['Znum'],
                weight=data["fields"]["weight"],
                time=data["fields"]["weight"]*0 )
        elif 'energy' in data['fieldNames']:
            # We have guiding centers (theta is random)
            print("Warning! Forcing time to zero and "
                  "randomizing zeta for all markers.")
            zeta = 2*np.pi*np.random.rand(data["fields"]["id"].size)
            mrk_gc.write_hdf5(
                fn=h5fn, n=data["fields"]["id"].size, ids=data["fields"]["id"],
                mass=data["fields"]["mass"], charge=data["fields"]["charge"],
                r=data["fields"]["R"], phi=data["fields"]["phi"],
                z=data["fields"]["z"],
                energy=data["fields"]["energy"], pitch=data["fields"]["pitch"],
                zeta=zeta,
                anum=data["fields"]['Anum'], znum=data["fields"]['Znum'],
                weight=data["fields"]["weight"],
                time=data["fields"]["weight"]*0 )

def read_bfield(a4folder, h5fn):
    fnamebkg = a4folder + "input.magn_bkg"
    fnamehdr = a4folder + "input.magn_header"
    fnameh5  = a4folder + "input.h5"
    if (os.path.isfile(fnamebkg)) and (os.path.isfile(fnamehdr)):
        data = a4magn_bkg.read_magn_bkg(fnamebkg, fnamehdr)
        if data["nPhi"] > 1:
            B_3DS.write_hdf5(
                fn=h5fn,
                b_rmin=data['r'][0], b_rmax=data['r'][-1], b_nr=data['r'].size,
                b_zmin=data['z'][0], b_zmax=data['z'][-1], b_nz=data['z'].size,
                b_phimin=0, b_phimax=360, b_nphi=data['nPhi'],
                axisr=data['axis_r'], axisz=data['axis_z'],
                psi=data['psi']/(2*np.pi),
                psi0=data['psi0']/(2*np.pi), psi1=data['psi1']/(2*np.pi),
                br=data['br'], bphi=data['bphi'], bz=data['bz'])
        else:
            B_2DS.write_hdf5(
                fn=h5fn,
                rmin=data['r'][0], rmax=data['r'][-1], nr=data['r'].size,
                zmin=data['z'][0], zmax=data['z'][-1], nz=data['z'].size,
                axisr=data['axis_r'], axisz=data['axis_z'],
                psi=data['psi']/(2*np.pi),
                psi0=data['psi0']/(2*np.pi), psi1=data['psi1']/(2*np.pi),
                br=data['br']*0, bphi=data['bphi'], bz=data['bz']*0)

    elif os.path.isfile(fnameh5):
        with h5py.File(fnameh5, 'r') as f:
            if (not "/bfield" in f):
                return
        data = a4magn_bkg.read_magn_bkg_stellarator(fnameh5)
        psilim = [0, 1]
        temp_B_name = B_STS.write_hdf5(
            fn=h5fn,
            b_rmin=data['r'][0], b_rmax=data['r'][-1], b_nr=data['r'].size,
            b_zmin=data['z'][0], b_zmax=data['z'][-1], b_nz=data['z'].size,
            b_phimin=data['phi'][0], b_phimax=data['phi'][-1],
            b_nphi=data['phi'].size - 1,
            psi0=psilim[0], psi1=psilim[1],
            br=data['br'], bphi=data['bphi'], bz=data['bz'], psi=data['s'],
            axis_phimin=data['axis_phi'][0], axis_phimax=data['axis_phi'][-1],
            axis_nphi=data['axis_phi'].size-1,
            axisr=data['axis_r'], axisz=data['axis_z'])
        print("Searching for psiaxis and psisepx.")
        try:
            psilim = psilims.bfield_psi_lims(data, h5fn)
        except OSError:
            print("Error: Ascotpy initialization failed. "
                  "Is libascot.so is in current folder?")
            return
        # psi1 > 1 breaks plasma evaluation, so we only keep the lower limit
        psilim = [psilim[0], 1]
        print("New limits: [" + str(psilim[0]) + ", " + str(psilim[1]) + "]")
        a5tools.removegroup(h5fn, temp_B_name)
        B_STS.write_hdf5(
            fn=h5fn,
            b_rmin=data['r'][0], b_rmax=data['r'][-1], b_nr=data['r'].size,
            b_zmin=data['z'][0], b_zmax=data['z'][-1], b_nz=data['z'].size,
            b_phimin=data['phi'][0], b_phimax=data['phi'][-1],
            b_nphi=data['phi'].size - 1,
            psi0=psilim[0], psi1=psilim[1],
            br=data['br'], bphi=data['bphi'], bz=data['bz'], psi=data['s'],
            axis_phimin=data['axis_phi'][0], axis_phimax=data['axis_phi'][-1],
            axis_nphi=data['axis_phi'].size-1,
            axisr=data['axis_r'], axisz=data['axis_z'])

def read_plasma(a4folder, h5fn):
    fname1d = a4folder + "input.plasma_1d"
    fname2d = a4folder + "input.plasma_2d"
    if (os.path.isfile(fname1d)):
        data = a4plasma.read_plasma(fname1d)
        plasma_1D.write_hdf5(
            fn=h5fn, nrho=data['nrho'], nion=data['nion'],
            anum=data['anum'], znum=data['znum'],
            mass=data['anum'], charge=data['znum'],
            rho=data['rho'], edensity=data['ne'], etemperature=data['te'],
            idensity=data['ni'], itemperature=data['ti'])
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
        E_1DS.write_hdf5(
            fn=h5fn, nrho=data['n_rho'],
            rhomin=np.amin(data['rho']), rhomax=np.amax(data['rho']),
            dvdrho=data['dV_drho'], reff=1.0)
    else:
        E = np.array([0, 0, 0])
        E_TC.write_hdf5(fn=h5fn, exyz=E)

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
            fn=h5fn, nelements=data["flag"].size,
            x1x2x3=data['x1x2x3'], y1y2y3=data['y1y2y3'], z1z2z3=data['z1z2z3'],
            desc='fromASCOT4',
            flag=np.reshape(data['flag'],(data["flag"].size,1)))
    elif (os.path.isfile(fnameh5)):
        with h5py.File(fnameh5, 'r') as f:
            if (not "/wall" in f):
                return
        data = a4wall_3d.read_wall_3d_hdf5(fnameh5)
        wall_3D.write_hdf5(
            fn=h5fn, nelements=data["flag"].size,
            x1x2x3=data['x1x2x3'], y1y2y3=data['y1y2y3'], z1z2z3=data['z1z2z3'],
            desc='fromASCOT4',
            flag=np.reshape(data['flag'],(data["flag"].size,1)))


def read_boozer(a4folder, h5fn):
    fname = a4folder + "boozer_maps.out"
    if (os.path.isfile(fname)):
        b = Boozermaps(fname)
        b.write_hdf5(h5fn)
    else:
        boozer.write_hdf5_dummy(h5fn)


def read_mhd(a4folder, h5fn):
    fname = a4folder + "input.alfven"
    if (os.path.isfile(fname)):
        data = a4mhdinput.read_alfven(fname)
        print("No MHD phase data in Ascot4. Assuming phase = 0.")
        mhd.write_hdf5(
            fn=h5fn,
            nmode     = data["nmode"],
            nmodes    = data["nmodes"],
            mmodes    = data["mmodes"],
            amplitude = data["amplitude"],
            omega     = data["omega"],
            phase     = np.zeros(data["omega"].shape),
            alpha     = data["alpha"],
            phi       = data["phi"],
            nrho      = data["nrho"],
            rhomin    = data["rhomin"],
            rhomax    = data["rhomax"])
    else:
        mhd.write_hdf5_dummy(h5fn)


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
        N0_3D.write_hdf5_dummy(h5fn)

    # Boozer data
    if overwrite or (not "boozer" in groups):
        read_boozer(a4folder, h5fn)

    # MHD input
    if overwrite or (not "mhd" in groups):
        read_mhd(a4folder, h5fn)

    # Wall.
    if overwrite or (not "wall" in groups):
        read_wall(a4folder, h5fn)

    # Options 
    if overwrite or (not "options" in groups):
        odict = options.generateopt()
        helpers.clean_opt(odict)
        #GCF
        odict["SIM_MODE"]                  = 2
        odict["FIXEDSTEP_USE_USERDEFINED"] = 1
        odict["FIXEDSTEP_USERDEFINED"]     = 1e-8
        odict["ENDCOND_SIMTIMELIM"]        = 1
        odict["ENDCOND_MAX_SIMTIME"]       = 5e-6
        odict["ENABLE_ORBIT_FOLLOWING"]    = 1
        odict["ENABLE_MHD"]                = 1
        odict["ENABLE_COULOMB_COLLISIONS"] = 1

        options.write_hdf5(h5fn, odict)       

