"""
Module for generating Poincare run inputs.

File: poincare.py
"""
import sys
import subprocess
import numpy as np

import a5py.ascot5io.mrk_fl      as marker
import a5py.ascot5io.mrk_gc      as markergc
import a5py.ascot5io.options     as options
import a5py.ascot5io.ascot5tools as tools

from a5py.physlib.species import species as speciesdata

from a5py.ascotpy import Ascotpy
from a5py.testascot.helpers import clean_opt

from scipy import interpolate

def gen_poincaremarkers(fn, tor, pol, rhomin, rhomax, ntrace, time, pitch,
                        desc=None, active=False, species=None, energy=None):
    """
    Generate Poincare tracers and return the marker struct.
    """
    if tor < 0:
        tor = 360 * np.random.rand(ntrace,)
    else:
        tor = tor * np.ones((ntrace,))

    if pol < 0:
        pol = 360 * np.random.rand(ntrace,)
    else:
        pol = pol * np.ones((ntrace,))

    # Get (R,z) coordinates of marker initial positions
    a5 = Ascotpy(fn)
    a5.init(bfield=True)
    R = np.zeros((ntrace,1))
    z = np.zeros((ntrace,1))
    for i in range(ntrace):
        rz = a5.get_rhotheta_rz( np.linspace(rhomin, rhomax, ntrace),
                                 pol[i], tor[i], time )
        R[i] = (rz[0])[i]
        z[i] = (rz[1])[i]

    a5.free(bfield=True)

    # Set field line marker input
    mrk = {}
    mrk["n"]      = R.size
    mrk["ids"]    = np.linspace(1, R.size,  R.size)
    mrk["r"]      = R
    mrk["z"]      = z
    mrk["phi"]    = tor
    mrk["pitch"]  = pitch * np.ones(R.shape)
    mrk["time"]   = time  * np.ones(R.shape)
    mrk["weight"] = 1     * np.ones(R.shape)

    if desc is not None:
        if species is not None:
            species       = speciesdata(species)
            mrk["energy"] = energy  * np.ones(R.shape)
            mrk["zeta"]   = 2*np.pi * np.random.rand(R.size,)
            mrk["anum"]   = species["anum"]   * np.ones(R.shape)
            mrk["znum"]   = species["znum"]   * np.ones(R.shape)
            mrk["mass"]   = species["mass"].to("amu").d   * np.ones(R.shape)
            mrk["charge"] = species["charge"] * np.ones(R.shape)
            mqid = markergc.write_hdf5(fn, **mrk, desc=desc)
        else:
            mqid = marker.write_hdf5(fn, **mrk, desc=desc)
        if active:
            tools.call_ascot5file(fn, "set_active", mqid)

    return mrk


def gen_poincareoptions(fn, tor, pol, porb, torb, rhoend, mhd, cputmax,
                        desc=None, active=False, simmode=4):
    """
    Generate Poincare options and return the options struct.
    """
    odict = options.generateopt()
    clean_opt(odict)

    odict["SIM_MODE"]        = simmode
    odict["ENABLE_ADAPTIVE"] = 1

    odict["ENABLE_ORBIT_FOLLOWING"] = 1
    odict["ENABLE_ORBITWRITE"]      = 1
    odict["ORBITWRITE_MODE"]        = 0
    odict["ORBITWRITE_TOROIDALANGLES"] = tor
    odict["ORBITWRITE_POLOIDALANGLES"] = pol
    odict["ORBITWRITE_NPOINT"] = porb*pol.size + torb*tor.size

    odict["ENDCOND_MAX_POLOIDALORBS"] = porb
    odict["ENDCOND_MAX_TOROIDALORBS"] = torb

    odict["ENDCOND_MAXORBS"] = 1
    if rhoend < 0:
        odict["ENDCOND_WALLHIT"] = 1
    else:
        odict["ENDCOND_RHOLIM"]  = 1
        odict["ENDCOND_MAX_RHO"] = rhoend

    if cputmax:
        odict["ENDCOND_CPUTIMELIM"]  = 1
        odict["ENDCOND_MAX_CPUTIME"] = cputmax

    odict["ENABLE_MHD"] = mhd

    if desc is not None:
        oqid = options.write_hdf5(fn, odict, desc=desc)
        if active:
            tools.call_ascot5file(fn, "set_active", oqid)

    return odict



if __name__ == "__main__":

    ## Interpret parameters and flags ##

    info = "\n\
    Generates marker and options input for a Poincare plot run.\n\
    \n\
    Parameters (default values in parenthesis):\n\
    -rhoend  : maximum rho at which markers are terminated. Negative if they \n\
               are stopped at the wall (1)\n\
    -torb    : number of toroidal plane crossings before termination (500)\n\
    -porb    : number of poloidal plane crossings before termination (500)\n\
    -cputmax : enable CPU time limit and use this value as the limit (0) \n\
    -ntrace  : number of field line tracers (100)\n\
    -rhomin  : minimum field line tracer initial rho (0)\n\
    -rhomax  : maximum field line tracer initial rho (1)\n\
    -tor0    : field line tracer initial toroidal angle. If negative, markers\n\
               are placed at random angles (-1)\n\
    -pol0    : field line tracer initial poloidal angle. If negative, markers\n\
               are placed at random angles (0)\n\
    -tor     : toroidal angles [deg] of poloidal Poincare planes. If negative,\n\
               poloidal Poincare planes are disabled (0)\n\
    -pol     : poloidal angles [deg] of toroidal Poincare planes. If negative\n\
               toroidal Poincare planes are disabled (0)\n\
    -pitch   : direction, co [1] or counter [-1], field lines are traced (1) \n\
    -mhd     : enable MHD (0) \n\
    -time    : time-slice at which field lines are traced (0) \n\
    -fn      : input HDF5 file (ascot.h5)\n\
    \n\
    Example:\n\
    python poincare.py -torb 2000 -pol \"10 20\" -fn input.h5\n\
    \n\
    To simulate particles instead of field lines, provide two additional \n\
    parameters:\n\
    -energy  : particle energy in eV.\n\
    -species : particle species (see a5py/physlib/species.py). \n\
    \n\
    "

    # Default values
    rhoend  = 1
    torb    = 500
    porb    = 500
    ntrace  = 100
    rhomin  = 0
    rhomax  = 1
    tor0    = -1
    pol0    = 0
    tor     = np.array([0],dtype="f8")
    pol     = np.array([0],dtype="f8")
    pitch   = 1
    mhd     = 0
    time    = 0
    cputmax = 0
    fn      = "ascot"

    energy  = None
    species = None

    if len(sys.argv) == 1:
        print(info)
        print("\n Continue with default parameters? (y/n)")
        while True:
            yn = input("")
            if yn == "y":
                break
            if yn == "n":
                exit(0)

    # Read user input
    unknownoption = False
    for j in range(1,len(sys.argv), 2):
        i = sys.argv[j]
        if i == "-rhoend":
            rhoend = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-torb":
            torb = np.fromstring(sys.argv[j+1], dtype=np.uint,sep=' ')
        elif i == "-porb":
            porb = np.fromstring(sys.argv[j+1], dtype=np.uint,sep=' ')
        elif i == "-ntrace":
            ntrace = int(np.fromstring(sys.argv[j+1], dtype=np.uint,sep=' '))
        elif i == "-rhomin":
            rhomin = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-rhomax":
            rhomax = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-tor0":
            tor0 = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-pol0":
            pol0 = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-tor":
            tor = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-pol":
            pol = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-pitch":
            pitch = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-mhd":
            mhd = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-cputmax":
            cputmax = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-time":
            time = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-fn":
            fn = sys.argv[j+1]
        elif i == "-energy":
            energy = np.fromstring(sys.argv[j+1], dtype=np.float,sep=' ')
        elif i == "-species":
            species = sys.argv[j+1]
        else:
            unknownoption = True

        if unknownoption:
            print(info)
            print("Error: Unknown option.")
            exit(0)

    if (species is not None and energy is None) or \
       (species is None and energy is not None):
        print("Error: Please provide both energy and species.")
        exit(0)

    simmode = 4
    if species is not None:
        simmode=2

    gen_poincaremarkers(fn, tor0, pol0, rhomin, rhomax, ntrace, time, pitch,
                        desc="Poincare", active=True,
                        species=species, energy=energy)
    gen_poincareoptions(fn, tor, pol, porb, torb, rhoend, mhd, cputmax,
                        desc="Poincare", active=True, simmode=simmode)

    print("Poincare run set.")
