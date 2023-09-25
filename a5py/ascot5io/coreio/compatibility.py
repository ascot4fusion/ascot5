"""Make older ASCOT5 HDF5 files compatible with newer versions.

Whenever a new version is released, add a new function here that converts
previous version to current version. Note that you might have to modify the HDF5
directly as this tool should be independent of ascot5io (which is subject to
changes).

Also update the CURRENT_VERSION number.
"""
import h5py
import unyt
import numpy as np
import subprocess

from .fileapi             import get_qid
from a5py.ascot5io.boozer import Boozer
from a5py.ascot5io.mhd    import MHD_STAT
from a5py.ascot5io.asigma import Asigma_loc

def convert(fn):
    """Convert an old file to the current version.

    The old file remains as it is, and the converted version will be named
    fnin_{CURRENT_VERSION}.h5.

    Parameters
    ----------
    fn : str
        Name of the input file.
    """
    fnout = fn[:-3] + "_" + "5" + ".h5"
    subprocess.call(["cp", fn, fnout])

    print("\nUpdating file to version 5.2 (if necessary)\n")
    _convert1to2(fnout)
    print("\nUpdating file to version 5.3 (if necessary)\n")
    _convert2to3(fnout)
    print("\nUpdating file to version 5.4 (if necessary)\n")
    _convert3to4(fnout)
    print("\nUpdating file to version 5.5 (if necessary)\n")
    _convert4to5(fnout)

def _convert4to5(fn):
    """Update version 4 HDF5 to version 5.

    - Adds constant of motion distribution settings.
    - Renames id -> ids and removes underscores from field names in NBI
      inputs.
    """
    with h5py.File(fn, "a") as h5:
        for opt in _loopchild(h5, "options"):
            grp = h5["options"][opt]
            if not "ENABLE_DIST_COM" in grp:
                print("Adding ENABLE_DIST_COM to %s" % opt)
                grp.create_dataset("ENABLE_DIST_COM", (1,), data=0, dtype='i8')
            if not "DIST_NBIN_MU" in grp:
                print("Adding DIST_NBIN_MU to %s" % opt)
                grp.create_dataset("DIST_NBIN_MU", (1,), data=100, dtype='i8')
            if not "DIST_MIN_MU" in grp:
                print("Adding DIST_MIN_MU to %s" % opt)
                grp.create_dataset("DIST_MIN_MU", (1,), data=0, dtype='f8')
            if not "DIST_MAX_MU" in grp:
                print("Adding DIST_MAX_MU to %s" % opt)
                grp.create_dataset("DIST_MAX_MU", (1,), data=2.5e-13,
                                   dtype='f8')
            if not "DIST_NBIN_EKIN" in grp:
                print("Adding DIST_NBIN_EKIN to %s" % opt)
                grp.create_dataset("DIST_NBIN_EKIN", (1,), data=100, dtype='i8')
            if not "DIST_MIN_EKIN" in grp:
                print("Adding DIST_MIN_EKIN to %s" % opt)
                grp.create_dataset("DIST_MIN_EKIN", (1,), data=0, dtype='f8')
            if not "DIST_MAX_EKIN" in grp:
                print("Adding DIST_MAX_EKIN to %s" % opt)
                grp.create_dataset("DIST_MAX_EKIN", (1,), data=1.0e-12,
                                   dtype='f8')
            if not "DIST_NBIN_PTOR" in grp:
                print("Adding DIST_NBIN_PTOR to %s" % opt)
                grp.create_dataset("DIST_NBIN_PTOR", (1,), data=200, dtype='i8')
            if not "DIST_MIN_PTOR" in grp:
                print("Adding DIST_MIN_PTOR to %s" % opt)
                grp.create_dataset("DIST_MIN_PTOR", (1,), data=-1e-18,
                                   dtype='f8')
            if not "DIST_MAX_PTOR" in grp:
                print("Adding DIST_MAX_PTOR to %s" % opt)
                grp.create_dataset("DIST_MAX_PTOR", (1,), data=1e-18,
                                   dtype='f8')

        for nbi in _loopchild(h5, "nbi"):
            for inj in h5["nbi"][nbi]:
                if inj == "ninj": continue
                grp = h5["nbi"][nbi][inj]
                if "id" in grp:
                    print("Removing id from %s %s" % (nbi, inj))
                    ids = grp["id"][:]
                    del grp["id"]
                    if "ids" not in grp:
                        print("Adding ids to %s %s" % (nbi, inj))
                        grp.create_dataset("ids", (1,), data=ids, dtype='i8')
                for q in ["div_h", "div_v", "div_halo_frac", "div_halo_h",
                          "div_halo_v"]:
                    if q in grp:
                        print("Removing %s from %s %s" % (q, nbi, inj))
                        val = grp[q][:]
                        del grp[q]
                        q = q.replace("_", "")
                        if q not in grp:
                            print("Adding %s to %s %s" % (q, nbi, inj))
                            grp.create_dataset(q, (1,), data=val, dtype='f8')
                if grp["mass"][()] < 1:
                    print("Converting mass to amu")
                    mass = (grp["mass"][()] / unyt.amu).v
                    grp["mass"][()] = mass
                if grp["energy"][()] < 1:
                    print("Converting energy to eV")
                    energy = (grp["energy"][()] / unyt.elementary_charge).v
                    grp["energy"][()] = energy

def _convert3to4(fn):
    """Update version 3 HDF5 to version 4.

    - Adds REVERSE_TIME, ENDCOND_IONIZED, ENDCOND_NEUTRALIZED, and ENABLE_ATOMIC
      options.
    - Renames ENDCOND_MAX_SIMTIME to ENDCOND_LIM_SIMTIME.
    - Adds dummy atomic input to existing runs.
    """
    with h5py.File(fn, "a") as h5:
        for opt in _loopchild(h5, "options"):
            grp = h5["options"][opt]
            if "ENDCOND_MAX_SIMTIME" in grp:
                print("Removing ENDCOND_MAX_SIMTIME from %s" % opt)
                simtime = grp["ENDCOND_MAX_SIMTIME"][:]
                del grp["ENDCOND_MAX_SIMTIME"]
                if not "ENDCOND_LIM_SIMTIME" in opt:
                    print("Adding ENDCOND_LIM_SIMTIME to %s" % opt)
                    grp.create_dataset("ENDCOND_LIM_SIMTIME", (1,),
                                       data=simtime, dtype='f8')
            if not "REVERSE_TIME" in grp:
                print("Adding REVERSE_TIME to %s" % opt)
                grp.create_dataset("REVERSE_TIME", (1,), data=0, dtype='f8')
            if not "ENABLE_ATOMIC" in grp:
                print("Adding ENABLE_ATOMIC to %s" % opt)
                grp.create_dataset("ENABLE_ATOMIC", (1,), data=0, dtype='i8')
            if not "ENDCOND_IONIZED" in grp:
                print("Adding ENDCOND_IONIZED to %s" % opt)
                grp.create_dataset("ENDCOND_IONIZED", (1,), data=0, dtype='i8')
            if not "ENDCOND_NEUTRALIZED" in grp:
                print("Adding ENDCOND_NEUTRALIZED to %s" % opt)
                grp.create_dataset("ENDCOND_NEUTRALIZED", (1,), data=0,
                                   dtype='i8')

        dummyasigma = False
        for run in _loopchild(h5, "results"):
            if run[:4] != "run": continue
            if "qid_asigma" not in h5[run].attrs:
                dummyasigma = True

    if dummyasigma:
        print("Adding dummy atomic input")
        dummyasigma = get_qid(Asigma_loc.write_hdf5_dummy(fn))
        with h5py.File(fn, "a") as h5:
            for run in _loopchild(h5, "results"):
                if run[:4] != "run": continue
                if "qid_asigma" not in h5[run].attrs:
                    print("Setting dummy atomic data as an input for %s" % run)
                    h5[run].attrs["qid_asigma"] = np.string_(dummyasigma)

def _convert2to3(fn):
    """Update version 2 HDF5 to version 3.

    - Adds dummy boozer and mhd groups for existing runs.
    - Adds flags to wall input.
    - Adds ENDCOND_MAX_MILEAGE from options.
    - Adds ENABLE_MHD from options.
    - Replace velocity with momentum in distribution options.
    """
    dummymhd    = False
    dummyboozer = False
    with h5py.File(fn, "a") as h5:
        for run in _loopchild(h5, "results"):
            if run[:4] != "run": continue
            if "qid_boozer" not in h5["results"][run].attrs:
                dummyboozer = True
            if "qid_mhd" not in h5["results"][run].attrs:
                dummymhd = True

    if dummyboozer:
        print("Adding dummy Boozer input")
        dummyboozer = get_qid(Boozer.write_hdf5_dummy())
        with h5py.File(fn, "a") as h5:
            for run in h5["results"]:
                if run[:4] != "run": continue
                print("Setting dummy boozer data as an input for %s" % run)
                if "qid_boozer" not in h5["results"][run].attrs:
                    h5["results"][run].attrs["qid_boozer"] = \
                        np.string_(dummyboozer)
    if dummymhd:
        print("Adding dummy MHD input")
        dummymhd = get_qid(MHD_STAT.write_hdf5_dummy())
        with h5py.File(fn, "a") as h5:
            for run in h5["results"]:
                if run[:4] != "run": continue
                if "qid_mhd" not in h5["results"][run].attrs:
                    print("Setting dummy MHD data as an input for %s" % run)
                    h5["results"][run].attrs["qid_mhd"] = np.string_(dummymhd)

    with h5py.File(fn, "a") as h5:
        for wall in _loopchild(h5, "wall"):
            if "wall_3D" in wall:
                g = h5["wall"][wall]
                if not "flag" in g:
                    print("Adding flag to %s" % wall)
                    nelements = int(g["nelements"][:])
                    flag = np.zeros(shape=(nelements,1), dtype=int)
                    g.create_dataset("flag", (nelements,1), data=flag,
                                     dtype="i4")

        for opt in _loopchild(h5, "options"):
            grp = h5["options"][opt]
            if "ENABLE_MHD" not in grp:
                print("Adding ENABLE_MHD to %s" % opt)
                grp.create_dataset("ENABLE_MHD", (1,), data=0, dtype='f8')
            if "ENDCOND_MAX_MILEAGE" not in grp:
                print("Adding ENDCOND_MAX_MILEAGE to %s" % opt)
                grp.create_dataset("ENDCOND_MAX_MILEAGE", (1,), data=100,
                                   dtype='f8')

            for coord in ["PA", "PE", "R", "PHI", "Z"]:
                helium_mass = 6.65e-27
                if "DIST_NBIN_V" + coord in grp:
                    print("Removing DIST_NBIN_V%s from %s" % (coord, opt))
                    nbin = grp["DIST_NBIN_V"+coord][:]
                    del grp["DIST_NBIN_V"+coord]
                    if not "DIST_NBIN_P"+coord:
                        print("Adding DIST_NBIN_P%s to %s" % (coord, opt))
                        grp.create_dataset("DIST_NBIN_P"+coord, (1,), data=nbin,
                                           dtype='f8')
                if "DIST_MIN_V"+coord in grp:
                    print("Removing DIST_MIN_V%s from %s" % (coord, opt))
                    mini = grp["DIST_MIN_V"+coord][:] * helium_mass
                    del grp["DIST_MIN_V"+coord]
                    if not "DIST_MIN_P"+coord:
                        print("Adding DIST_MIN_P%s to %s" % (coord, opt))
                        grp.create_dataset("DIST_MIN_P"+coord, (1,), data=mini,
                                           dtype='f8')
                if "DIST_MAX_V"+coord in grp:
                    print("Removing DIST_MAX_V%s from %s" % (coord, opt))
                    maxi = grp["DIST_MAX_V"+coord][:] * helium_mass
                    del grp["DIST_MAX_V"+coord]
                    if not "DIST_MAX_P"+coord:
                        print("Adding DIST_MAX_P%s to %s" % (coord, opt))
                        grp.create_dataset("DIST_MAX_P"+coord, (1,), data=maxi,
                                           dtype='f8')

def _convert1to2(fn):
    """Update version 1 HDF5 to version 2.

    - Adds dummy NBI input to existing runs.

    (Later it was decided not to include NBI input to existing runs
    so this does nothing.)
    """
    pass

def _loopchild(h5, parent):
    """Return all children or empty list if parent doesn't exist.

    The purpose of this function is to get replace following

    .. code-block:: python

       if "bfield" in h5:
           for b in h5["bfield"]:
               ...

    with

    .. code-block:: python

       for b in _loopchild(h5, "bfield"):
           ...

    """
    if parent not in h5:
        return []
    return h5[parent]
