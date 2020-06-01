"""
Make older ASCOT5 HDF5 files compatible with newer versions of ASCOT5.

Whenever a new version is released, add a new function here that converts
previous version to current version. Note that you might have to modify the HDF5
directly as this tool should be independent of ascot5io (which is subject to
changes).

Also update the CURRENT_VERSION number.

File compatibility.py
"""
import h5py
import numpy as np
import tempfile
import subprocess

## KEEP THIS UPDATED ##
CURRENT_VERSION = 3

def convert_oldtonew(fnin, origin, keeptemp=False):
    """
    Convert an old file to the current version.

    The old file remains as it is, and the converted version will be named
    fnin_{CURRENT_VERSION}.h5.

    If there are several versions between the current and the old version, the
    conversion is done incrementally e.g. 1->2, 2->3, 3->4 etc.

    ***********DEPRECATED*******************

    Args:
        fnin : str <br>
            Name of the input file.
        origin : int <br>
            Version number of the old file.
        keeptemp : bool, optional <br>
            Flag whether the intermediate files are kept.
    """

    version = origin
    fntemp = fnin
    while version < CURRENT_VERSION:
        fnout = fnin[:-3] + "_" + str(version+1) + ".h5"

        # Make the version conversion.
        funname = "convert_" + str(version) + "to" + str(version+1)
        globals()[funname](fntemp, fnout)

        if (not keeptemp) and (origin < version):
            # TODO remove temporary files
            pass

        version += 1
        fntemp  = fnout


def convert(fnin):
    """
    Update version 2 HDF5 to version 3.

    - Adds dummy boozer and mhd groups.
    - Adds flags to wall input.
    - Adds ENDCOND_MAX_MILEAGE from options.
    - Replace velocity with momentum in distribution options.
    """
    from a5py.ascot5io.ascot5file import get_qid
    from a5py.ascot5io.boozer import write_hdf5_dummy as boozer_write_hdf5_dummy
    from a5py.ascot5io.mhd    import write_hdf5_dummy as mhd_write_hdf5_dummy

    fnout = fnin[:-3] + "_" + str(CURRENT_VERSION) + ".h5"

    subprocess.call(["cp", fnin, fnout])

    print("Adding a dummy Boozer and MHD input groups.")
    boozer_qid = get_qid(boozer_write_hdf5_dummy(fnout))
    mhd_qid    = get_qid(mhd_write_hdf5_dummy(fnout))

    def wrapper():
        with h5py.File(fnout, "a") as h5:
            if not "results" in h5:
                return

            print("Adding dummy Boozer and MHD inputs for existing runs.")
            for run in h5["results"]:
                h5["results"][run].attrs["qid_mhd"] = np.string_(mhd_qid)
                h5["results"][run].attrs["qid_boozer"] = np.string_(boozer_qid)

            print("Adding flags to wall inputs.")
            walls = list(h5["wall"].keys())
            for wall in walls:
                if "wall_3D" in wall:
                    g = h5["wall"][wall]
                    nelements = int(g["nelements"][:])
                    flag = np.zeros(shape=(nelements,1), dtype=np.int)
                    g.create_dataset('flag',(nelements,1), data=flag,
                                     dtype='i4')

            print("Adding ENDCOND_MAX_MILEAGE=100 and setting distribution")
            print("momentum limits (DIST_MIN/MAX_P* ) to")
            print("(DIST_MIN/MAX_V*)*helium mass. Bin number is kept same and")
            print("the old velocity options are removed.")
            for opt in h5["options"]:
                opt = h5["options"][opt]
                opt.create_dataset("ENDCOND_MAX_MILEAGE", (1,), data=100,
                                   dtype='f8')

                for coord in ["PA", "PE", "R", "PHI", "Z"]:
                    helium_mass = 6.65e-27
                    nbin = opt["DIST_NBIN_V"+coord][:]
                    mini = opt["DIST_MIN_V"+coord][:] * helium_mass
                    maxi = opt["DIST_MAX_V"+coord][:] * helium_mass
                    opt.create_dataset("DIST_NBIN_P"+coord, (1,), data=nbin,
                                       dtype='f8')
                    opt.create_dataset("DIST_MIN_P"+coord, (1,), data=mini,
                                       dtype='f8')
                    opt.create_dataset("DIST_MAX_P"+coord, (1,), data=maxi,
                                       dtype='f8')

                    del opt["DIST_NBIN_V"+coord]
                    del opt["DIST_MIN_V"+coord]
                    del opt["DIST_MAX_V"+coord]

    wrapper()

    print("Conversion complete.")


def convert_1to2(fnin):
    """
    Update version 1 HDF5 to version 2.

    - Adds dummy NBI input to existing runs.
    """
    from a5py.ascot5io.ascot5file import get_qid
    from a5py.ascot5io.nbi import write_hdf5_dummy

    fnout = fnin[:-3] + "_" + str(CURRENT_VERSION) + ".h5"

    subprocess.call(["cp", fnin, fnout])

    print("Adding a dummy NBI input group.")
    qid = get_qid(write_hdf5_dummy(fnout))

    print("Adding the dummy NBI input as a pseudo-input for existing runs.")
    def wrapper():
        with h5py.File(fnout, "a") as h5:
            if not "results" in h5:
                return

            for run in h5["results"]:
                h5["results"][run].attrs["qid_nbi"] = np.string_(qid)

    wrapper()

    print("Conversion complete.")
