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
CURRENT_VERSION = 2

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


def convert_1to2(fnin, fnout):
    """
    Convert version 1 to 2.

    Changes:
        - Add dummy NBI input to each run.

    **** DEPRECATED ****
    """
    pass


def convert(fnin):
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
