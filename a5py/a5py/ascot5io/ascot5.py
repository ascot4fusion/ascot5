"""
Main module for reading ASCOT5 HDF5 files.
"""
import numpy as np
import h5py
from . import B_2D
from . import B_3D
from . import B_ST
from . import B_TC
from . import B_GS

from . import E_TC
from . import E_1D

from . import wall_2D
from . import wall_3D

from . import plasma_1D

from . import metadata

from . import markers

from . import orbits
from . import dists
from . import states

def read_hdf5(fn, groups="all"):
    """
    Read all or specified data groups that are present in ASCOT5 HDF5 file.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
        Full path to the HDF5 file to be read.
    groups: str list, optional
        List of groups to be read. Default is all.

    Returns
    -------

    Dictionary filled with data. Structure is similar as
    the HDF5 file.
    """

    if groups == "all":
        groups = ["bfield", "efield", "options", "wall", "plasma"
                  "markers", "metadata", "states", "orbits", "dists"]

    f = h5py.File(fn, "r")

    # Read the requested input if present.
    out = {}

    out["options"] = {}
    if "options" in f and "options" in groups:
        out["options"] = options.read_hdf5(f["options"])

    out["bfield"] = {}
    if "bfield" in f and "bfield" in groups:
        if "B_2D" in f["bfield"]:
            out["bfield"]["B_2D"] = B_2D.read_hdf5(fn)
        if "B_3D" in f["bfield"]:
            out["bfield"]["B_3D"] = B_3D.read_hdf5(fn)
        if "B_GS" in f["bfield"]:
            out["bfield"]["B_GS"] = B_GS.read_hdf5(fn)
        if "B_ST" in f["bfield"]:
            out["bfield"]["B_ST"] = B_ST.read_hdf5(fn)
        if "B_TC" in f["bfield"]:
            out["bfield"]["B_TC"] = B_TC.read_hdf5(fn)

    out["efield"] = {}
    if "efield" in f and "efield" in groups:
        if "E_1D" in f["efield"]:
            out["efield"]["E_1D"] = E_1D.read_hdf5(fn)
        if "E_TC" in f["efield"]:
            out["efield"]["E_TC"] = E_TC.read_hdf5(fn)

    out["wall"] = {}
    if "wall" in f and "wall" in groups:
        if "2D" in f["wall"]:
            out["wall"]["2D"] = wall_2D.read_hdf5(fn)
        if "3D" in f["wall"]:
            out["wall"]["3D"] = wall_3D.read_hdf5(fn)

    out["plasma"] = {}
    if "plasma" in f and "plasma" in groups:
        if "P_1D" in f["plasma"]:
            out["plasma"] = plasma_1D.read_hdf5(fn)

    out["metadata"] = {}
    if "metadata" in f and "metadata" in groups:
        out["metadata"] = metadata.read_hdf5(fn)

    out["markers"] = {}
    if "markers" in f and "markers" in groups:
        out["markers"] = markers.read_hdf5(fn)

    out["dists"] = {}
    if "distributions" in f and "dists" in groups:
        out["dists"] = dists.read_hdf5(fn)
        
    out["orbits"] = {}
    if "orbits" in f and "orbits" in groups:
        out["orbits"] = orbits.read_hdf5(fn)
        
    out["states"] = {}
    if "states" in f and "states" in groups:
        out["states"] = states.read_hdf5(fn)
        

    f.close()

    return out


def write_hdf5(fn, a5):
    """
    Generate ASCOT5 HDF5 file.

    TODO implement

    Parameters
    ----------
    
    fn : str
        Full path to HDF5 file.
    a5 : dictionary
        ASCOT5 HDF5 file in dictionary format (as given by the 
        read_hdf5 function).
    """

    f = h5py.File(fn,"a")


    f.close()
