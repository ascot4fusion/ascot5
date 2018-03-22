"""
Initialize groups for storing ASCOT5 data.
"""
import numpy as np
import h5py
import random
import datetime

def setactive(hdf5file, path):
    """
    Set given group as active.

    Parameters
    ----------
    df5file : HDF5 file
        Target HDF5 file.
    path : HDF5 group
        Path to the group.
    """
    mastergroup = hdf5file[path].parent.name
    qid = path[-10:]
    hdf5file[mastergroup].attrs["active"] = np.string_(qid)

def generatemetadata():
    """
    Generate metadata
    
    Metdata consists of unique (not quaranteed but almost 
    certainly) id and date when field was created.

    """
    # Generate random unsigned 32 bit integer and convert it to string.
    qid = str(np.uint32(random.getrandbits(32)))

    # Left-padding with zeroes so that qid is always 10 characters long.
    while len(qid) < 10:
        qid = "0" + qid

    # Get date
    date = str(datetime.datetime.now())

    # Last digits are milliseconds which we don't need
    date = date[0:20]
    
    return (qid, date)

    
def removegroup(hdf5file, path):
    """
    Remove a subgroup.

    Parameters
    ----------
    
    hdf5file : HDF5 file
        Target HDF5 file.
    path : str
        Path to group e.g. "bfield/B_2D-XXXXXXXXXX".
    """

    if path in hdf5file:
        del hdf5file[path]

def copygroup(hdf5fs, hdf5ft, path):
    """
    Copy a subgroup and set it as active.

    Parameters
    ----------
    hdf5fs : HDF5 file
        Source HDF5 file.
    hdf5ft : HDF5 file
        Target HDF5 file.
    path : str
        Group path e.g. "bfield/B_2D-XXXXXXXXXX".
    """

    # Create target mastergroup if none exists.
    mastergroup = hdf5fs[path].parent.name
    group_id = hdf5ft.require_group(mastergroup)

    # Copy
    hdf5fs.copy(path, group_id)

    # Set as active
    qid = path[-10:]
    hdf5fs[mastergroup].attrs["active"] = np.string_(qid)

def setdescription(hdf5file, path, desc):
    """
    Set group description.

    Parameters
    ----------
    hdf5file : HDF5 file
        Target HDF5 file.
    path : str
        Path to group e.g. "bfield/B_2D-XXXXXXXXXX".
    desc : str
        Description.
    """
    hdf5file[path].attrs["description"] = np.string_(desc)   


def creategroup(hdf5file, mastergroup, subgroup, desc="-"):
    """
    Create a new subgroup.

    Creates a new subgroup (e.g. B_2D). Associated mastergroup (e.g. bfield)
    is created if one does not exists. The new group is set as being active.

    The new group will be named with unique id and metadata includes date
    when the group was created.

    Parameters
    ----------
    
    hdf5file : HDF5 file
        Target HDF5 file from h5py.File().
    mastergroup : str
        Name of the mastergroup.
    subgroup : str
        Name of the subgroup.
    desc : str, optional
        Description or other documentation related to this subgroup.

    Returns
    -------

    str : path
        Path to new group (e.g. "bfield/B_2D-0123456789")
    """

    # Crete a mastergroup if one does not exists yet.
    if not mastergroup in hdf5file:
        hdf5file.create_group(mastergroup)

    # Generate metadata and include qid in subgroup's name.
    qid, date = generatemetadata()
    subgroup = subgroup + "-" + qid

    # Create subgroup and set it as the active group
    hdf5file[mastergroup].create_group(subgroup)
    hdf5file[mastergroup].attrs["active"] = np.string_(qid)

    # Path to new group.
    path = mastergroup + "/" + subgroup

    # Set date.
    hdf5file[path].attrs["date"] = np.string_(date)

    # Set description.
    setdescription(hdf5file, path, desc)
    
    return path
