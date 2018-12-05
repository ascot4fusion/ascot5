"""
Module for creating and modifying the HDF5 file and for accessing meta data.

ASCOT5 HDF5 file is built as follows. At top level are parent groups, one (and
only one) for each type of input (e.g. bfield for magnetic field input) and one
parent group that contains all the results.

Parent groups contain groups that hold the actual datasets. Single parent can
have multiple groups e.g. bfield can store multiple magnetic field inputs.
The result group is a different in a sense that it contains run groups, one
for each simulation, which then contain the groups that store inistate,
endstate, and whatever diagnostics were used. However, here we refer to the
top level groups as parents, and the groups that are directly below them
as groups.

A typical structure of a file can be like this:

> /                              <br>
> /bfield                        <br>
> /bfield/B_2DS-1234567890       <br>
> /bfield/B_2DS-2345678901       <br>
> /bfield/B_3DS-3456789012       <br>
> /efield                        <br>
> /efield/E_TC-1234567890        <br>
> /neutral                       <br>
> /neutral/N0_3D-1234567890      <br>
> /plasma                        <br>
> /plasma/plasma_1D-1234567890   <br>
> /wall                          <br>
> /wall/wall_2D-1234567890       <br>
> /options                       <br>
> /options/opt-1234567890        <br>
> /marker                        <br>
> /marker/particle-1234567890    <br>
> /results                       <br>
> /results/run-1234567890        <br>
> /results/run-2345678901        <br>

The number in group name is unique identifier (QID) which is generated when a
group is created. It is generated from a 32 bit random integer and as such it
is almost quaranteed to be unique.

QID is used to determine which input group is active, i.e., to be used in
a simulation, and to which run group the simulation results are written.
QID of the active group is stored as an attribute (in a string format) in the
parent. For the results parent this is the QID of the most recent run (unless
manually changed). Run groups also store QIDs of all input fields used in that
particular run.

All groups also store date and time at which they where created and also a
description field that user can modify at will (see ascot5.py).

The groups should be created, removed, and their metadata accessed via this
module. The actual datasets are accessed via corresponding modules. This module
assumes the HDF5 file is open when these routines are called, and closed
afterwards.

File: ascot5file.py
"""

import numpy as np
import h5py
import random
import datetime

def set_active(f, group):
    """
    Set given group as active.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group object.

    Raise:
        ValueError if the group or its parent does not exist.
    """
    mastergroup = hdf5file[path].parent.name
    qid = path[-10:]
    hdf5file[mastergroup].attrs["active"] = np.string_(qid)

    raise ValueError("Requested group does not exist.")

def get_active(f, parent):
    """
    Get active group.

    Args:
        f: h5py file.
        parent: Either parent's name or its h5py group.

    Returns:
        h5py group object.

    Raise:
        ValueError if the parent or the active group does not exist.
    """

    raise ValueError("Requested group does not exist.")

def set_activeqid(f, qid):
    """
    Set a group with the given QID as the active one.

    Args:
        f: h5py file.
        parent: Parent's name.

    Raise:
        ValueError if a group with given QID does not exist.
    """
    hdf5file[mastergroup].attrs["active"] = np.string_(qid)

    raise ValueError("Requested group does not exist.")

def get_activeqid(f, parent):
    """
    Get QID of the group which is the active one.

    Args:
        f: h5py file.
        parent: Parent's name.

    Returns:
        QID string.

    Raise:
        ValueError if parent or active group does not exist.
    """

    raise ValueError("Requested group does not exist.")


def set_desc(f, group, desc):
    """
    Set group description.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.
        desc: Description as a string.

    Raise:
        ValueError if group does not exist.
    """
    hdf5file[path].attrs["description"] = np.string_(desc)

def get_desc(f, group):
    """
    Get group description.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.

    Returns:
        Date as a string.

    Raise:
        ValueError if group does not exist.
    """

def get_date(f, group):
    """
    Get date (as a string) in a YYYY-MM-DD hh:mm:ss format.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.

    Returns:
        Date as a string.

    Raise:
        ValueError if group does not exist.
    """

def get_qid(group):
    """
    Get QID from a given group or from its name.

    Args:
        group: Either the group's name or its h5py group.

    Returns:
        QID string.

    Raise:
        ValueError if group does not have a valid QID.
    """

def add_group(f, parent, group, desc=None):
    """
    Create a new group. A parent is created if need be.

    If parent is created, then this new group is set as active.

    Args:
        f: h5py file.
        parent: Either the parent's name or its h5py group.
        group: Name of the group.
        desc: optional, Description for the group.

    Raise:
        ValueError if group does not have a valid QID.
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
    if desc:
        setdescription(hdf5file, path, desc)
    else:
        setdescription(hdf5file, path, "-")

    return path

def remove_group(f, group):
    """
    Remove a group.

    If this was an active group, a most recent group is set as an active
    instead. If no other groups exist the parent is also removed.

    Note that to reclaim the disk space which the group occupied, one needs
    to call h5repack in a terminal.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.

    Raise:
        ValueError if group does not have a valid QID.
    """

    if path in hdf5file:
        del hdf5file[path]

def copy_group(fs, ft, group):
    """
    Copy group from one file to another. A parent is also created if need be.

    The new group is set as active if the parent on the target file has no other
    groups.

    Args:
        fs: h5py file from which group is copied.
        ft: h5py file to which group is copied.
        group: Either the group's name or its h5py group.

    Raise:
        ValueError if group does not have a valid QID.
    """

    # Create target mastergroup if none exists.
    mastergroup = hdf5fs[path].parent.name
    group_id = hdf5ft.require_group(mastergroup)

    # Copy
    hdf5fs.copy(path, group_id)

    # Set as active
    qid = path[-10:]
    hdf5fs[mastergroup].attrs["active"] = np.string_(qid)

def _generate_meta():
    """
    Generate QID, date and default description "No description.".

    Calls random number generator to create 32 bit string which is then
    converted as a QID string using left-padding with zeroes if necessary.

    Returns:
        tuple (QID, date, description) where all elements are strings.
    """
    # Generate random unsigned 32 bit integer and convert it to string.
    qid = str(np.uint32(random.getrandbits(32)))

    # Left-padding with zeroes so that QID is always 10 characters long.
    while len(qid) < 10:
        qid = "0" + qid

    # Get date
    date = str(datetime.datetime.now())

    # Last digits are milliseconds which we don't need
    date = date[0:20]

    desc = "No description."

    return (qid, date, desc)
