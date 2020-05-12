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
import unyt
import random
import datetime

## Names of input parent groups.
INPUT_PARENTS = ["options", "bfield", "efield", "marker", "plasma", "neutral",
                 "wall", "boozer", "mhd", "nbi"]

## Current version
VERSION = "0"

def set_active(f, group):
    """
    Set given group as active.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group object.

    Raise:
        ValueError if the group or its parent does not exist.
    """

    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    qid = get_qid(group)
    group.parent.attrs["active"] = np.string_(qid)

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

    if str(parent) == parent:
        if parent in f:
            parent = f[parent]
        else:
            raise ValueError("Parent " + parent + " does not exist.")

    qid = parent.attrs["active"].decode('utf-8')
    group = get_group(f, qid)

    if group == None:
        raise ValueError("Active group does not exist.")
    else:
        return group

def set_activeqid(f, qid):
    """
    Set a group with the given QID as the active one.

    Args:
        f: h5py file.
        qid: Group's QID.

    Raise:
        ValueError if a group with given QID does not exist.
    """

    group = get_group(f, qid)
    set_active(f, group)

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

    group = get_active(f, parent)
    return group.parent.attrs["active"].decode('utf-8')


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
    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    group.attrs["description"] = np.string_(desc)

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

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    return group.attrs["description"].decode('utf-8')

def _set_date(f, group, date):
    """
    Set group date.

    Note that this function should only be called when the group is created.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.
        date: Date as a string in a YYYY-MM-DD hh:mm:ss format.

    Raise:
        ValueError if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    group.attrs["date"] = np.string_(date)

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

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    return group.attrs["date"].decode('utf-8')

def _set_version(f, group, date):
    """
    Set group version.

    Note that this function should only be called when the group is created.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.
        vers: Version number as a string.

    Raise:
        ValueError if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    group.attrs["version"] = np.string_(vers)

def get_version(f, group):
    """
    Get input version.

    Args:
        f: h5py file.
        group: Either the group's name or its h5py group.

    Returns:
        Version as a string.

    Raise:
        ValueError if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group" + group)
        else:
            group = grp

    return group.attrs["version"].decode('utf-8')

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
    if(str(group) != group):
        group = group.name

    if len(group) > 10:
        qid = group[-10:]
        if qid.isdigit():
            # Seems like a valid QID
            return qid

    # Not a valid QID
    raise ValueError(group + " is not a valid QID.")

def get_type(group):
    """
    Get type from a given group or from its name.

    Args:
        group: Either the group's name or its h5py group.

    Returns:
        Type string.

    Raise:
        ValueError if group does not have a valid type.
    """
    if(str(group) != group):
        group = group.name

    group = group.split("/")[-1]
    if len(group) > 12:
        type_ = group[:-11]
        if not type_.isdigit():
            # Seems like a valid QID
            return type_

    # Not a valid type
    raise ValueError(group + " does not contain a valid type.")

def get_group(f, qid):
    """
    Scan the file and return the group the QID corresponds to.

    Args:
        f: h5py file.
        qid: QID string.

    Returns:
        h5py group or None if the group was not present.
    """

    for parent in f.keys():
        for group in f[parent].keys():
            if qid in group:
                return f[parent][group]

    return None

def get_qids(f, parent):
    """
    Get QIDs of all the groups that a given parent contains.

    The QIDs are returned as list of no specific order.

    Args:
        f: h5py file
        parent Either the parent's name or its h5py group.
    Returns:
        A list of QID strings.
    Raise:
        ValueError if parent group does not exist.
    """

    # Check the parent exists and access it
    if(str(parent) == parent):
        if not parent in f:
            raise ValueError("Could not find parent " + parent)
        parent = f[parent]

    qids = []
    for group in parent:
        qids.append(get_qid(group))

    return qids

def get_inputqids(f, rungroup, ignore=[]):
    """
    Get all QIDs that tell which input was used in the given run group.

    The QIDs are returned as list, where the order is same as in
    ascot5file.INPUT_PARENTS.

    Args:
        f: h5py file.
        rungroup: Either the run group's name or its h5py group.
        ignore: List of names of inputs that are not relevant for this run.
    Returns:
        A list of QID strings.
    Raise:
        ValueError if run group does not exist.
    """

    # Check the group exists and access it.
    if(str(rungroup) == rungroup):
        qid = get_qid(rungroup)
        grp = get_group(f, qid)
        if grp is None:
            raise ValueError("Could not find group " + rungroup)
        else:
            rungroup = grp

    qids = [];
    for inp in range(0, len(INPUT_PARENTS)):
        if INPUT_PARENTS[inp] in ignore:
            continue

        try:
            qid = rungroup.attrs["qid_" + INPUT_PARENTS[inp]].decode('utf-8')
        except KeyError as err:
            print(err)
        else:
            qids.append(qid)

    return qids

def add_group(f, parent, group, desc=None):
    """
    Create a new group. A parent is created if need be.

    If parent is created, then this new group is set as active.

    Args:
        f: h5py file.
        parent: Either the parent's name or its h5py group.
        group: Name of the group (without QID of course).
        desc: optional, Description for the group.

    Returns:
        h5py group of the new group.
    """

    # Create a parent if one does not exists yet.
    if str(parent) == parent:
        parent = f.require_group(parent)

    # Generate metadata and include qid in group's name.
    qid, date, defdesc, vers = _generate_meta()
    group = group + "_" + qid
    if desc == None:
        desc = defdesc

    # Create group and set it as the active group if it is the only group.
    # Also set date and description.
    parent.create_group(group)
    set_desc(f, group, desc)
    _set_date(f, group, date)

    if len(parent.keys()) == 1:
        set_active(f, group)

    return get_group(f, qid)


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
        ValueError if group could not be found.

    """

    # Check the group exists and access it.
    if(str(group) == group):
        try:
            #group is not a parent group
            qid = get_qid(group)
            grp = get_group(f, qid)
            if grp is None:
                raise ValueError("Could not find group " + group)
            else:
                group = grp

        except ValueError:
            #group is a parent group
            group = f[group]

    # Remove the group
    parent = group.parent
    if parent.name!='/':
        was_active = get_active(f, parent) == group
        del f[group.name]
    else:
        was_active=False
        del f[group.name]

    # Set next group active (if removed group was) or remove the parent if no
    # other groups exist
    if was_active:
        if len(parent.keys()) == 0:
            del f[parent.name]
        else:
            date = "0"
            for grp in parent:
                grpdate = get_date(f, grp)
                if grpdate > date:
                    date  = grpdate
                    group = grp

            set_active(f, group)


def copy_group(fs, ft, group, newgroup=False):
    """
    Copy group from one file to another. A parent is also created if need be.

    The new group is set as active if the parent on the target file has no other
    groups. The copied group retains its original QID and date of creation.

    Args:
        fs: h5py file from which group is copied.
        ft: h5py file to which group is copied.
        group: Either the group's name or its h5py group.
        newgroup : bool, optional Flag indicating if copied group should be
        given new QID and creation date. Default is false.

    Returns:
        h5py group of the copied group.

    Raise:
        ValueError if group cannot be found or if it exists on the target file.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(fs, qid)
        if grp is None:
            raise ValueError("Could not find group " + group)
        else:
            group = grp

    # Create target parent if none exists.
    parentname = group.parent.name

    newparent  = ft.require_group(parentname)
    if group.name in newparent:
        raise ValueError("Target already has the group " + group.name)

    # Copy
    if newgroup:
        qid, date, defdesc, vers = _generate_meta()
        newname = group.name[:-11] + "_" + qid
        fs.copy(group, newparent, name=newname)
        _set_date(ft, newname, date)
        newgroupobj = ft[parentname][newname]
    else:
        fs.copy(group, newparent)
        newgroupobj = ft[parentname][group.name]

    if(len(newparent)) == 1:
        set_active(ft, newgroupobj)

    return newgroupobj


def write_data(group, name, data, dtype="f8", unit=None):
    """
    Write a dataset.

    The shape of the written dataset is deduced from the given dataset.

    Args:
        group : HDF5 group where the dataset will be written.
        name : Name of the new dataset.
        data : Data to be written.
        dtype : Data type.
        unit : Unit string if the data has units.
    """
    g = group.create_dataset(
        name  = name,
        shape = name.shape,
        data  = data,
        dtype = dtype
    )
    if unit is not None:
        g.attrs.create("unit", np.string_(unit))


def read_data(group, name):
    """
    Read a dataset.
    """
    if "unit" in group[name].attrs.keys():
        unit_str = group[name].attrs["unit"]
        unit     = unyt.Unit(unit_str)

        return group[name][:] * unit
    else:
        return group[name][:]


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
    date = date[0:19]

    desc = "No description."

    vers = VERSION

    return (qid, date, desc, vers)
