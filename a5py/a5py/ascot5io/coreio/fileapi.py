"""Module for creating and modifying the HDF5 file and for accessing meta data.

ASCOT5 HDF5 file is built as follows. At top level are parent groups, one (and
only one) for each type of input (e.g. bfield for magnetic field input) and one
parent group that contains all the results.

Parent groups contain groups that hold the actual datasets. Single parent can
have multiple groups e.g. bfield can store multiple magnetic field inputs.
The result group is a different in a sense that it contains run groups, one
for each simulation, which then contain the groups that store inistate,
endstate, and whatever diagnostics were used. However, here we refer to the
top level groups as parents, and the groups that are directly below them
as data containers.

A typical structure of a file can be like this:

> /
> /bfield
> /bfield/B_2DS_1234567890
> /bfield/B_2DS_2345678901
> /bfield/B_3DS_3456789012
> ...
> <other input groups>
> ...
> /results
> /results/run_4567890123
> /results/run_5678901234

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
"""

import numpy as np
import h5py
import unyt
import random
import datetime

from collections import OrderedDict
from a5py.exceptions import AscotNoDataException, AscotIOException

INPUTGROUPS = ["options", "bfield", "efield", "marker", "plasma", "neutral",
               "wall", "boozer", "mhd", "asigma", "nbi", "marker_shined"]
"""Names of the input parent groups.
"""

OUTPUTGROUPS = ["inistate", "endstate", "dist5d", "distrho5d", "dist6d",
                "distrho6d", "orbit", "transcoef"]
"""Names of the output data containers in runs.
"""

VERSION = "5.5"
"""Current version of the code."""

def set_active(f, group):
    """Set given group as active.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.

    Raises
    ------
    AscotNoDataException
        Raised if the group or its parent does not exist.
    """

    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    qid = get_qid(group)
    group.parent.attrs["active"] = np.string_(qid)

def get_active(f, parent):
    """Get active group.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    parent : str or `h5py.Group`
        Group object or its name.

    Returns
    -------
    group : `h5py.Group`
        The active group.

    Raises
    ------
    AscotNoDataException
        Raised if the parent or the active group does not exist.
    """

    if str(parent) == parent:
        if parent in f:
            parent = f[parent]
        else:
            raise AscotNoDataException(
            "Parent " + parent + " does not exist.")

    qid = parent.attrs["active"].decode('utf-8')
    group = get_group(f, qid)

    if group == None:
        raise AscotNoDataException("Active group does not exist.")
    else:
        return group

def set_activeqid(f, qid):
    """Set a group with the given QID as the active one.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    qid : str
        QID.
    """
    group = get_group(f, qid)
    set_active(f, group)

def get_activeqid(f, parent):
    """Get QID of the currently active group.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    parent : str or `h5py.Group`
        Parent whose active group is sought.

    Returns
    -------
    qid : str
        QID string.
    """
    group = get_active(f, parent)
    return group.parent.attrs["active"].decode('utf-8')


def set_desc(f, group, desc):
    """Set group description.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.
    desc : str
        Description as a string.

    Raises
    ------
    AscotNoDataException
        Raised if group does not exist.
    """
    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    group.attrs["description"] = np.string_(desc)

def get_desc(f, group):
    """Get group description.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.

    Returns
    -------
    desc : str
        Description as a string.

    Raises
    ------
    AscotNoDataException
        Raised if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    return group.attrs["description"].decode('utf-8')

def _set_date(f, group, date):
    """Set group date.

    Note that this function should only be called when the group is created.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.
    date : str
        Date as a string in a YYYY-MM-DD hh:mm:ss format.

    Raises
    ------
    AscotNoDataException
        Raised if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    group.attrs["date"] = np.string_(date)

def get_date(f, group):
    """Get date.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.

    Returns
    -------
    date : str
        Date as a string in format YYYY-MM-DD hh:mm:ss.

    Raises
    ------
    AscotNoDataException
        Raised if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    return group.attrs["date"].decode('utf-8')

def _set_version(f, group):
    """Set group version.

    Note that this function should only be called when the group is created.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.

    Raises
    ------
    AscotNoDataException
        Raised if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    group.attrs["version"] = np.string_(vers)

def get_version(f, group):
    """Get input version.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        Group object or its name.

    Returns
    -------
    version : str
        Version as a string.

    Raises
    ------
    AscotNoDataException
        Raised if group does not exist.
    """

    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group" + group)
        else:
            group = grp

    return group.attrs["version"].decode('utf-8')

def get_qid(group):
    """Get QID from a given group.

    Parameters
    ----------
    group : str or `h5py.Group`
        Group object or its name.

    Returns
    -------
    qid : str
        QID string.

    Raises
    ------
    AscotNoDataException
        Raised if group does not have a valid QID.
    """
    if(str(group) != group):
        group = group.name

    if len(group) >= 10:
        qid = group[-10:]
        if qid.isdigit():
            # Seems like a valid QID
            return qid

    # Not a valid QID
    raise AscotNoDataException(group + " is not a valid QID.")

def get_type(group):
    """Get type from a given group or from its name.

    Parameters
    ----------
    group : str or `h5py.Group`
        Group object or its name.

    Returns
    -------
    type : str
         Group's type.

    Raises
    ------
    AscotNoDataException
        Raised if group does not have a valid type.
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
    raise AscotNoDataException(group + " does not contain a valid type.")

def get_group(f, qid):
    """Scan the file and return the group the QID corresponds to.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    qid : str
        QID string.

    Returns
    -------
    group : `h5py.Group`
        Group or None if the group was not present.
    """
    for parent in f.keys():
        for group in f[parent].keys():
            if qid in group:
                return f[parent][group]

    return None

def get_qids(f, parent):
    """Get QIDs of all the groups that a given parent contains.

    The QIDs are returned as list of no specific order.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    parent : str or `h5py.Group`
        Parent object or its name.

    Returns
    -------
    qids : list [str]
        A list of QID strings.

    Raises
    ------
    AscotNoDataException
        Raised if parent group does not exist.
    """
    # Check the parent exists and access it
    if(str(parent) == parent):
        if not parent in f:
            raise AscotNoDataException("Could not find parent " + parent)
        parent = f[parent]

    qids = []
    for group in parent:
        qids.append(get_qid(group))

    return qids

def get_inputqids(f, rungroup):
    """Get all QIDs that tell which input was used in the given run group.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    rungroup : str or `h5py.Group`
        Either the run group object or its name.

    Returns
    -------
    qids : `collections.OrderedDict` [str, str]
        Ordered dictionary with "parent name" - "qid" value pairs with the order
        being same as in `INPUTGROUPS`.

    Raises
    ------
    AscotNoDataException
        Raised if run group does not exist.
    """
    # Check the group exists and access it.
    if(str(rungroup) == rungroup):
        qid = get_qid(rungroup)
        grp = get_group(f, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group " + rungroup)
        else:
            rungroup = grp

    qids = OrderedDict();
    for inp in INPUTGROUPS:
        try:
            qids[inp] = rungroup.attrs["qid_" + inp].decode("utf-8")
        except KeyError as err:
            pass

    return qids

def add_group(f, parent, group, desc=None):
    """Create a new group. A parent is created if need be.

    If parent is created, then this new group is set as active.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    parent : str or `h5py.Group`
        Parent of the created group.
    group : str
        The group name.
    desc : str, optional
        Description for the group.

    Returns
    -------
    newgroup : `h5py.Group`
        The new group.
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
    """Remove a group.

    If this was an active group, a most recent group is set as an active
    instead. If no other groups exist the parent is also removed.

    Input groups are not removed if they have been used in an existing run.

    Note that to reclaim the disk space which the group occupied, one needs
    to call h5repack in a terminal.

    Parameters
    ----------
    f : `h5py.File`
        Open HDF5 file.
    group : str or `h5py.Group`
        The group to be removed.

    Raises
    ------
    AscotNoDataException
        Raised if group could not be found.
    AscotIOException
        Raised if the group is input used by a run.
    """
    # Check the group exists and access it.
    if(str(group) == group):
        if group in f.keys():
            # Group is a parent group
            group = f[group]
        else:
            # assume group is data group
            qid = get_qid(group)
            grp = get_group(f, qid)
            if grp is None:
                raise AscotNoDataException("Could not find group" + group)
            else:
                group = grp

    # Check if this is an input group and whether it has been used in run.
    parent = group.parent
    if "results" in f.keys() and parent.name != "results":
        if parent.name == "/":
            for run in f["results"].keys():
                for q in f["results"][run].attrs.keys():
                    if q[:4] == "qid_" and q[4:] == group.name[1:]:
                        raise AscotIOException(
                        "Input \"" + group.name[1:] + "\" is used by " + run)
        else:
            qid = get_qid(group)
            for run in f["results"].keys():
                for q in f["results"][run].attrs.keys():
                    t = f["results"][run].attrs[q].decode('utf-8')
                    if q[:4] == "qid_" and t == qid:
                        raise AscotIOException(
                            "Input " + qid + " is used by " + run)

    # Remove the group
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
    """Copy group from one file to another. A parent is also created if need be.

    The new group is set as active if the parent on the target file has no other
    groups. The copied group retains its original QID and date of creation.

    Parameters
    ----------
    fs : `h5py.File`
        File from which the group is copied.
    ft : `h5py.File`
        File to which the group is copied.
    group : str pr `h5py.Group`
        The group to be copied.
    newgroup : bool, optional
        Flag indicating if copied group should be given new QID and date.

    Returns
    -------
    groupt : `h5py.Group`
        The new group which is a copy of group at ft.

    Raises
    ------
    AscotNoDataException
        Raised if the copied group cannot be found.
    AscotIOException
        Raised if the group already exists on the target file.
    """
    # Check the group exists and access it.
    if(str(group) == group):
        qid = get_qid(group)
        grp = get_group(fs, qid)
        if grp is None:
            raise AscotNoDataException("Could not find group " + group)
        else:
            group = grp

    # Create target parent if none exists.
    parentname = group.parent.name

    newparent  = ft.require_group(parentname)
    if group.name in newparent:
        raise AscotIOException("Target already has the group " + group.name)

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


def write_data(group, name, data, shape, dtype, unit=None):
    """Write a dataset.

    The shape of the written dataset is same as the input array.

    Parameters
    ----------
    group : `h5py.Group`
        HDF5 group where the dataset will be written.
    name : str
        Name of the new dataset.
    data : `np.array`
        Data to be written.
    dtype : str
        Data type.
    unit : str, optional
        Unit string if the data has units.
    """
    g = group.create_dataset(
        name  = name,
        shape = shape,
        data  = data,
        dtype = dtype
    )
    if unit is not None:
        g.attrs.create("unit", np.string_(unit))


def read_data(group, name):
    """Read a dataset and add units if present.

    Parameters
    ----------
    group : `h5py.Group`
        HDF5 group containing the dataset.
    name : str
        Name of the dataset.

    Returns
    -------
    data : `np.array` or `unyt.array`, `(N,)`
        Dataset with units read from the dataset attribute "unit".

        If dataset has no "unit" attribute, ordinary `np.array` is returned.
    """
    if "unit" in group[name].attrs.keys():
        unit_str = group[name].attrs["unit"]
        unit     = unyt.Unit(unit_str)

        return group[name][:].ravel() * unit
    else:
        return group[name][:].ravel()



def _generate_meta():
    """Generate QID, date and default description/tag.

    Calls random number generator to create 32 bit string which is then
    converted as a QID string using left-padding with zeroes if necessary.

    Returns
    -------
    meta : `tuple` [str, str, str, str]
        QID, date, description, and version number.
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

    desc = "TAG"

    vers = VERSION

    return (qid, date, desc, vers)
