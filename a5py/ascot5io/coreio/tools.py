"""Tools to modify HDF5 files.

The purpose of these tools is to work on files that are possibly corrupted,
so they can't necessarily be processed with a5py. Therefore, the methods
here use the low-level API only (with the exception for the combine tool).
"""
import numpy as np
import h5py
import warnings
from a5py import Ascot, AscotNoDataException
from . import fileapi

def call_fileapi(fn, method, *args):
    """Wrapper for calling fileapi methods.

    This wrapper handles opening and closing of HDF5 file.

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    method : str
        Name of the method to be called.
    *args
        Arguments the called method requires.

    Returns
    -------
    val
        The value the called method returns.

    Raises
    ------
    ValueError
        If the method does not exist.
    """
    if hasattr(fileapi, method):
        method_to_call = getattr(fileapi, method)
        with h5py.File(fn, "a") as f:
            return method_to_call(f, *args)
    else:
        raise ValueError(method + " is not a valid method.")

def removegroup(fn, group, force=False):
    """Remove a group or a parent.

    Tries to remove a group. If the group to be removed is an input group which
    has been used as an input in some of the run groups, the group won't be
    removed and exception is raised.

    Note that to reclaim the disk space which the group occupied, one needs
    to call h5repack in a terminal.

    Parameters
    ----------
    fn : str
        Full path to HDF5 file.
    group: str
        Name of the group to be removed.
    force: bool, optional
        Remove group and don't raise an exception even if the group has been
        used as an input.

    Raises
    ------
    RuntimeError
        If the group has been used as an input and ```force`` = False.
    """
    with h5py.File(fn, "a") as f:
        isparent = False
        try:
            # This will raise an exception if group is not a data group.
            qid = fileapi.get_qid(group)
            filegroup = fileapi.get_group(f, qid)

            # This is a data group. If it is a run group or removal is forced,
            # we can remove it directly.
            if force or fileapi.get_type(filegroup) == "run":
                fileapi.remove_group(f, group)
                return

            # For input data groups we have to check if any run group refers to
            # it.
            runqids = fileapi.get_qids(f, "results")

            for runqid in runqids:
                rungroup = fileapi.get_group(f, runqid)
                inqids   = fileapi.get_inputqids(f, rungroup)
                if qid in inqids:
                    raise RuntimeError("Run " + runqid
                                       + " has used group " + qid
                                       + " as an input. Removal aborted.")

            # No references, the group can be removed.
            fileapi.remove_group(f, group)

        except AscotNoDataException:
            isparent = True

        if isparent:
            # The group is a parent group. If it is a results group or removal
            # is forced, we can remove it directly.
            if force or group == "results":
                fileapi.remove_group(f, group)
                return

            try:
                runqids  = fileapi.get_qids(f, "results")
            except AscotNoDataException:
                # There is no results group, so we can safely remove the group
                pass
            else:
                # For input parent groups we have to check if any run group
                # refers to any of its data groups.
                dataqids = fileapi.get_qids(f, group)

                for dataqid in dataqids:
                    for runqid in runqids:
                        rungroup = fileapi.get_group(f, runqid)
                        inqids   = fileapi.get_inputqids(f, rungroup)
                        if dataqid in inqids:
                            raise RuntimeError(
                                "Run %s has used group %s as an input. Removal"\
                                + " aborted." % (runqid, dataqid))

            # No references, the group can be removed.
            fileapi.remove_group(f, group)

def copygroup(fns, fnt, group, newgroup=False):
    """Copy a group or a parent to a different HDF5 file.

    The copied field is set as the active field in target HDF5 file. The copied
    group retains it qid and description.

    Parameters
    ----------
    fns : str
        Full path to source file.
    fnt : str
        Full path to target file.
    group : str
        Name of the group to be copied.
    newgroup : bool, optional
        Flag indicating if copied group should be given new QID and creation
        date.

    Returns
    -------
    name : str
        Name of the copied group or None if the group was a parent.
    """

    # Get the target field and type from source
    try:
        with h5py.File(fnt, "r") as ft:
            pass
    except FileNotFoundError:
        warnings.warn("Creating file %s" % fnt)
        ft = h5py.File(fnt, "w")
        ft.close()

    with h5py.File(fns, "r") as fs, h5py.File(fnt, "a") as ft:

        try:
            # This will raise exception if group is not a data group
            fileapi.get_qid(group)

            # This is a data group
            grp = fileapi.copy_group(fs, ft, group, newgroup=newgroup)
            return grp.name
        except AscotNoDataException:
            # This is a parent group
            qids = fileapi.get_qids(fs, group)
            for qid in qids:
                grp = fileapi.get_group(fs, qid)
                fileapi.copy_group(fs, ft, grp, newgroup=newgroup)
            return None

def combineoutput(fnt, addfns=None, contfns=None):
    """Combine outputs of two HDF5 files.

    Depending on which argument is given, this either combines the output by
    assuming ``fnt`` and ``addfns`` are both part of a larger simulation,
    or markers in ``contfns`` were initialized from the endstate of ``fnt``.

    The combined outputs will be stored on ``fnt`` on a run group that is
    active. The outputs are read from ``addfns`` (or ``contfns``) from the group
    that is active.

    Parameters
    ----------
    fnt : str
        Full path to HDF5 file where combined output will be stored.
    addfns : str, optional
        Name of the HDF5 file where results are read when assuming that they are
        all part of a large simulation.
    contfns : str, optional
        Name of the HDF5 file where results are read when assuming that they
        continue the simulation from the endstate of active run in ``fnt``.
    """
    if addfns is not None:
        fns = addfns
    elif contfns is not None:
        fns = contfns
    else:
        raise ValueError("Specify either addfns or contfns but not both.")

    # Find the active groups
    source = Ascot(fns).data.active
    target = Ascot(fnt).data.active

    # Combine states
    if addfns is not None:
        if hasattr(target, "inistate") and hasattr(source, "inistate"):
            with target.inistate as tdata, source.inistate as sdata:
                for field in tdata:
                    tsize = tdata[field].size
                    ssize = sdata[field].size

                    tdata[field].resize((tsize+ssize, ))
                    tdata[field][tsize:] = sdata[field][:]

        if hasattr(target, "endstate") and hasattr(source, "endstate"):
            with target.endstate as tdata, source.endstate as sdata:
                for field in tdata:
                    tsize = tdata[field].size
                    ssize = sdata[field].size

                    tdata[field].resize((tsize+ssize, ))
                    tdata[field][tsize:] = sdata[field][:]
    else:
        # Sort target inistate by id.
        if hasattr(target, "inistate"):
            with target.inistate as tdata:
                idx = np.argsort(tdata["id"][:])
                for field in tdata:
                    tdata[field][:] = tdata[field][:][idx]

        # Replace target endstate with sorted by id source endstate.
        if hasattr(target, "endstate") and hasattr(source, "endstate"):
            with target.endstate as tdata, source.endstate as sdata:
                idx = np.argsort(sdata["id"][:])
                for field in tdata:
                    tdata[field][:] = sdata[field][:][idx]

    # Combine dists
    for d in ["5d", "6d", "rho5d", "rho6d", "com"]:
        dname = "dist" + d
        if hasattr(target, d) and \
           hasattr(source, d):
            with getattr(target, dname) as tdata, \
                 getattr(source, dname) as sdata:
                tdata["ordinate"][:] += sdata["ordinate"][:]

    # Combine orbits
    if hasattr(target, "orbit") and hasattr(source, "orbit"):
        with target.orbit as tdata, source.orbit as sdata:
            for field in tdata:
                tsize = tdata[field].size
                ssize = sdata[field].size

                tdata[field].resize((tsize+ssize, ))
                tdata[field][tsize:] = sdata[field][:]

    # Combine transport coefficients
    if hasattr(target, "transcoef") and hasattr(source, "transcoef"):
        with target.transcoef as tdata, source.transcoef as sdata:
            for field in tdata:
                tsize = tdata[field].size
                ssize = sdata[field].size

                tdata[field].resize((tsize+ssize, ))
                tdata[field][tsize:] = sdata[field][:]
