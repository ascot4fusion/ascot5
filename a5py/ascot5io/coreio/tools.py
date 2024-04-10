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

def combineoutput(fnt, fns, add=True):
    """Combine outputs of two HDF5 files.

    Depending on ``add`` value, this either combines the output by assuming
    ``fnt`` and ``fns`` are both part of a larger simulation, or markers in
    ``fns`` were initialized from the endstate of ``fnt``.

    The combined outputs will be stored on ``fnt`` on a run group that is
    active. The outputs are read from ``addfns`` (or ``contfns``) from the group
    that is active.

    Parameters
    ----------
    fnt : str
        Full path to HDF5 file where combined output will be stored.
    fns : str, optional
        Name of the HDF5 file where results are read.
    add : bool, optional
        If True, outputs are combined assuming they are part of a larger
        simulation.

        If False, it is assumed that the source simulation continues the target
        simulation.
    """
    with h5py.File(fns, "r") as fs, h5py.File(fnt, "a") as ft:
        source = fileapi.get_active(fs, "results")
        target = fileapi.get_active(ft, "results")

        # Combine states
        if add:
            if "inistate" in target and "inistate" in source:
                tdata = target["inistate"]
                sdata = source["inistate"]
                for field in tdata:
                    tsize = tdata[field].size
                    ssize = sdata[field].size

                    tdata[field].resize((tsize+ssize, ))
                    tdata[field][tsize:] = sdata[field][:]

            if "endstate" in target and "endstate" in source:
                tdata = target["endstate"]
                sdata = source["endstate"]
                for field in tdata:
                    tsize = tdata[field].size
                    ssize = sdata[field].size

                    tdata[field].resize((tsize+ssize, ))
                    tdata[field][tsize:] = sdata[field][:]
        else:
            # Sort both target and source data by id and replace
            if "endstate" in target and "endstate" in source:
                tdata = target["endstate"]
                sdata = source["endstate"]
                t_sorted = np.argsort(tdata["ids"][:])
                s_sorted = np.argsort(sdata["ids"][:])
                idx = np.argwhere(np.in1d(tdata['ids'][t_sorted],
                                          sdata['ids'][s_sorted])).ravel()
                for field in tdata:
                    data = tdata[field][:][t_sorted]
                    if field in ['time', 'mileage', 'cputime']:
                        data[idx] += sdata[field][:][s_sorted]
                    else:
                        data[idx] = sdata[field][:][s_sorted]
                    tdata[field][:] = data

        # Combine dists
        for d in ["5d", "6d", "rho5d", "rho6d", "com"]:
            dname = "dist" + d
            if dname in target and dname in source:
                tdata = target[dname]
                sdata = source[dname]
                tdata["ordinate"][:] += sdata["ordinate"][:]

        # Combine orbits
        if "orbit" in target and "orbit" in source:
            tdata = target["orbit"]
            sdata = source["orbit"]
            for field in tdata:
                tsize = tdata[field].size
                ssize = sdata[field].size

                tdata[field].resize((tsize+ssize, ))
                tdata[field][tsize:] = sdata[field][:]

        # Combine transport coefficients
        if "transcoef" in target and "transcoef" in source:
            tdata = target["transcoef"]
            sdata = source["transcoef"]
            for field in tdata:
                tsize = tdata[field].size
                ssize = sdata[field].size

                tdata[field].resize((tsize+ssize, ))
                tdata[field][tsize:] = sdata[field][:]
