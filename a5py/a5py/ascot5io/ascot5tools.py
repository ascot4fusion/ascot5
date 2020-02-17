"""
Tools to modify HDF5 files.

File: ascot5tools.py
"""
import numpy as np
import h5py
from . import ascot5
from . import ascot5file

def call_ascot5file(fn, method, *args):
    """
    Wrapper for calling ascot5file methods.

    This wrapper handles opening and closing of HDF5 file.

    Args:
        fn : str <br>
            Full path to HDF5 file.
        method: str <br>
            Name of the method to be called.
        *args: <br>
            Arguments the called method requires.

    Returns:
        The value the called method returns.

    Raise:
        ValueError if the method does not exist.
    """
    if hasattr(ascot5file, method):
        method_to_call = getattr(ascot5file, method)
        with h5py.File(fn, "a") as f:
            return method_to_call(f, *args)
    else:
        raise ValueError(method + " is not a valid method.")


def removegroup(fn, group, force=False):
    """
    Remove a group or a parent.

    Tries to remove a group. If the group to be removed is an input group which
    has been used as an input in some of the run groups, the group won't be
    removed and exception is raised.

    Note that to reclaim the disk space which the group occupied, one needs
    to call h5repack in a terminal.

    Args:
        fn : str <br>
            Full path to HDF5 file.
        group: str <br>
            Name of the group to be removed.
        force: bool, optional <br>
            Remove group and don't raise an exception even if the group has been
            used as an input. Default is false.
    Raise:
        RuntimeError if the group has been used as an input.
    """
    with h5py.File(fn, "a") as f:

        try:
            # This will raise an exception if group is not a data group.
            qid = ascot5file.get_qid(group)
            filegroup = ascot5file.get_group(f, qid)

            # This is a data group. If it is a run group or removal is forced,
            # we can remove it directly.
            if force or ascot5file.get_type(filegroup) == "run":
                ascot5file.remove_group(f, group)
                return

            # For input data groups we have to check if any run group refers to
            # it.
            runqids = ascot5file.get_qids(f, "results")

            for runqid in runqids:
                rungroup = ascot5file.get_group(f, runqid)
                inqids   = ascot5file.get_inputqids(f, rungroup)
                if qid in inqids:
                    raise RuntimeError("Run " + runqid
                                       + " has used group " + qid
                                       + " as an input. Removal aborted.")

            # No references, the group can be removed.
            ascot5file.remove_group(f, group)

        except ValueError:
            # The group is a parent group. If it is a results group or removal
            # is forced, we can remove it directly.
            if force or group == "results":
                ascot5file.remove_group(f, group)
                return

            try:
                runqids  = ascot5file.get_qids(f, "results")
            except ValueError:
                # There is no results group, so we can safely remove the group
                pass
            else:
                # For input parent groups we have to check if any run group
                # refers to any of its data groups.
                dataqids = ascot5file.get_qids(f, group)

                for dataqid in dataqids:
                    for runqid in runqids:
                        rungroup = ascot5file.get_group(f, runqid)
                        inqids   = ascot5file.get_inputqids(f, rungroup)
                        if dataqid in inqids:
                            raise RuntimeError("Run " + runqid
                                               + " has used group " + dataqid
                                               + " as an input. Removal aborted.")

            # No references, the group can be removed.
            ascot5file.remove_group(f, group)


def copygroup(fns, fnt, group, newgroup=False):
    """
    Copy a group or a parent to a different HDF5 file.

    The copied field is set as the active field in target HDF5 file. The copied
    group retains it qid and description.

    Args:
        fns : str <br>
            Full path to source file.
        fnt : str <br>
            Full path to target file.
        group : str <br>
            Name of the group to be copied.
        newgroup : bool, optional <br>
            Flag indicating if copied group should be given new QID and creation
            date. Default is false.

    Returns:
        Name of the copied group or none if the group was parent.
    """

    # Get the target field and type from source
    with h5py.File(fns, "r") as fs, h5py.File(fnt, "a") as ft:

        try:
            # This will raise exception if group is not a data group
            ascot5file.get_qid(group)

            # This is a data group
            grp = ascot5file.copy_group(fs, ft, group, newgroup=newgroup)
            return grp.name
        except ValueError:
            # This is a parent group
            qids = ascot5file.get_qids(fs, group)
            for qid in qids:
                grp = ascot5file.get_group(fs, qid)
                ascot5file.copy_group(fs, ft, grp, newgroup=newgroup)


def combineoutput(fnt, addfns=None, contfns=None):
    """
    Combine outputs of two HDF5 files.

    Depending on which argument is given, this either combines the output by
    assuming fnt and addfns are both part of a larger simulation ("add"), or
    markers in contfns were initialized from the endstate of fnt ("continue").

    The combined outputs will be stored on fnt on a run group that is active.
    The outputs are read from addfns (or contfns) from the group that is active.

    Args:
        fnt : str <br>
            Full path to HDF5 file where combined output will be stored.
        addfns : str, optional <br>
            Name of the HDF5 file where results are read if "add" mode is being
            used.
        contfns : str, optional <br>
            Name of the HDF5 file where results are read if "continue" mode is
            being used.
    """

    if addfns is not None:
        fns = addfns
    elif contfns is not None:
        fns = contfns
    else:
        return

    # Find the active groups
    source = (ascot5.Ascot(fns)).active
    target = (ascot5.Ascot(fnt)).active

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
    if hasattr(target, "dist5d") and \
       hasattr(source, "dist5d"):
        with target.dist5d as tdata, \
             source.dist5d as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    if hasattr(target, "dist6d") and \
       hasattr(source, "dist6d"):
        with target.dist6d as tdata, \
             source.dist6d as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    if hasattr(target, "distrho5d") and \
       hasattr(source, "distrho5d"):
        with target.distrho5d as tdata, \
             source.distrho5d as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    if hasattr(target, "distrho6d") and \
       hasattr(source, "distrho6d"):
        with target.distrho6d as tdata, \
             source.distrho6d as sdata:
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
