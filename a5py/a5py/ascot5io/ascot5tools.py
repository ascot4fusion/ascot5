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
    Wrapper for calling ascot4file methods.

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

            # This is a data group. If it is a run group or removal is forced,
            # we can remove it directly.
            if force or ascot5file.get_type(group) is "run":
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
            if force or group is "results":
                del f[group]
                return

            # For input parent  groups we have to check if any run group refers
            # to any of its data groups.
            runqids  = ascot5file.get_qids(f, "results")
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
            del f[group]


def copygroup(fns, fnt, group):
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
    """

    # Get the target field and type from source
    with h5py.File(fns, "r") as fs, h5py.File(fnt, "a") as ft:

        try:
            # This will raise exception if group is not a data group
            ascot5file.get_qid(group)

            # This is a data group
            ascot5file.copy_group(fs, ft, group)
        except ValueError:
            # This is a parent group
            qids = ascot5file.get_qids(fs, group)
            for qid in qids:
                grp = ascot5file.get_group(fs, qid)
                ascot5file.copy_group(fs, ft, grp)


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
    if hasattr(target, "R_phi_z_vpa_vpe_t_q") and \
       hasattr(source, "R_phi_z_vpa_vpe_t_q"):
        with target.R_phi_z_vpa_vpe_t_q as tdata, \
             source.R_phi_z_vpa_vpe_t_q as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    if hasattr(target, "R_phi_z_vr_vphi_vz_t_q") and \
       hasattr(source, "R_phi_z_vr_vphi_vz_t_q"):
        with target.R_phi_z_vr_vphi_vz_t_q as tdata, \
             source.R_phi_z_vr_vphi_vz_t_q as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    if hasattr(target, "rho_pol_phi_vpa_vpe_t_q") and \
       hasattr(source, "rho_pol_phi_vpa_vpe_t_q"):
        with target.rho_pol_phi_vpa_vpe_t_q as tdata, \
             source.rho_pol_phi_vpa_vpe_t_q as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    if hasattr(target, "rho_pol_phi_vr_vphi_vz_t_q") and \
       hasattr(source, "rho_pol_phi_vr_vphi_vz_t_q"):
        with target.rho_pol_phi_vr_vphi_vz_t_q as tdata, \
             source.rho_pol_phi_vr_vphi_vz_t_q as sdata:
            tdata["ordinate"][:] += sdata["ordinate"][:]

    # Combine orbits
    if hasattr(target, "orbits") and hasattr(source, "orbits"):
        with target.orbits as tdata, source.orbits as sdata:
            for field in tdata:
                tsize = tdata[field].size
                ssize = sdata[field].size

                tdata[field].resize((tsize+ssize, ))
                tdata[field][tsize:] = sdata[field][:]
