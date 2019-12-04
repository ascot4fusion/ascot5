"""
Make older ASCOT5 HDF5 files compatible with newer versions of ASCOT5.

File compatibility.py
"""
import h5py
import numpy as np
import tempfile

from . ascot5file  import add_group, get_qid, get_type, get_activeqid
from . ascot5tools import call_ascot5file
from . ascot5      import create_inputobject, create_outputobject

def convert_oldtonew(fnin, fnout, oldver):

    fnfinal = fnout
    newver = 2
    while oldver < newver:
        if oldver + 1 == newver:
            fnout = fnfinal
        else:
            fnout, _ = tempfile.mkstemp()

        globals()["convert_" + str(oldver) + "to" + str(oldver+1)](fnin, fnout)
        fnin = fnout
        oldver += 1


def convert_1to2(fnin, fnout):

    # Store old and corresponding new QIDs in a dictionary
    qids = {}

    # It might not be possible to initialize Ascot object so we access HDF5 file
    # directly.
    with h5py.File(fnin, "r") as h5in:

        # Create run groups in the target file if those exists in h5in
        if "results" in h5in:
            with h5py.File(fnout, "a") as h5out:
                h5out.create_group("results")

                for run in h5in["results"]:
                    g = add_group(h5out, h5out["results"], "run")
                    qids[get_qid(run)] = get_qid(g)

                    ## Update metadata ##
                    h5out["results"][g.name].attrs["description"] = \
                        h5in["results"][run].attrs["description"]
                    h5out["results"][g.name].attrs["date"] = \
                        h5in["results"][run].attrs["date"]

        for parent in h5in:
            if parent == "results":
                for run in h5in["results"]:
                    newrun = "run_" + qids[get_qid(run)]

                    ## Output groups ##
                    for group in h5in["results"][run]:
                        if group == "inistate" or group == "endstate":
                            o = create_outputobject(group,
                                                    h5in["results"][run][group],
                                                    None)
                            data = o.read()
                            o.write(fnout, newrun, group, data)
                        else:
                            o = create_outputobject(group,
                                                    h5in["results"][run][group],
                                                    None)
                            data = o.read()
                            o.write(fnout, newrun, data)

            else:
                ## Input groups ##
                for group in h5in[parent]:
                    o = create_inputobject(get_type(group), h5in[parent][group])
                    data = o.read()
                    g = o.write(fnout, data)
                    qids[o.get_qid()] = get_qid(g)

                    with h5py.File(fnout, "a") as h5out:
                        ## Update metadata ##
                        h5out[parent][g].attrs["description"] = \
                            h5in[parent][group].attrs["description"]
                        h5out[parent][g].attrs["date"] = \
                            h5in[parent][group].attrs["date"]


        ## Set relationships correctly ##
        for parent in h5in:
            oldqid = get_activeqid(h5in, parent)
            newqid = qids[oldqid]
            call_ascot5file(fnout, "set_activeqid", newqid)

        if "results" in h5in:
            with h5py.File(fnout, "a") as h5out:
                for run in h5in["results"]:
                    newqid = qids[get_qid(run)]
                    new = "results/run_"+newqid
                    old = "results/" + run

                    for field in ["CC", "CFLAGS", "repository"]:
                        h5out[new].attrs[field] = h5in[old].attrs[field]

                    for field in ["bfield", "efield", "options", "plasma",
                                  "wall", "neutral", "marker"]:
                        h5out[new].attrs["qid_"+field] = \
                            np.string_( qids[h5in[old].attrs["qid_"+field].decode('utf-8')] )
