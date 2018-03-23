"""
Tools to modify HDF5 files.
"""
import numpy as np
import h5py
from . import ascot5
from . import ascot5group
from . import states
from . import orbits
from . import dists

def setactivegroup(fn, mastergroup, qid):
    """
    Set a group active.

    Parameters
    ----------

    fn : str
         Full path to HDF5 file.
    mastergroup : str
         Mastergroup where the group belongs to.
    qid: str
         Id for the active group.
    """
    f = h5py.File(fn, "a")
    exists = False
    if mastergroup in f:
        for r in f[mastergroup]:
            if r[-10:] == qid:
                ascot5group.setactive(f, mastergroup + "/" + r)
                exists = True
                
    if exists == None:
        print("Error: group not found.")
    
    f.close()
    

def removerun(fn, qid=None):
    """
    Remove all or specific runs.

    Parameters
    ----------

    fn : str
         Full path to HDF5 file.
    qid : str, optional
         Id for for the run to be removed. By default all runs
         are removed.
    """

    f = h5py.File(fn, "a")

    for r in f["results"]:
        if qid == None:
            del f["results"][r]
        elif r[-10:] == qid:
            del f["results"][r]
    
    f.close()

def removegroup(fn, mastergroup, qid):
    """
    Remove a group.

    Parameters
    ----------

    fn : str
         Full path to HDF5 file.
    mastergroup : str
         Mastergroup where the group belongs to.
    qid: str
         Id for the group to be removed
    """
    f = h5py.File(fn, "a")
    if not mastergroup in f:
        print("Error: Source file does not contain the group to be removed.")
        return

    for r in f[mastergroup]:
        if r[-10:] == qid:
            del f[mastergroup][r]
            f.close()
            return
            
    print("Error: Source file does not contain the group to be removed.")


def copygroup(fns, fnt, mastergroup, qid):
    """
    Copy input field.

    The copied field is set as the active field in target HDF5 file.

    Parameters
    ----------
    fns : str
          Full path to source file.
    fnt : str
          Full path to target file.
    mastergroup : str 
          Mastergroup where the group belongs to.
    qid : str
          Id for the group to be removed.
    """

    # Get the target field and type from source
    fs = h5py.File(fns, "r")
    
    if not mastergroup in fs:
        print("Error: Source file does not contain the group to be copied")
        return

    field = None
    for r in fs[mastergroup]:
        if r[-10:] == qid:
            field = r
            break
    if field == None:
        print("Error: Source file does not contain the group to be copied")
        return

    # Check if the mastergroup in target exists (otherwise it is created)
    ft = h5py.File(fnt, "a")
    if not mastergroup in ft:
        ft.create_group(mastergroup)
    
    # Copy group into target
    ft[mastergroup].attrs["active"] = qid
    fs.copy(mastergroup + "/" + field, ft[mastergroup], name=field)    

    fs.close()
    ft.close()


def combineresults(fnt, fns, mode="add"):
    """
    Combine output of multiple HDF5 files into one.

    Since files can have multiple run, this function affects lthose
    that are defined as being "active" or most recent in the results group.

    Parameters
    ----------

    fnt : str
        Full path to HDF5 file where combined output is added.
    fns: str list
        List of HDF5 filenames from which the output is read.
    mode, str=["add", "continue"]
        Specifies how output is combined.
        "add" adds the output assuming "fnt" and "fns" are independent 
        simulations.
        "continue" assumes "fns" is continued simulation of "fnt".
    """

    print("Combining output to " + fnt + " with mode " + "\"" + mode + "\"")

    f = h5py.File(fnt, "a")
    # If target does not have results, we get them from the first
    # source file.
    if "results" not in f:
        qid = ascot5.get_qids(fns[0], "results")[0][0]
        path = "results/run-" + qid
        print("Creating 'results' field for target.")
        fs = h5py.File(fns[0])
        fs.copy("results", f)
        for state in ["inistate","endstate"]: # f[path]:
            for field in f[path][state]:
                del f[path][state][field]
                f[path][state].create_dataset(field, (0,))
                for a in fs[path][state][field].attrs.items():
                    f[path][state][field].attrs.create(a[0],a[1])
        del f[path]["dists"]
        fs.close()
    f.close()
        
    qid = ascot5.get_qids(fnt, "results")[0][0]
    path = "run-" + qid

    for state in ["inistate", "endstate"]:
        print("Combining " + state)

        target = ascot5.read_hdf5(fnt,"results")
        # Check whether target has the desired state.
        # Init empty state if necessary.
        if state not in target["results"][path]:
            target["results"][path][state] = {}

        target = target["results"][path]

        # Iterate over source files
        for fnind, fn in enumerate(fns):

            source = ascot5.read_hdf5(fn,"results")

            # Check that source contains the state
            if state in source["results"][path]:                
                source = source["results"][path][state]
                # Iterate over all fields in a state
                for field in source:
                    # If target does not have this field, add it.
                    # Otherwise extend or replace it depending on mode.
                    if field not in target[state]:
                        target[state][field] = source[field]
                    # No need to combine unit specifiers (fields ending with "unit")
                    elif field[-4:] != "unit" and field != "N" and field != "uniqueId":
                        if mode == "add":
                            target[state][field] = np.concatenate((target[state][field],source[field]))
                        elif mode == "continue" and state != "inistate":
                            idx = np.argsort(target[state]["id"][:])
                            target[state]["id"][:] = target[state]["id"][idx]
                            idx = np.argsort(source["id"][:])
                            fld = source[field][idx]

                            idx = np.isin(target[state]["id"][:], source["id"][:], assume_unique=True)
                            target[state][field][idx] = fld[:]
                            print(np.sum(idx==True))

        # Target now contains all data from combined runs, so we just need to write it.
        states.write_hdf5(fnt,target,qid,state)
    
    print("Combining distributions.")

    target = ascot5.read_hdf5(fnt,"results")
    
    # If target does not have distributions, we get them from the first
    # source file.
    source = ascot5.read_hdf5(fns[0],"results")
    for dist in source["results"][path]["dists"]:
        if not "dists" in target["results"][path]:
            target["results"][path]["dists"] = {}
        if not dist in target["results"][path]["dists"]:
            target["results"][path]["dists"][dist] = source["results"][path]["dists"][dist]
            target["results"][path]["dists"][dist]["ordinate"] = target["results"][path]["dists"][dist]["ordinate"] * 0

    # Sum ordinates in all distributions (same for both modes)
    for fn in fns:
        source = ascot5.read_hdf5(fn,"results")
        for dist in source["results"][path]["dists"]:
            target["results"][path]["dists"][dist]["ordinate"] += source["results"][path]["dists"][dist]["ordinate"]

    dists.write_hdf5(fnt,target["results"][path]["dists"],qid)
    
    f = h5py.File(fnt, "a")
    if "orbits" in f["results/" + path]:
        print("Combining orbits.")        
        target = ascot5.read_hdf5(fnt,"orbits")["orbits"]
        
        # Iterate over source files
        for fn in fns:
            source = ascot5.read_hdf5(fn,"orbits")["orbits"]
            
            # Check whether target has the desired state.
            # Init empty state if necessary.
            for orbgroup in source:
                if orbgroup not in target:
                    target[orbgroup] = {}

                # Iterate over all fields in a orbit
                for field in source[orbgroup]:
                    
                    # If target does not have this field, add it.
                    # Otherwise extend it.
                    if field not in target[orbgroup]:
                        target[orbgroup][field] = source[orbgroup][field]
                        
                    # No need to combine unit specifiers (fields ending with "unit")
                    elif field[-4:] != "unit" and field != "N" and field != "uniqueId":
                        target[orbgroup][field] = np.concatenate((target[orbgroup][field],source[orbgroup][field]))
                        
        # Target now contains all data from combined runs, so we just need to write it.
        orbits.write_hdf5(fnt,target,qid)
