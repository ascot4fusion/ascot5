"""
Tools to modify HDF5 files.
"""
import numpy as np
import h5py
from . import ascot5
from . import states
from . import orbits
from . import dists

def remove(fn,runid=0):
    """
    Remove all or specific runs.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
         Full path to HDF5 file.
    runid : int, optional
         Id for for the run to be removed. By default all runs
         are removed.
    """

    f = h5py.File(fn, "a")

    
    if  "inistate" in f:
        del f["inistate"]
    if  "endstate" in f:
        del f["endstate"]
    if  "orbits" in f:
        del f["orbits"]
    if  "distributions" in f:
        del f["distributions"]
    
    f.close()

def rmgroup(fn,group):
    """
    Remove a group.

    TODO Not compatible with new HDF5 format.

    Parameters
    ----------

    fn : str
         Full path to HDF5 file.
    group : str
         Group to be removed
    """
    f = h5py.File(fn, "a")
    if not group in f:
        print("Error: Source file does not contain the field to be remove")
        return

    del f[group]
    f.close()


def copy(fns,fnt,field,subfield):
    """
    Copy input field.

    The copied field is set as the active field in target HDF5 file.

    TODO not compatible with new HDF5 format.

    Parameters
    ----------
    fns : str
          Full path to source file.
    fnt : str
          Full path to target file.
    field : str 
          Master field where copying is done e.g. bfield.
    subfield : str
          Subfield to be copied e.g. B2D.
    """

    # Get the target field and type from source
    fs = h5py.File(fns, "r")
    o = fs[field]
    types = o.attrs["type"]
    
    if not field in fs:
        print("Error: Source file does not contain the field to be copied")
        return

    # Check if the field in target exists (otherwise it is created)
    # Delete possible existing types and sub fields of same type as copied
    ft = h5py.File(fnt, "a")
    if not field in ft:
        ot = ft.create_group(field)
    else:
        ot = ft[field]
        del ot.attrs["type"]
        if subfield in ot:
            del ot[subfield]

    # Do the copying and set the type
    fs.copy(field + "/" + subfield, ot, name=subfield)
    ft[field].attrs["type"] = np.string_(subfield)
    ft[field][subfield].attrs["qid"] = fs[field][subfield].attrs["qid"]
    ft[field][subfield].attrs["date"] = fs[field][subfield].attrs["date"]
    

    fs.close()
    ft.close()


def combine(fnt, fns, mode="add"):
    """
    Combine output of multiple HDF5 files into one.

    TODO Not compatible with new HDF5 format.

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
    

    for state in ["inistate", "endstate"]:
        print("Combining " + state)

        target = ascot5.read_hdf5(fnt,"states")

        # Check whether target has the desired state.
        # Init empty state if necessary.
        if state not in target["states"]:
            target["states"][state] = {}

        target = target["states"]

        # Iterate over source files
        for fn in fns:
            source = ascot5.read_hdf5(fn,"states")

            # Check that source contains the state
            if state in source["states"]:
                source = source["states"][state]
                
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
        states.write_hdf5(fnt,target,0)

    
    print("Combining distributions.")

    target = ascot5.read_hdf5(fnt,"dists")

    # If target does not have distributions, we get them from the first
    # source file.
    source = ascot5.read_hdf5(fns[0],"dists")
    for dist in source["dists"]:
        if not dist in target["dists"]:
            target["dists"][dist] = source["dists"][dist]
            target["dists"][dist]["ordinate"] = target["dists"][dist]["ordinate"] * 0


    # Sum ordinates in all distributions (same for both modes)
    for fn in fns:
        source = ascot5.read_hdf5(fn,"dists")
        for dist in source["dists"]:
            target["dists"][dist]["ordinate"] += source["dists"][dist]["ordinate"]

    dists.write_hdf5(fnt,target["dists"],0)

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
    orbits.write_hdf5(fnt,target,0)
