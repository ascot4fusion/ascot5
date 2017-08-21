"""
Tools to modify HDF5 files.
"""
import numpy as np
import h5py
from __future__ import a5py.ascot5io.ascot5 as ascot5

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
    ot.attrs["type"] = np.string_(subfield)
    fs.copy(field + "/" + subfield, ot, name=subfield)

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

    print("Combining output to " + fnt + "with mode " + mode)

    # Open target file.
    ft = h5py.File(fnt, "a")

    for state in ["inistate", "endstate"]:
        print("Combining " + state)

        target = ascot5.read(fnt,state)
        target = target["states"][state]

        for fn in fns:
            source = ascot5.read(fn,state)
            
            for field in target:
                if target[field][-4:] != "unit": # No need to combine unit specifiers
                    target[field].extend(source[field])
            

        ft["states"][state] = target
        #states.write(ft,fnt,state)


    print("Combining distributions.")

    target = ascot5.read(fnt,"dists")
    for fn in fns:
        source = ascot5.read(fn,"dists")
        target["dists"]["ordinate"] += source["dists"]["ordinate"]
        

    ft["distributions"]["ordinate"] = target["dists"]["ordinate"]
    #dists.write(ft,fnt)

    print("Combining orbits.")
    
    target = ascot5.read(fnt,"orbits")
    for fn in fns:
        source = ascot5.read(fn,"orbits")

        for orbgroup in ft["orbits"]:
            target = target["orbits"][state]
            
            
            for field in target[orbgroup]:
                if target[orbgroup][field][-4:] != "unit": # No need to combine unit specifiers
                    target[orbgroup][field].extend(source[orbgroup][field])
            

        ft["orbits"][orbgroup] = target["orbits"][orbgroup]
        #orbits.write(ft,fnt)

    # Clean.
    ft.close()
