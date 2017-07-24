import numpy as np

def read(orbits, mode):
    out = {}
    for orbgroup in orbits:
        out[orbgroup] = {}
        for field in orbits[orbgroup]:
            if mode != 0:
                out[orbgroup][field]           = orbits[orbgroup][field][:]
                out[orbgroup][field + "_unit"] = orbits[orbgroup][field].attrs["unit"]
            else:
                out[orbgroup][field]           = 0
                out[orbgroup][field + "_unit"] = orbits[orbgroup][field].attrs["unit"]

    # Read-raw ends here
    if mode == 3:
        return out

    for orbgroup in orbits:
        if "id" in orbits[orbgroup]:
            out[orbgroup]["N"]        = (np.unique(orbits[orbgroup]["id"][:])).size
            out[orbgroup]["uniqueId"] = np.unique(orbits[orbgroup]["id"][:])
        else:
            out[orbgroup]["N"] = 0
    
    # Show ends here
    if mode == 0:
        return out 

    # Sort fields by id (major) and time (minor), both ascending
    for orbgroup in orbits:
        if out[orbgroup]["N"] > 0:
            ind = np.lexsort((out[orbgroup]["id"], out[orbgroup]["time"]))
            for field in orbits[orbgroup]:
                out[orbgroup][field] = out[orbgroup][field][ind]
    
    # Read ends here
    if mode == 1:
        return out

    # Non-SI fields to SI units
    return out
