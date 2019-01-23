"""
Module containing widely used synonyms for different quantities.

Examples:
    Check if quantity is mass:
    check("m", "mass")
    True

    Get standard name of the quantity:
    get("q")
    "charge"

Note that the synonyms are case insensitive and ignore spaces, lowercases and
dashes.

File: alias.py
"""

## Standard name : list of synonyms dictionary.
synonyms = {
    "mass"   : ["m"],
    "charge" : ["q"],
    "vnorm"  : ["v", "velocity", "velocitynorm"],
    "bnorm"  : ["b", "magneticfield", "magneticfieldnorm"],
    "energy" : ["ekin", "e"],
    "pitch"  : ["xi"],
    "mu"     : ["magneticmoment"],
    "gamma"  : ["lorentzfactor", "relativisticfactor"],
    "phi"    : ["toroidalangle"],
    "pol"    : ["poloidalangle"],
    "theta"  : ["gyroangle"],
    "vpar"   : ["parallelvelocity", "vparallel", "vpa", "vpara"],
    "ppar"   : ["parallelmomentum", "pparallel", "ppa", "ppara"],
    "vperp"  : ["perpendicularvelocity", "vperpendicular", "vpe"],
    "pparp"  : ["perpendicularmomentum", "pperpendicular", "ppe"],
    "pnorm"  : ["p", "momentum", "momentumnorm"],
    "ptor"   : ["canonicaltoroidalmomentum"],
    "wall"   : ["lost"],
    "therm"  : ["thermalized"]
}

def check(name, standardname):
    """
    Check if name is alias of the standardname
    """

    name = name.lower()
    name = name.replace(" ", "")
    name = name.replace("-", "")
    name = name.replace("_", "")

    return name in synonyms[standardname]

def get(name):
    """
    Get standard name of the given quantity.
    """

    name = name.lower()
    name = name.replace(" ", "")
    name = name.replace("-", "")
    name = name.replace("_", "")

    for key in synonyms.keys():
        if name in synonyms[key]:
            return key

    return name
