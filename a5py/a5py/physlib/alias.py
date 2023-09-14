"""
Module containing widely used synonyms for different quantities.

Examples:
    Check if quantity is mass:
    isalias("m", "mass")
    True

    Get standard name of the quantity:
    getalias("q")
    "charge"

Note that the synonyms are case insensitive and ignore spaces, lowercases and
dashes.

File: alias.py
"""

## Standard name : list of synonyms dictionary. Flag indicates whether quantity
## is common for both gyro-orbit and guiding-center.
synonyms = {
    "ids"      : ["id"],
    "mass"     : ["m"],
    "charge"   : ["q"],
    "vnorm"    : ["v", "velocity", "velocitynorm"],
    "bnorm"    : ["b", "magneticfield", "magneticfieldnorm"],
    "energy"   : ["ekin", "e"],
    "pitch"    : ["xi"],
    "mu"       : ["magneticmoment"],
    "gamma"    : ["lorentzfactor", "relativisticfactor"],
    "phi"      : ["tor", "toroidalangle"],
    "phimod"   : ["tormod"],
    "theta"    : ["pol", "poloidalangle"],
    "thetamod" : ["polmod"],
    "zeta"     : ["gyroangle"],
    "vpar"     : ["parallelvelocity", "vparallel", "vpa", "vpara"],
    "ppar"     : ["parallelmomentum", "pparallel", "ppa", "ppara"],
    "vperp"    : ["perpendicularvelocity", "vperpendicular", "vpe"],
    "pparp"    : ["perpendicularmomentum", "pperpendicular", "ppe"],
    "pnorm"    : ["p", "momentum", "momentumnorm"],
    "ptor"     : ["toroidalcanonicalangularmomentum", "ctor", "torcanangmom"]
}


def isalias(name, standardname):
    """
    Check if name is an alias of the standardname.
    """

    name = name.lower()
    name = name.replace(" ", "")
    name = name.replace("-", "")
    name = name.replace("_", "")

    suffix = ""
    if len(name) > 3 and name[-3:] == "prt":
        suffix = "prt"
        name   = name[:-3]

    return name in synonyms[standardname]


def getalias(name):
    """
    Get standard name of the given quantity.
    """

    name = name.lower()
    name = name.replace(" ", "")
    name = name.replace("-", "")
    name = name.replace("_", "")

    suffix = ""
    if len(name) > 3 and name[-3:] == "prt":
        suffix = "prt"
        name   = name[:-3]

    for key in synonyms.keys():
        if name in synonyms[key]:
            return key + suffix

    return name + suffix
