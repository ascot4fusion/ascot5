"""
Contains class representation of Ascot5 endconditions.

File: endcond.py
"""

import numpy as np

from .alias import get as alias

endconds = [
    ("aborted",  -1, 0x2),
    ("none",      0, 0x1),
    ("tmax",      1, 0x4),
    ("emin",      2, 0x8),
    ("therm",     3, 0x10),
    ("wall",      4, 0x20),
    ("rhomin",    5, 0x40),
    ("rhomax",    6, 0x80),
    ("polmax",    7, 0x100),
    ("tormax",    8, 0x200),
    ("cpumax",    9, 0x400)
]


def getval(name):
    for ec in endconds:
        if ec[0] == name:
            return ec[1]

    return None


def getname(binr):
    estr = ""
    for ec in endconds:
        if ec[2] & binr:
            if estr == "":
                estr = ec[0]
            else:
                estr += " + " + ec[0]

    return estr


def getbin(name):
    for ec in endconds:
        if ec[0] == name:
            return ec[2]

    return None


def parse(endcond):
    names = []
    vals  = []

    if endcond == 0:
        names.append("none")
        vals.append(getval("none"))
        return (names, vals)

    for ec in endconds:
        if endcond & ec[2]:
            names.append(ec[0])
            vals.append(ec[1])

    return (names, vals)
