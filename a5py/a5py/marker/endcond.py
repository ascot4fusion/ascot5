"""
Contains class representation of Ascot5 endconditions.

File: endcond.py
"""

import numpy as np

from .alias import get as alias

endconds = {
    "none"    : 0,
    "aborted" : -1,
    "emin"    : 1,
    "therm"   : 2,
    "wall"    : 3,
    "rhomin"  : 4,
    "rhomax"  : 5,
    "polmax"  : 6,
    "tormax"  : 7,
    "cpumax"  : 8,
}

class Endcond():
    """
    Class for representing endconditions.

    Each Endcond represents a single marker's endcondition value in Ascot5
    ini/endstate. Since a marker may have multiple endconditions active,
    this class provides a functionality for comparing endconditions with a
    simple equal operator e.g. Endcond == "wall" is a valid comparison.
    """

    def __init__(self, rawendcond, errmsg):
        """
        Init endcondition based on endcond and errmod fields in output data.
        """
        self._endconds = []
        if errmsg > 0:
            self._endconds.append(endconds["aborted"])

        if rawendcond == 0:
            self._endconds.append(endconds["none"])

        ec = rawendcond - np.power(2, endconds["cpumax"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["cpumax"])

        ec = rawendcond - np.power(2, endconds["tormax"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["tormax"])

        ec = rawendcond - np.power(2, endconds["polmax"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["polmax"])

        ec = rawendcond - np.power(2, endconds["rhomax"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["rhomax"])

        ec = rawendcond - np.power(2, endconds["rhomin"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["rhomin"])

        ec = rawendcond - np.power(2, endconds["wall"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["wall"])

        ec = rawendcond - np.power(2, endconds["therm"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["therm"])

        ec = rawendcond - np.power(2, endconds["emin"])
        if ec >= 0:
            rawendcond = ec
            self._endconds.append(endconds["emin"])


    def __eq__(self, other):
        """
        Test if this endcondition is equal to the given value.
        """
        if isinstance(other, str):
            other = endconds[alias(other)]
            for key in self._endconds:
                if key == other:
                    return True

        return False

    def __neq__(self, other):
        """
        Implementation for != operator.
        """
        return not (self == other)
