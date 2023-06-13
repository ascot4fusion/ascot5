"""Tools for assessing physical quantities and units.
"""
from .coords import *
from .gamma import *
from .particle import *

from .alias import isalias, getalias

# Make some frequently needed physical constants easily accessible
import unyt
import scipy.constants as constants

e       = constants.elementary_charge * unyt.C
m_e     = constants.physical_constants["electron mass"][0] * unyt.kg
m_p     = constants.physical_constants["proton mass"][0] * unyt.kg
m_a     = constants.physical_constants["alpha particle mass"][0] * unyt.kg
c       = unyt.c
eps_0   = unyt.eps_0
