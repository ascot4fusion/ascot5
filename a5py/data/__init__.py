"""Contains all data classes for simulation inputs and outputs, and the main
data structure class.
"""
from .access import InputVariant

from .base import AscotData
from .bfield import Bfield
#from .efield import Efield
#from .plasma import Plasma
#from .neutral import Neutral
#from .wall import Wall
#from .mhd import Mhd
#from .boozer import BoozerMap



class MetaData():
    pass

class InputGroup(InputVariant):
    """Node containing input data groups.
    """


#class RunGroup(RunMixin):
    """Node containing results and methods to process them.
    """


#class AfsiGroup(AfsiMixin):
    """Node containing AFSI results and methods to process them.
    """


#class BBNBIGroup(BBNBIMixin):
    """Node containing BBNBI results and methods to process them.
    """
