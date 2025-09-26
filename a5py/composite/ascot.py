"""Contains the definition of the Ascot class which is the main interface.

Ascot class acts as an interface to other tools e.g. AscotData (interface to the
data), BeamRun (interface to the beam simulation), etc. It also acts as an
interface to input data interpolation. While each input data itself has the
tools to process it's own data, the logical place for processing data that
requires multiple inputs is here.
"""
import os

from a5py.data import AscotData
from a5py.engine.simulate import Simulate

class Ascot(Simulate):
    """Primary tool for executing and processing ASCOT5 simulations and data.

    Attributes
    ----------
    data : :class:`.AscotData`
        Contains the simulation inputs and outputs.
    """

    def __init__(self, inputfile=None, create=False):
        """Initialize Ascot instance.

        Parameters
        ----------
        inputfile : str, optional
            Name of the HDF5 file or `None` to create an empty instance.
        create : bool, optional
            Create a new HDF5 file with given name.
        """
        super().__init__()

        self.data: AscotData = None
        #self.biosaw: BioSaw = BioSaw(self)
        #self.afsi       = Afsi(self)
        #self.markergen  = MarkerGenerator(self)
        if inputfile is None:
            self.data = AscotData()
        else:
            self.data = AscotData((inputfile, not create))
