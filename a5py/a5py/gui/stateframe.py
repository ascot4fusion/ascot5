"""
Contains definition of StateFrame class.

File: stateframe.py
"""
import tkinter

from .plotframe import PlotFrame


class StateFrame(PlotFrame):
    """
    An opening frame where other frames can be accessed and HDF5 file modified.

    File: indexframe.py
    """

    def __init__(self, gui, inistate, endstate=None):
        """
        Initialize
        """
        super().__init__(gui)
        self._inistate = inistate
        self._endstate = endstate

        self._plot_inistate()

    def _plot_inistate(self):
        pass
