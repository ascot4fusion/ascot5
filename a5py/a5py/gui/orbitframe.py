"""
Contains definition of StateFrame class.

File: orbitframe.py
"""
import tkinter

from .plotframe import PlotFrame


class OrbitFrame(PlotFrame):
    """
    An opening frame where other frames can be accessed and HDF5 file modified.

    File: indexframe.py
    """

    def __init__(self, gui, orbits):
        """
        Initialize
        """
        super().__init__(gui)
        self._orbits = orbits

        self._plot("R", "z")

    def _plot(self, *args):
        fig = self.get_fig()
        ax = fig.add_subplot(1,1,1)
        self._orbits.plot(args[0], args[1], equal=True, axes=ax)
