"""
Contains definition of WallFrame class.

File: wallframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from .plotframe import PlotFrame

class WallFrame(PlotFrame):
    """
    A frame for plotting wall inputs.
    """

    def __init__(self, gui, wall):
        """
        Initialize and show default plot.
        """

        super().__init__(gui)
        self._gui  = gui
        self._wall = wall

        self._plot()


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig  = self.get_fig()
        axes = fig.add_subplot(1,1,1)
        self._wall.plotRz(axes=axes)

        self.draw()
