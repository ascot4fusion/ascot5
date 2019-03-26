"""
Contains definition of BoozerFrame class.

File: boozerframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from a5py.ascotpy.libboozer import LibBoozer

from .fieldframe import FieldFrame

class BoozerFrame(FieldFrame):
    """
    A frame for plotting boozer data.
    """

    def __init__(self, gui, ascotpy, walldata=None):
        """
        Initialize and show default plot.
        """
        ascotpy.init(bfield=True, boozer=True)
        super().__init__(gui, ascotpy, LibBoozer.quantities)

        self._walldata = walldata

        self._plot()


    def _plot(self, *args):
        """
        Read control states and plot.
        """
        fig = self.get_fig()

        r = np.linspace( float(self._binxminchoice.get()),
                         float(self._binxmaxchoice.get()),
                         float(self._nbinxchoice.get()) )

        z = np.linspace( float(self._binyminchoice.get()),
                         float(self._binymaxchoice.get()),
                         float(self._nbinychoice.get()) )

        phi  = float(self._phichoice.get()) * np.pi / 180
        time = 0

        axes = fig.add_subplot(1,1,1)
        self.ascotpy.plotRz(r, phi, z, time, self._qchoice.get(), axes)
        self.ascotpy.plotseparatrix(r, phi, z, time, axes)

        if self._walldata is not None:
            self._walldata.plotRz(axes=axes)

        self.draw()
