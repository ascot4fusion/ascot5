"""
Contains definition of BfieldFrame class.

File: bfieldframe.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np

from a5py.ascotpy.libbfield import LibBfield

from .fieldframe import FieldFrame

class BfieldFrame(FieldFrame):
    """
    A frame for plotting magnetic field data.
    """

    def __init__(self, gui, ascotpy, walldata=None):
        """
        Initialize and show default plot.
        """
        ascotpy.init(bfield=True)
        super().__init__(gui, ascotpy, LibBfield.quantities)

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

        phi  = 0
        time = 0

        axes = fig.add_subplot(1,1,1)
        self.ascotpy.plotRz(r, phi, z, time, self._qchoice.get(), axes)
        self.ascotpy.plotseparatrix(r, phi, z, time, axes)

        if self._walldata is not None:
            self._walldata.plotRz(axes=axes)

        self.draw()
