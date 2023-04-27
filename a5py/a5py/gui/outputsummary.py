import tkinter as tk
from tkinter import ttk
import numpy as np

from matplotlib.gridspec import GridSpec
from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class SummaryFrame(ttk.Frame):

    def init(self, gui, canvas):
        self.gui = gui
        self.canvas = canvas

        return self


    def display(self):
        run = self.gui.ascot.hdf5.active
        self.canvas.fig.clear()

        run.inistate.histogram(
            x="rho", y=None, xbins=(0,1,100),
            ybins=None, weight=True, logx=False,
            logy=False, logscale=False,
            endcond=None, axes=self.canvas.fig.axis[0],
            iniend=[0,0])

        run.endstate.histogram(
            x="mileage", y=None, xbins=(1e-6, 1.0, 100),
            ybins=None, weight=True, logx=True,
            logy=False, logscale=False,
            endcond=None, axes=self.canvas.fig.axis[1],
            iniend=[1,1])

        self.canvas.fig.draw()


class SummaryCanvas(ttk.Frame):

    def init(self):

        class Plot3DFrame(PlotFrame):

            def set_axes(self):
                gs = GridSpec(2,2)
                ax = []
                ax += [self.fig.add_subplot(gs[0,1])]
                ax += [self.fig.add_subplot(gs[1,1])]
                return ax

        fig = Plot3DFrame(self)
        fig.place(relheight=0.8, anchor="nw")

        self.fig = fig
        return self
