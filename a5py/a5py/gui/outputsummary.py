"""
Summary (frame) for the output content.

File: outputsummary.py
"""
import tkinter as tk
from tkinter import ttk
import numpy as np

from matplotlib.gridspec import GridSpec
from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class SummaryFrame(ttk.Frame):
    """
    Settings frame summarizing simulation output.

    The summary canvas is not interactive, so this frame does not contain any
    settings. Instead, it has a text box that shows sum of markers by
    end condition and possible error messages.
    """

    def init(self, gui, canvas):
        """
        Initializes frame by creating all the widgets.
        """
        self.markersummary = tk.Text(self, width=50)
        self.markersummary.pack(fill="both", expand=True)
        self.markersummary.config(state="disabled")

        self.gui = gui
        self.canvas = canvas

        return self


    def viewsummarytext(self, run):
        self.markersummary.config(state="normal")
        self.markersummary.delete("1.0", "end")

        text = ""
        msg = run.getstate_markersummary()
        for m in msg:
            text += m + "\n"

        self.markersummary.insert("end", text)
        self.markersummary.config(state="disabled")


    def display(self):
        """
        Display this frame and redraws canvas.
        """
        run = self.gui.ascot.data.active
        self.viewsummarytext(run)
        self.canvas.clear()

        run.plotstate_summary(axes_inirho=self.canvas.axes_inirho,
                              axes_endrho=self.canvas.axes_endrho,
                              axes_mileage=self.canvas.axes_mileage,
                              axes_energy=self.canvas.axes_energy,
                              axes_rz=self.canvas.axes_rz,
                              axes_rhophi=self.canvas.axes_rhophi)

        self.canvas.draw()


class SummaryCanvas(PlotFrame):
    """
    Canvas for showing summary plots.

    This canvas creates several axes where the plots can be placed.
    """

    def init(self):
        """
        Name axes just so that they are easier to reference when plotting.
        """
        self.axes_inirho  = None
        self.axes_endrho  = None
        self.axes_mileage = None
        self.axes_energy  = None
        self.axes_rz      = None
        self.axes_rhophi  = None
        return self


    def set_axes(self):
        """
        Override the PlotFrame method in order to create multiple axes.
        """
        gs = GridSpec(2,4)
        self.axes_inirho  = self.fig.add_subplot(gs[0,0])
        self.axes_endrho  = self.fig.add_subplot(gs[0,1])
        self.axes_mileage = self.fig.add_subplot(gs[1,0])
        self.axes_energy  = self.fig.add_subplot(gs[1,1])
        self.axes_rz      = self.fig.add_subplot(gs[:,2])
        self.axes_rhophi  = self.fig.add_subplot(gs[:,3])

        ax = [self.axes_inirho, self.axes_endrho, self.axes_mileage,
              self.axes_energy, self.axes_rz, self.axes_rhophi]
        return ax
