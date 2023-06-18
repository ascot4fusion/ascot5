import tkinter as tk
from tkinter import ttk

from .outputsummary import SummaryFrame, SummaryCanvas
from .outputstate   import StateFrame, StateCanvas
from .outputorbit   import OrbitFrame, OrbitCanvas
from .outputloss    import LossFrame, LossCanvas
from .outputdist    import DistFrame, DistCanvas

class ContentOutput:
    """
    Class for viewing and analyzing results.

    Settings panel consists of tabs that open different plot views on the
    canvas.

    Overview:
    - Provides overview of the simulation results.
    - Summary text of endstates and aborted markers.
    - Ini-rho histogram colored with endstate.
    - Final energy histogram (log) colored with endstate.
    - Mileage histogram (log) colored with endstate
    - Rz-plot of marker endstates colored with endstate.

    State:
    - Plot any marker quantity with respect to another marker quantity
    - Scatter (2D/3D + color) and histogram (1D and 2D) plots.
    - Can choose between ini/end coordinates.
    """

    def __init__(self, gui, settings, canvas):
        self.gui = gui

        settings = ttk.Notebook(settings)

        canvassummary = SummaryCanvas(canvas).init()
        canvasstate   = StateCanvas(canvas).init()
        canvasorbit   = OrbitCanvas(canvas).init()
        canvasloss    = LossCanvas(canvas).init()
        canvasdist    = DistCanvas(canvas).init()

        self.framesummary = SummaryFrame(settings).init(gui, canvassummary)
        self.framestate   = StateFrame(settings).init(gui, canvasstate)
        self.frameorbit   = OrbitFrame(settings).init(gui, canvasorbit)
        self.framedist    = DistFrame(settings).init(gui, canvasdist)
        self.frameloss    = LossFrame(settings).init(gui, canvasloss)

        settings.add(self.framesummary, text="Summary")
        settings.add(self.framestate,   text="Ini/Endstate")
        settings.add(self.frameorbit,   text="Orbit")
        settings.add(self.framedist,    text="Dists")
        settings.add(self.frameloss,    text="Losses")

        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            if tab == "Summary":
                self.active_canvas.pack_forget()
                self.framesummary.display()
                self.active_canvas = self.framesummary.canvas
                self.active_canvas.pack(fill="both", expand=True)
            elif tab == "Ini/Endstate":
                self.active_canvas.pack_forget()
                self.framestate.display()
                self.active_canvas = self.framestate.canvas
                self.active_canvas.pack(fill="both", expand=True)
            elif tab == "Orbit":
                self.active_canvas.pack_forget()
                self.frameorbit.display()
                self.active_canvas = self.frameorbit.canvas
                self.active_canvas.pack(fill="both", expand=True)
            elif tab == "Dists":
                self.active_canvas.pack_forget()
                self.framedist.display()
                self.active_canvas = self.framedist.canvas
                self.active_canvas.pack(fill="both", expand=True)
            elif tab == "Losses":
                self.active_canvas.pack_forget()
                self.frameloss.display()
                self.active_canvas = self.frameloss.canvas
                self.active_canvas.pack(fill="both", expand=True)

        settings.bind('<<NotebookTabChanged>>', on_tab_change)
        settings.pack(fill="both", expand=True)

        self.active_canvas = ttk.Frame(canvas)
        self.active_canvas.pack()

        self.canvas = canvas

    def display(self):

        self.gui.ascot.init_from_run(
            self.gui.ascot.data.active,
            bfield=True, efield=True, neutral=True,
            plasma=True, boozer=True, mhd=True)
