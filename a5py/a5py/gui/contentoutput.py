import tkinter as tk
from tkinter import ttk

from .outputstate import StateFrame

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

        framestate = StateFrame(settings).init()
        frameorbit = ttk.Frame(settings)
        framelost  = ttk.Frame(settings)
        framedist  = ttk.Frame(settings)

        settings.add(framestate, text="Ini/Endstate")
        settings.add(frameorbit, text="Orbit")
        settings.add(framelost,  text="Losses")
        settings.add(framedist,  text="Dists")

        settings.pack(fill="both", expand=True)

        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            if tab == "Ini/Endstate":
                pass

        settings.bind('<<NotebookTabChanged>>', on_tab_change)

        self.canvas = canvas

    def display(self):
        pass

        self.gui.ascot.init_from_run(
            self.gui.ascot.hdf5.active,
            bfield=True, efield=True, neutral=True,
            plasma=True, boozer=True, mhd=True)
