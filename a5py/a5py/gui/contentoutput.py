import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

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

        settingsframe = ttk.Notebook(settings)

        class StateFrame(tk.Frame):

            def init(self, canvas):
                frame = tk.Frame(self)

                self.sctr_xcrd = DropdownMenu(frame, log=True, label="x: ")
                self.sctr_ycrd = DropdownMenu(frame, log=True, label="y: ")
                self.sctr_zcrd = DropdownMenu(frame, log=True, label="z: ")
                self.sctr_ccrd = DropdownMenu(frame, log=True, label="c: ")
                self.sctr_endc = DropdownMenu(frame, label="Endcond: ")
                self.sctr_axeq = Tickbox(frame, label="Axis equal")

                self.sctr_xcrd.pack()
                self.sctr_ycrd.pack()
                self.sctr_zcrd.pack()
                self.sctr_ccrd.pack()
                self.sctr_endc.pack()
                self.sctr_axeq.pack()

                self.scatterframe = frame

                self.visible = None

                self.plotframe = PlotFrame(canvas)
                

            def show(self, frame):
                if frame != self.visible:
                    self.visible = frame
                    if frame == "scatter":
                        self.histframe.pack_forget()
                        self.scatterframe.pack()
                    if frame == "hist":
                        self.scatterframe.pack_forget()
                        self.histframe.pack()

            def plot(self):
                self.sctr_xcrd.get_val()

        stateframe = StateFrame(settingsframe)
        plotcanvas = tk.Frame(canvas)
        fig_dist2d = PlotFrame(plotcanvas)
        fig_dist2d.place(relheight=0.8, anchor="nw")

        settingsframe.add(stateframe, text="Ini/Endstate")


        def on_tab_change(event):
            tab = event.widget.tab('current')['text']
            if tab == "Ini/Endstate":
                pass

        settingsframe.bind('<<NotebookTabChanged>>', on_tab_change)

        self.settingsframe = settingsframe
        self.plotcanvas    = plotcanvas

    def display(self):


        self.gui.ascot.init_from_run(
            self.gui.ascot.hdf5.active,
            bfield=True, efield=True, neutral=True,
            plasma=True, boozer=True, mhd=True)
