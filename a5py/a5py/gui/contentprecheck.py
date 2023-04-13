import tkinter as tk
import matplotlib.pyplot as plt
from a5py.ascotpy.preflight import *
from tkinter import ttk
from .components import PlotFrame

class ContentPrecheck:

    def __init__(self, gui, settings, canvas):
        self.gui = gui

        ## Settings frame ##
        settingsframe = tk.Frame(settings)
        figcanvas     = tk.Frame(canvas)

        fig_topview = PlotFrame(figcanvas)
        fig_topview.place(relheight=0.5, relwidth=0.5, anchor="nw")

        class PrecheckResult(tk.Text):

            def show(self):
                """
                Run and display results of the preflight tests.
                """
                msg = []
                msg += check_inputs_present(gui.ascot)
                msg += check_options_consistent(gui.ascot)
                if len(msg) == 0:
                    msg += ["Preflight checks completed with no issues."]

                msg0 = ""
                for m in msg:
                    msg0 += m + "\n"

                self.delete("1.0", "end")
                self.insert("end", msg0)


        precheckresults = PrecheckResult(settingsframe, height=25, width=50)
        precheckresults.pack()

        self.settingsframe   = settingsframe
        self.figcanvas       = figcanvas
        self.fig_topview     = fig_topview
        self.precheckresults = precheckresults


    def display(self):

        self.figcanvas.pack_forget()

        self.precheckresults.show()

        self.gui.ascot.init(bfield=True, ignorewarnings=True)
        plot_top_view(self.gui.ascot, axes=self.fig_topview.axis)
        self.fig_topview.draw()

        self.settingsframe.pack()
        self.figcanvas.pack(fill="both", expand=True)
