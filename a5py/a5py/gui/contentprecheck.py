import tkinter as tk
import matplotlib.pyplot as plt
from a5py.ascotpy.preflight import *
from tkinter import ttk
from .components import PlotFrame

class ContentPrecheck:

    def display_precheck(self, gui, frame, canvas):

        ## First add widgets ##

        fig_topview = PlotFrame(canvas)
        fig_topview.place(relheight=0.5, relwidth=0.5, anchor="nw")

        text_msg = tk.Text(frame, height=25, width=50)
        text_msg.pack()

        ## Run preflight checks ##
        msg = []
        msg += check_inputs_present(gui.ascot)
        msg += check_options_consistent(gui.ascot)
        if len(msg) == 0:
            msg += ["Preflight checks completed with no issues."]

        msg0 = ""
        for m in msg:
            msg0 += m + "\n"

        ## Display checks ##

        text_msg.insert("end", msg0)

        gui.ascot.init(bfield=True)
        plot_top_view(gui.ascot, axes=fig_topview.axis)
        fig_topview.draw()
