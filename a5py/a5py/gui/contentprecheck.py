import tkinter as tk
import matplotlib.pyplot as plt
from a5py.ascotpy.preflight import *
from tkinter import ttk
from .components import PlotFrame

class ContentPrecheck:

    def __init__(self, gui, settings, canvas):
        self.gui = gui

        class PrecheckFrame(ttk.Frame):

            def init(self):
                self.text = tk.Text(self, height=25, width=50)
                self.text.pack()
                return self

            def view(self):
                """
                Run and display results of the preflight tests.
                """
                msg = []
                msg += check_inputs_present(gui.ascot)
                msg += check_options_consistent(gui.ascot)
                msg += check_bfield_psi0(gui.ascot)
                if len(msg) == 0:
                    msg += ["Preflight checks completed with no issues."]

                msg0 = ""
                for m in msg:
                    msg0 += m + "\n"

                self.text.configure(state='normal')
                self.text.delete("1.0", "end")
                self.text.insert("end", msg0)
                self.text.configure(state='disabled')

        class PrecheckCanvas(ttk.Frame):

            def init(self):
                #self.fig_topview     = PlotFrame(self)
                class PlotFrames(PlotFrame):

                    def set_axes(self):
                        ax1 = self.fig.add_subplot(2,2,(1,3))
                        ax2 = self.fig.add_subplot(2,2,2)
                        ax3 = self.fig.add_subplot(2,2,4)
                        return [ax1, ax2, ax3]


                self.fig_topview     = PlotFrames(self)
                #self.fig_energypitch = PlotFrame(self)
                #self.fig_rhophi      = PlotFrame(self)
                self.fig_topview.place(relheight=0.9, relwidth=0.9, anchor="nw")
                #self.fig_topview.grid(column=0, row=0, rowspan=2)
                #self.fig_energypitch.grid(column=1, row=0)
                #self.fig_rhophi.grid(column=1, row=1)

                return self

            def view(self, ascot):
                self.fig_topview.clear()
                #self.fig_topview.axis.remove()
                plot_top_view(ascot, axes=self.fig_topview.axis[0])
                #self.fig_topview.draw()

                #self.fig_energypitch.clear()
                plot_energypitch(ascot, axes=self.fig_topview.axis[1])
                #self.fig_energypitch.draw()

                #self.fig_rhophi.clear()
                plot_rhophi(ascot, axes=self.fig_topview.axis[2])
                self.fig_topview.draw()
                #self.fig_rhophi.draw()


        frameprecheck = PrecheckFrame(settings).init()
        frameprecheck.pack(fill="both", expand=True)
        canvasprecheck = PrecheckCanvas(canvas).init()
        canvasprecheck.pack(fill="both", expand=True)

        self.canvas         = canvas
        self.frameprecheck  = frameprecheck
        self.canvasprecheck = canvasprecheck


    def display(self):
        self.gui.ascot.init(bfield=True, ignorewarnings=True)

        self.frameprecheck.view()
        self.canvasprecheck.view(self.gui.ascot)
