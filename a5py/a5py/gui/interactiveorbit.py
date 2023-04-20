import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry
from .livelink import initandrun, clean

class OrbitFrame(ttk.Frame):

    def init(self, canvas):
        f1 = ttk.Frame(self)
        f1.pack()
        r_entry = NumEntry(f1, labeltext="R [m] = ", entrywidth=5,
                           labelwidth=8, anchor="w", defval=6.7)
        r_entry.pack()
        plotbutton = tk.Button(f1, text="Run and plot", width=8)
        plotbutton.pack()

        self.canvas     = canvas
        self.r_entry    = r_entry
        self.plotbutton = plotbutton

        self.first = True

    def display(self, gui):

        def plot():
            r = self.r_entry.getval()
            M = initandrun(gui.ascot.h5fn, r, self.first)
            if self.first:
                self.first = False

            self.canvas.fig_rzview.clear()
            idx=1;
            self.canvas.fig_rzview.axis.plot(M.get_orbit(idx,"r"), M.get_orbit(idx,"z"))
            clean(M)
            self.canvas.fig_rzview.draw()

        self.plotbutton.configure(command=plot)


class OrbitCanvas(ttk.Frame):

    def init(self):
        fig_rzview = PlotFrame(self)
        fig_rzview.place(relheight=0.8, anchor="nw")

        self.fig_rzview = fig_rzview
        return self
