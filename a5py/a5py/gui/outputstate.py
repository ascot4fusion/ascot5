import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class StateFrame(ttk.Frame):

    def init(self):

        class ScatterFrame(ttk.Frame):

            def init(self):
                xcrd = DropdownMenu(self, log=True, label="x: ")
                ycrd = DropdownMenu(self, log=True, label="y: ")
                zcrd = DropdownMenu(self, log=True, label="z: ")
                ccrd = DropdownMenu(self, log=True, label="c: ")
                endc = DropdownMenu(self, label="Endcond: ")
                axeq = Tickbox(self, label="Axis equal")

                xcrd.pack()
                ycrd.pack()
                zcrd.pack()
                ccrd.pack()
                endc.pack()
                axeq.pack()

                return self


        class HistFrame(ttk.Frame):

            def init(self):
                return self


        master = ttk.Notebook(self)
        master.pack(fill="both", expand=True)

        scatterframe = ScatterFrame(master).init()
        histframe    = HistFrame(master).init()
        master.add(scatterframe, text="Scatter")
        master.add(histframe, text="Histogram")

        return self
