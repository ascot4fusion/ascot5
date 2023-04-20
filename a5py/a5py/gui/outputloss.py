import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class LossFrame(ttk.Frame):

    def init(self, gui, canvas):

        class GeneralFrame(ttk.Frame):

            def init(self, canvas):
                return self

        class WallloadFrame(ttk.Frame):

            def init(self, canvas):
                return self

        class View3DFrame(ttk.Frame):

            def init(self, canvas):
                return self

        master = ttk.Notebook(self)

        framegeneral  = GeneralFrame(master).init(canvas)
        framewallload = WallloadFrame(master).init(canvas)
        frameview3d   = View3DFrame(master).init(canvas)
        master.add(framegeneral,  text="General")
        master.add(framewallload, text="Wall loads")
        master.add(frameview3d,   text="View 3D")
        master.pack(fill="both", expand=True)

        self.canvas = canvas

        return self


    def display(self):
        pass


class LossCanvas(ttk.Frame):

    def init(self):
        return self
