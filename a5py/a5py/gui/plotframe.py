"""
Contains definition of PlotFrame class.

File: plotframe.py
"""
import tkinter
import importlib.util as util

plt = util.find_spec("matplotlib")
if plt:
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

#from .gui import GUIMINWIDTH
#from .gui import GUIMINHEIGHT

from . import gui as guimod

# Size specifiers
INDEXPANELHEIGHT = 100
INDEXPANELWIDTH  = 250


class PlotFrame(tkinter.Frame):
    """
    An opening frame where other frames can be accessed and HDF5 file modified.
    """

    def __init__(self, gui):
        """
        Initialize index.
        """
        super().__init__(gui._root)
        self._gui = gui
        self._fig = plt.figure()

        #toprow = tkinter.Frame(self, height=100)
        #botrow = tkinter.Frame(self)

        indexpanel = tkinter.Frame(self,
                                   width=INDEXPANELWIDTH,
                                   height=INDEXPANELHEIGHT)
        toppanel   = tkinter.Frame(self)
        sidepanel  = tkinter.Frame(self)
        backbutton  = tkinter.Button(indexpanel, text="Back")
        canvas      = FigureCanvasTkAgg(self._fig, master=self)

        # Navigation toolbars but without the coordinate display
        class Toolbar(NavigationToolbar2Tk):
            def set_message(self, msg):
                pass
        navbar = Toolbar(canvas, indexpanel)

        backbutton.config(command=self._backtoindex)

        backbutton.pack()
        navbar.pack()

        indexpanel.grid(row=0, column=0)
        toppanel.grid(row=0, column=1, sticky="NEWS")
        sidepanel.grid(row=1, column=0, sticky="NEWS")
        canvas.get_tk_widget().grid(row=1, column=1, sticky="NEWS",
                                    padx=10, pady=10)

        indexpanel.pack_propagate(0)
        toppanel.grid_propagate(0)
        sidepanel.grid_propagate(0)
        canvas.get_tk_widget().pack_propagate(0)

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._toppanel  = toppanel
        self._sidepanel = sidepanel
        self._canvas    = canvas
        self.draw()


    def _backtoindex(self):
        from .indexframe import IndexFrame
        self._gui.displayframe(IndexFrame(self._gui))


    def get_fig(self):
        self._fig.clear()
        return self._fig


    def get_toppanel(self):
        self._toppanel.destroy()
        toppanel = tkinter.Frame(self)
        toppanel.grid(row=0, column=1, sticky="NEWS")
        toppanel.pack_propagate(0)
        toppanel.grid_propagate(0)
        self._toppanel = toppanel
        return self._toppanel


    def get_sidepanel(self):
        self._sidepanel.destroy()
        sidepanel = tkinter.Frame(self)
        sidepanel.grid(row=1, column=0, sticky="NEWS")
        sidepanel.pack_propagate(0)
        sidepanel.grid_propagate(0)
        self._sidepanel = sidepanel
        return self._sidepanel


    def draw(self):
        self._canvas.draw()
