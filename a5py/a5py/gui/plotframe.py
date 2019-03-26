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

class PlotFrame(tkinter.Frame):
    """
    An opening frame where other frames can be accessed and HDF5 file modified.
    """

    def __init__(self, gui):
        """
        Initialize index .
        """
        super().__init__(gui._root)
        self._gui = gui
        self._fig = plt.figure()

        #toprow = tkinter.Frame(self, height=100)
        #botrow = tkinter.Frame(self)

        indexpanel = tkinter.Frame(self, width=250, height=100)
        toppanel   = tkinter.Frame(self, width=gui.width-250, height=100)
        sidepanel  = tkinter.Frame(self, height=gui.height-100, width=200)
        backbutton = tkinter.Button(indexpanel, text="Back")
        canvas     = FigureCanvasTkAgg(self._fig, master=self)

        # Navigation toolbars but without the coordinate display
        class Toolbar(NavigationToolbar2Tk):
            def set_message(self, msg):
                pass
        navbar = Toolbar(canvas, indexpanel)

        backbutton.config(command=self._backtoindex)

        backbutton.pack()
        navbar.pack()

        indexpanel.grid(row=0, column=0)
        toppanel.grid(row=0, column=1)
        sidepanel.grid(row=1, column=0)
        canvas.get_tk_widget().grid(row=1, column=1, sticky="NEWS",
                                    padx=10, pady=10)

        indexpanel.pack_propagate(0)
        toppanel.grid_propagate(0)
        sidepanel.grid_propagate(0)
        canvas.get_tk_widget().pack_propagate(0)

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
        toppanel = tkinter.Frame(self, width=self._gui.width-200, height=100)
        toppanel.grid(row=0, column=1)
        toppanel.pack_propagate(0)
        toppanel.grid_propagate(0)
        self._toppanel = toppanel
        return self._toppanel

    def get_sidepanel(self):
        self._sidepanel.destroy()
        sidepanel = tkinter.Frame(self, width=200, height=self._gui.height-100)
        sidepanel.grid(row=1, column=0)
        sidepanel.pack_propagate(0)
        sidepanel.grid_propagate(0)
        self._sidepanel = sidepanel
        return self._sidepanel

    def draw(self):
        self._canvas.draw()
