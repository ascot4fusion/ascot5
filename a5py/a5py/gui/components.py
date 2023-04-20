"""
Contains definition of various GUI component classes.

File: components.py
"""
import tkinter
import tkinter.ttk as ttk

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

class NumEntry(tkinter.Frame):
    """
    A component to receive numerical input.
    """

    def __init__(self, master, isint=False, labeltext=None, defval=None,
                 entrywidth=6, labelwidth=8, anchor="w"):
        """
        Initialize a frame which accepts inputs.

        The frame consists of label and entrybox side-by-side. The entry only
        accepts valid numerical input. The value is obtained with getval()
        method.
        """
        super().__init__(master)
        self.isint = isint

        self.choice = tkinter.StringVar(self)

        if defval is not None:
            self.choice.set(defval)

        label = tkinter.Label(self, text=labeltext, anchor=anchor,
                              width=labelwidth, font=("Calibri 10"))

        vcmd  = self.register(self._validate)
        entry = tkinter.Entry(self, validate = "all",
                              validatecommand=(vcmd, "%P"), width=entrywidth,
                              textvariable=self.choice, font=("Calibri 10"))

        label.grid(row=0, column=0, sticky="W")
        entry.grid(row=0, column=1, sticky="E")


    def getval(self):
        """
        Parse entry to a valid float or int.
        """
        val = self.choice.get()

        val = val.split("e")

        s1 = val[0]
        if s1 in ["", ".", "+", "-"]:
            s1 = 0
        else:
            if self.isint:
                s1 = int(s1)
            else:
                s1 = float(s1)

        if len(val) == 2:
            s2 = val[1]
            if s2 in ["", "+", "-"]:
                s2 = 0
            else:
                s2 = int(s2)

            s1 = s1 * np.power(10.0, s2)

        if self.isint:
            return int(s1)
        else:
            return float(s1)


    def isempty(self):
        """
        Check if entry is just an empty string.
        """
        val = self.choice.get()
        return val == ""


    def _validate(self, P):
        """
        Check if entry is a valid float or int.
        """

        # Check if there is (exactly one) exponent
        if P.count("e") > 1:
            return False
        P = P.split("e")

        # Check that the prefix is valid
        s = P[0]
        if not ( (s == "" or s in "+-") or (not self.isint and s in ".") ) :
            try:
                if self.isint:
                    int(s)
                else:
                    float(s)
            except ValueError:
                return False

        # Check that exponent is valid
        if len(P) == 2:
            s = P[1]
            if not ( (s == "" or s in "+") or (not self.isint and s in "-") ):
                try:
                    int(s)
                except ValueError:
                    return False

        return True


class DropdownMenu(tkinter.Frame):
    """
    Dropdown menu where user can choose from given options.
    """

    def __init__(self, master, defval=None, values=None, log=False, trace=None,
                 label=None, width=6, labelwidth=2):
        super().__init__(master)
        self.var = tkinter.StringVar(self)

        if defval is not None:
            self.var.set(defval)

        if trace is not None:
            self.var.trace('w', trace)

        if label is not None:
            label = tkinter.Label(self, text=label, width=labelwidth)
            label.grid(row=0, column=0)

        self.menu = ttk.Combobox(self, width=width, textvariable=self.var,
                                 state="readonly")
        self.menu.grid(row=0, column=1)
        if values is not None:
            self.menu["values"] = values

        self.log = None
        if log:
            self.log = tkinter.IntVar(self)
            logtick = tkinter.Checkbutton(
                self, text=" log10", onvalue=1,
                offvalue=0, variable=self.log, height=1, width=5)
            logtick.grid(row=0, column=2)

            if trace is not None:
                self.log.trace('w', trace)


    def islog(self):
        if self.log is None:
            return False

        return self.log.get()


    def getval(self):
        return self.var.get()

    def setvals(self, values, defval, log=None):
        self.menu["values"] = values
        self.var.set(defval)
        if log is not None:
            self.log.set(log)


    def settrace(self, trace):
        self.var.trace('w', trace)
        if self.log is not None:
            self.log.trace('w', trace)


class Tickbox(tkinter.Frame):
    """
    A tickbox and label.
    """

    def __init__(self, master, defval=None, trace=None, label=None, width=8):
        super().__init__(master)
        self.var = tkinter.IntVar(self)
        if defval is not None:
            self.var.set(defval)

        if trace is not None:
            self.var.trace('w', trace)

        tick = tkinter.Checkbutton(self, text=label,
                                   variable=self.var,
                                   onvalue=1, offvalue=0,
                                   height=1, width=width)
        tick.pack()


    def getval(self):
        return self.var.get()


class NavToolbarNocoord(NavigationToolbar2Tk):
    """
    Navigation toolbar but without the coordinate display
    """
    def set_message(self, msg):
        pass


class PlotFrame(tkinter.Frame):
    """
    Frame containing matplotlib plot and NavToolbarNocoord
    """

    def __init__(self, master, tight_layout=True):
        super().__init__(master)
        fig = plt.figure(tight_layout=tight_layout)
        figcanvas = FigureCanvasTkAgg(fig, self)

        toolbar = NavToolbarNocoord(figcanvas, self, pack_toolbar=False)
        toolbar.config(background="white")
        for button in toolbar.winfo_children():
            button.config(background="white")
        toolbar.update()

        toolbar.pack(side=tkinter.TOP, fill=tkinter.X)
        figcanvas.get_tk_widget().pack(fill=tkinter.BOTH, expand=True)

        self.fig       = fig
        self.axis      = self.set_axes()
        self.figcanvas = figcanvas

    def draw(self):
        self.figcanvas.draw()

    def clear(self):
        if self.fig.axes:
            for ax in self.fig.axes:
                self.fig.delaxes(ax)

        self.axis = self.set_axes()

    def set_axes(self):
        return self.fig.add_subplot(1,1,1)


class ScrollableFrame(ttk.Frame):

    def __init__(self, container, *args, framewidth=100, **kwargs):
        super().__init__(container, *args, **kwargs)
        canvas = tkinter.Canvas(self, width=framewidth)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")


class ToggleButton(tkinter.Canvas):
    def __init__(self, *args, command=None, fg='black', bg='gray', width=35,
                 height=18, **kwargs):
        super().__init__(*args, **kwargs)

        self.configure(width=width, height=height, borderwidth=0,
                       highlightthickness=0)

        self.back_ground = self.create_arc((0, 0, 0, 0), start=90, extent=180,
                                           fill="green", outline='')
        self.back_ground1 = self.create_arc((0, 0, 0, 0), start=-90, extent=180,
                                            fill="red", outline='')
        self.rect = self.create_rectangle(0, 0, 0, 0, fill=bg, outline='')

        self.btn = self.create_oval(0, 0, 0, 0, fill=fg, outline='')

        self.bind('<Configure>', self._resize)
        self.bind('<Button>', self._animate, add='+')
        self.bind('<Button>', command, add='+')

        self.state = 0

    def _resize(self, event):
        self.coords(self.back_ground, 5, 5, event.height-5, event.height-5)
        self.coords(self.back_ground1, 5, 5, event.height, event.height-5)

        factor = event.width-(self.coords(self.back_ground1)[2] -
                              self.coords(self.back_ground1)[0])-10
        self.move(self.back_ground1, factor, 0)

        self.coords(self.rect, self.bbox(self.back_ground)[2]-2, 5,
                    self.bbox(self.back_ground1)[0]+2, event.height-5)

        self.coords(self.btn, 5, 5, event.height-5, event.height-5)

        if self.state:
            self.moveto(self.btn, self.coords(self.back_ground1)[0]+4, 4)

    def _animate(self, event):
        x, y, w, h = self.coords(self.btn)
        x = int(x-1)
        y = int(y-1)

        if x == self.coords(self.back_ground1)[0]+3:
            self.moveto(self.btn, 4, 4)
            self.state = 0

        else:
            self.moveto(self.btn, self.coords(self.back_ground1)[0]+4, 4)
            self.state = 1

    def get_state(self):
        return self.state
