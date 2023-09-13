"""Contains definition of various GUI component classes.
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

        self.label = tkinter.Label(self, text=labeltext, anchor=anchor,
                                   width=labelwidth, font=("Calibri 10"))

        vcmd  = self.register(self._validate)
        self.entry = tkinter.Entry(
            self, validate = "all", validatecommand=(vcmd, "%P"),
            width=entrywidth, textvariable=self.choice, font=("Calibri 10"))

        self.label.grid(row=0, column=0, sticky="W")
        self.entry.grid(row=0, column=1, sticky="E")


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


    def setval(self, val):
        self.choice.set(val)


    def isempty(self):
        """
        Check if entry is just an empty string.
        """
        val = self.choice.get()
        return val == ""


    def disable(self):
        self.label.configure(fg="gray")
        self.entry.config(state="disabled")


    def enable(self):
        self.label.configure(fg="black")
        self.entry.config(state="normal")

    def setlabel(self, text):
        self.label.configure(text=text)


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
                 label=None, width=6, labelwidth=2, labelanchor="c"):
        super().__init__(master)
        self.var = tkinter.StringVar(self)

        if defval is not None:
            self.var.set(defval)

        if trace is not None:
            self.var.trace('w', trace)

        if label is not None:
            label = tkinter.Label(self, text=label, width=labelwidth)
            label.pack(side="left", anchor=labelanchor)

        self.menu = ttk.Combobox(self, width=width, textvariable=self.var,
                                 state="readonly")
        self.menu.pack(side="left")
        if values is not None:
            self.menu["values"] = values

        self.log = None
        if log:
            self.log = tkinter.IntVar(self)
            logtick = tkinter.Checkbutton(
                self, text=" log10", onvalue=1,
                offvalue=0, variable=self.log, height=1, width=5)
            logtick.pack(side="left")

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


class Tickbox(tkinter.Checkbutton):
    """
    A tickbox and label.
    """

    def __init__(self, master, defval=None, trace=None, label=None, width=8,
                 justify=None):
        self.var = tkinter.IntVar()
        super().__init__(master, text=label, variable=self.var,
                         onvalue=1, offvalue=0, height=1, width=width,
                         justify=justify)
        if defval is not None:
            self.var.set(defval)

        if trace is not None:
            self.var.trace('w', trace)


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
        self.axes      = self.set_axes()
        self.figcanvas = figcanvas

    def draw(self):
        self.figcanvas.draw()

    def clear(self):
        if self.fig.axes:
            for ax in self.fig.axes:
                self.fig.delaxes(ax)

        self.axes = self.set_axes()

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


class Switch(tkinter.Canvas):
    def __init__(self, *args, command=None, fg='black', bg='gray', width=35,
                     height=18, **kwargs):
        super().__init__(*args, **kwargs)

        self.configure(width=width, height=height, borderwidth=0,
                       highlightthickness=0)

        self.back_ground = self.create_arc((0, 0, 0, 0), start=90,
                                           extent=180, fill=bg, outline='')
        self.back_ground1 = self.create_arc((0, 0, 0, 0), start=270,
                                            extent=180, fill=bg, outline='')
        self.rect = self.create_rectangle(0, 0, 0, 0, fill=bg, outline='')

        self.btn = self.create_oval(0, 0, 0, 0, fill=fg, outline='')

        self.bind('<Configure>', self._resize)
        self.bind('<Button>', self._animate, add='+')
        self.bind('<Button>', command, add='+')

        self.var = tkinter.IntVar(self)
        self.var.set(0)
        self.isdisabled = False
        self.fillcolor = bg


    def _resize(self, event):
        self.coords(self.back_ground, 5, 5, event.height-5, event.height-5)
        self.coords(self.back_ground1, 5, 5, event.height-5, event.height-5)

        factor = event.width-(self.coords(self.back_ground1)[2] -
                              self.coords(self.back_ground1)[0])-10
        self.move(self.back_ground1, factor, 0)

        self.coords(self.rect, self.bbox(self.back_ground)[2]-1, 6,
                    self.bbox(self.back_ground1)[0]+2, event.height-5)

        self.coords(self.btn, 4, 5, event.height-6, event.height-5)

        if self.var.get():
            self.moveto(self.btn, self.coords(self.back_ground1)[0]+1, 4)


    def _animate(self, event):
        if self.isdisabled:
            return

        x, y, w, h = self.coords(self.btn)
        x = int(x-1)
        y = int(y-1)

        if x == self.coords(self.back_ground1)[0]:
            self.toggle("off")

        else:
            self.toggle("on")


    def toggle(self, onoff="on"):
        if onoff == "off":
            self.moveto(self.btn, 4, 4)
            self.var.set(0)
        elif onoff == "on":
            self.moveto(self.btn, self.coords(self.back_ground1)[0]+1, 4)
            self.var.set(1)

    def get_state(self):
        return self.var.get()


    def disable(self, onoff="off"):
        self.toggle(onoff)
        self.isdisabled = True
        self.itemconfigure(self.rect, fill="red")
        self.itemconfigure(self.back_ground, fill="red")
        self.itemconfigure(self.back_ground1, fill="red")


    def enable(self, onoff="off"):
        self.toggle(onoff)
        self.isdisabled = False
        self.itemconfigure(self.rect, fill=self.fillcolor)
        self.itemconfigure(self.back_ground, fill=self.fillcolor)
        self.itemconfigure(self.back_ground1, fill=self.fillcolor)


class ToggleButton(ttk.Frame):

    def __init__(self, *args, command=None, fg='black', bg='gray', width=35,
                 height=18, label1text="Off", label2text="On", **kwargs):
        super().__init__(*args, **kwargs)
        self.switch = Switch(self, command=command, fg=fg, bg=bg, width=width,
                             height=height)

        self.labeloff = tkinter.Label(self, text=label1text)
        self.labelon  = tkinter.Label(self, text=label2text)

        self.labeloff.pack(side="left")
        self.switch.pack(side="left")
        self.labelon.pack(side="left")
        self.var = self.switch.var


    def disable(self, onoff="off"):
        self.switch.disable(onoff)
        if onoff == "off":
            self.labelon.configure(fg="gray")
        elif onoff == "on":
            self.labeloff.configure(fg="gray")


    def enable(self, onoff="off"):
        self.switch.enable(onoff)
        self.labelon.configure(fg="black")
        self.labelon.configure(fg="black")

class ContentTab(ttk.Frame):
    """Leaf node in the NestedNotebook that is a frame widget.
    """

    def __init__(self, frame, canvas, gui, *args, tabselected=None, tabdisabled=None, **kwargs):
        super().__init__(frame, *args, tabselected=tabselected, tabdisabled=tabdisabled, **kwargs)
        self.frame  = frame
        self.canvas = canvas
        self.gui    = gui

class NestedNotebook(ttk.Notebook):
    """Notebook that can has other notebooks or frames in it's tabs.
    """

    def __init__(self, frame, *args, tabselected=None, **kwargs):
        super().__init__(frame, *args, **kwargs)
        self._children = {}
        self._sleep = True

        if tabselected is None:
            self.tabselected = lambda : 0
        else:
            self.tabselected = tabselected

        def eventhandler(event):
            #print(event.widget.tab("current")["text"])
            self.tabchanged()
        self.bind("<<NotebookTabChanged>>", eventhandler)

    def add(self, text, tabselected=None, tab=None):
        if tab == None:
            tab = NestedNotebook(self, tabselected=tabselected)
        self._children[text] = tab
        super().add(tab, text=text)

    def traverse(self, text):
        for txt, frame in self._children.items():
            if txt == text: return frame
        for c in self._children.values():
            if not isinstance(c, NestedNotebook): continue
            frame = c.traverse(text)
            if frame is not None: return frame
        return None

    def currenttab(self):
        tab = self.nametowidget(self.select())
        if isinstance(tab, NestedNotebook):
            return tab.currenttab()
        for name, item in self._children.items():
            if item == tab:
                return tab, name

    def tabchanged(self):
        if self._sleep: return
        tab = self.nametowidget(self.select())
        if isinstance(tab, NestedNotebook):
            tab.tabselected()
            tab.tabchanged()
        else:
            tab.selecttab()

    def wakeup(self):
        self._sleep = False
        for c in self._children.values():
            if not isinstance(c, NestedNotebook): continue
            c.wakeup()
