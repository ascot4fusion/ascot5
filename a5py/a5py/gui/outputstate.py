import tkinter as tk
from tkinter import ttk

from collections import OrderedDict
from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox, ToggleButton

class StateFrame(ttk.Frame):

    def init(self, gui, canvas):
        self.gui = gui

        quantities = OrderedDict([
            ("R",         "m"),
            ("phi (mod)", "deg"),
            ("z",         "m"),
            ("Time",      "s"),
            ("Energy",    "eV"),
            ("Pitch",     "vpa/|v|"),
            ("ID",        "1"),
            ("x",         "m"),
            ("y",         "m"),
            ("rho",       "1"),
            ("pol (mod)", "deg"),
            ("phi",       "deg"),
            ("pol",       "deg"),
            ("mu",        "eV/T"),
            ("|p|",       "kg m/s"),
            ("ppa",       "kg m/s"),
            ("ppe",       "kg m/s"),
            ("pR",        "kg m/s"),
            ("pphi",      "kg m/s"),
            ("pz",        "kg m/s"),
            ("|B|",       "T"),
            ("BR",        "T"),
            ("Bphi",      "T"),
            ("Bz",        "T"),
            ("charge",    "T")
        ])

        class ScatterFrame(ttk.Frame):

            def init(self, canvas):
                f1 = ttk.Frame(self)
                f2 = ttk.Frame(self)
                f3 = ttk.Frame(self)

                f1.grid(row=0, column=0, sticky="e")
                f2.grid(row=0, column=1, sticky="ne")
                f3.grid(row=1, column=0, sticky="ew")

                self.columnconfigure(1, weight=1)

                self.xcrd = DropdownMenu(f1, log=True, label="x: ", width=12)
                self.ycrd = DropdownMenu(f1, log=True, label="y: ", width=12)
                self.zcrd = DropdownMenu(f1, log=True, label="z: ", width=12)
                self.ccrd = DropdownMenu(f1, log=True, label="c: ", width=12)

                self.xcrd.grid(row=0, column=0)
                self.ycrd.grid(row=1, column=0)
                self.zcrd.grid(row=2, column=0)
                self.ccrd.grid(row=3, column=0)

                self.xbtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.ybtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.zbtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.cbtn = ToggleButton(f1, label1text="Ini", label2text="End")

                self.xbtn.grid(row=0, column=1)
                self.ybtn.grid(row=1, column=1)
                self.zbtn.grid(row=2, column=1)
                self.cbtn.grid(row=3, column=1)

                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                self.endc = DropdownMenu(f3, label="Endcond: ", width=12,
                                         labelwidth=8)
                self.axeq = Tickbox(f3, label=" Axis equal")

                self.endc.pack(side="left", anchor="w")
                self.axeq.pack(side="left", anchor="w")

                self.xcrd.setvals(list(quantities.keys()), "R", 0)
                self.ycrd.setvals(list(quantities.keys()), "z", 0)
                qwithnone = quantities.copy()
                qwithnone["None"] = ("None","")
                self.zcrd.setvals(list(qwithnone.keys()), "None", 0)
                self.ccrd.setvals(list(qwithnone.keys()), "None", 0)

                self.canvas = canvas
                return self

            def plot(self, run):
                xcoord  = self.xcrd.getval()
                ycoord  = self.ycrd.getval()
                zcoord  = self.zcrd.getval()
                ccoord  = self.ccrd.getval()
                equal   = self.axeq.getval() == 1
                xlog    = self.xcrd.islog()
                ylog    = self.ycrd.islog()
                zlog    = self.zcrd.islog()
                clog    = self.ccrd.islog()
                endcond = self.endc.getval()

                if endcond == "all":
                    endcond = None

                log = (xlog==1, ylog==1, zlog==1, clog==1)

                if ccoord == "None":
                    ccoord = None

                iniend = ["e" if self.xbtn.var.get() else "i",
                          "e" if self.ybtn.var.get() else "i",
                          "e" if self.zbtn.var.get() else "i",
                          "e" if self.cbtn.var.get() else "i"]

                if zcoord == "None":
                    self.canvas.fig_rzview.clear()
                    axes = self.canvas.fig_rzview.axis
                    run.plotstate_scatter(xcoord, ycoord, c=ccoord,
                                          axesequal=equal, log=log,
                                          endcond=endcond, axes=axes,
                                          iniend=iniend)
                else:
                    self.canvas.fig_rzview.make3d()
                    self.canvas.fig_rzview.clear()
                    axes = self.canvas.fig_rzview.axis
                    run.plotstate_scatter(xcoord, ycoord, zcoord, ccoord,
                                          axesequal=equal, log=log,
                                          iniend=iniend,
                                          endcond=endcond, axes=axes)
                self.canvas.fig_rzview.draw()


        class HistFrame(ttk.Frame):

            def init(self, canvas):

                f1 = ttk.Frame(self)
                f2 = ttk.Frame(self)
                f4 = tk.Frame(self)
                f3 = ttk.Frame(self)

                f1.grid(row=0, column=0, sticky="e")
                f2.grid(row=0, column=1, sticky="ne")
                f3.grid(row=2, column=0, columnspan=2, sticky="ew")
                f4.grid(row=1, column=0, columnspan=2, sticky="ew")

                self.columnconfigure(1, weight=1)

                self.xcrd = DropdownMenu(f1, log=True, label="x: ", width=18,
                                         trace=self.newxcoord)
                self.ycrd = DropdownMenu(f1, log=True, label="y: ", width=18,
                                         trace=self.newycoord)
                self.xbtn = ToggleButton(f1, label1text="Ini", label2text="End")
                self.ybtn = ToggleButton(f1, label1text="Ini", label2text="End")

                self.xcrd.grid(row=0, column=0)
                self.ycrd.grid(row=1, column=0)

                self.xbtn.grid(row=0, column=1)
                self.ybtn.grid(row=1, column=1)

                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                self.wght = Tickbox(f3, label=" With weights", width=10)
                self.axeq = Tickbox(f3, label=" Axis equal", width=10)

                self.wght.pack(side="left", anchor="w")
                self.axeq.pack(side="left", anchor="w")

                self.xmin_entry = NumEntry(
                    f4, labeltext="rho [1] = ", entrywidth=5, labelwidth=26,
                    anchor="e", defval=0.0)
                self.xmax_entry = NumEntry(f4, labeltext="–", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=1.0)
                self.xnum_entry = NumEntry(f4, labeltext="x", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=50,
                                           isint=True)

                self.ymin_entry = NumEntry(
                    f4, labeltext="* = ", entrywidth=5, labelwidth=26,
                    anchor="e", defval=0.0)
                self.ymax_entry = NumEntry(f4, labeltext="–", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=1.0)
                self.ynum_entry = NumEntry(f4, labeltext="x", entrywidth=5,
                                           labelwidth=2, anchor="c", defval=50,
                                           isint=True)


                self.xmin_entry.grid(row=0, column=0)
                self.xmax_entry.grid(row=0, column=1)
                self.xnum_entry.grid(row=0, column=2)
                self.ymin_entry.grid(row=1, column=0)
                self.ymax_entry.grid(row=1, column=1)
                self.ynum_entry.grid(row=1, column=2)

                f4.columnconfigure(0, weight=1)

                self.xcrd.setvals(list(quantities.keys()), "rho", 0)
                qwithnone = quantities.copy()
                qwithnone["None"] = ("None","")
                self.ycrd.setvals(list(qwithnone.keys()), "None", 0)

                #gui.params.add(
                    #output_stateplot_minr=xmin_entry.choice,
                    #output_stateplot_maxr=xmax_entry.choice,
                    #output_stateplot_numr=xnum_entry.choice,
                    #output_stateplot_minz=ymin_entry.choice,
                    #output_stateplot_maxz=ymax_entry.choice,
                    #output_stateplot_numz=ynum_entry.choice,
                    #output_stateplot_qnt=xcrd.var
                    #)
                self.canvas = canvas

                return self

            def newxcoord(self, *args):
                val = self.xcrd.var.get()
                unit = quantities[val]
                self.xmin_entry.setlabel(val + " [" + unit + "]" + " = ")

            def newycoord(self, *args):
                val = self.ycrd.var.get()
                if val == "None":
                    self.ymin_entry.disable()
                    self.ymax_entry.disable()
                    self.ynum_entry.disable()
                    self.ymin_entry.setlabel("* = ")
                else:
                    unit = quantities[val]
                    self.ymin_entry.enable()
                    self.ymax_entry.enable()
                    self.ynum_entry.enable()
                    self.ymin_entry.setlabel(val + " [" + unit + "]" + " = ")

            def plot(self, run):
                iniend = [self.xbtn.var.get(), self.ybtn.var.get()]
                xcoord = self.xcrd.getval()
                ycoord = self.ycrd.getval()
                if ycoord == "None":
                    ycoord = None

                #endcond = self.endc.getval()
                #if endcond == "all":
                #    endcond = None
                endcond = None

                logx     = self.xcrd.islog()
                logy     = self.ycrd.islog()
                logscale = False#self.zlogtick.getval() == 1

                xbins = [ self.xmin_entry.getval(),
                          self.xmax_entry.getval(),
                          self.xnum_entry.getval() + 1 ]

                ybins = [ self.ymin_entry.getval(),
                          self.ymax_entry.getval(),
                          self.ynum_entry.getval() + 1 ]

                weight = self.wght.getval() == 1

                self.canvas.fig_rzview.clear()
                run.inistate.histogram(x=xcoord, y=ycoord, xbins=xbins,
                                       ybins=ybins, weight=weight, logx=logx,
                                       logy=logy, logscale=logscale,
                                       endcond=endcond, axes=self.canvas.fig_rzview.axis,
                                       iniend=iniend)
                self.canvas.fig_rzview.draw()


        master = ttk.Notebook(self)

        self.framescatter = ScatterFrame(master).init(canvas)
        self.framehist    = HistFrame(master).init(canvas)
        master.add(self.framescatter, text="Scatter")
        master.add(self.framehist, text="Histogram")
        master.pack(fill="both", expand=True)

        self.canvas = canvas
        return self


    def display(self):
        run = self.gui.ascot.hdf5.active
        try:
            endconds, counts = run.endstate.listendconds()
            if self.framescatter.xbtn.switch.isdisabled:
                self.framescatter.xbtn.enable()
                self.framescatter.ybtn.enable()
                self.framescatter.zbtn.enable()
                self.framescatter.cbtn.enable()
                self.framehist.xbtn.enable()
                self.framehist.ybtn.enable()

        except AttributeError:
            endconds, counts = run.inistate.listendconds()
            self.framescatter.xbtn.disable()
            self.framescatter.ybtn.disable()
            self.framescatter.zbtn.disable()
            self.framescatter.cbtn.disable()
            self.framehist.xbtn.disable()
            self.framehist.ybtn.disable()

        endconds = ["all"] + endconds

        self.framescatter.endc.setvals(endconds, "all")

        def plotscatter():
            self.framescatter.plot(run)

        def plothist():
            self.framehist.plot(run)

        self.framescatter.plotbutton.configure(command=plotscatter)
        self.framehist.plotbutton.configure(command=plothist)


class StateCanvas(ttk.Frame):

    def init(self):

        class Plot3DFrame(PlotFrame):

            def __init__(self, master, tight_layout=True):
                self.axis3d = False
                super().__init__(master, tight_layout=tight_layout)

            def make3d(self):
                self.axis3d = True

            def set_axes(self):
                if self.axis3d:
                    self.axis3d = False
                    return self.fig.add_subplot(1,1,1, projection="3d")

                return self.fig.add_subplot(1,1,1)

        fig_rzview = Plot3DFrame(self)
        fig_rzview.place(relheight=0.8, anchor="nw")

        self.fig_rzview = fig_rzview
        return self
