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
            ("|v|",       "m/s"),
            ("vpa",       "m/s"),
            ("vpe",       "m/s"),
            ("vR",        "m/s"),
            ("vphi",      "m/s"),
            ("vz",        "m/s"),
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

                return self


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

                wght = Tickbox(f3, label=" With weights", width=10)
                axeq = Tickbox(f3, label=" Axis equal", width=10)

                wght.pack(side="left", anchor="w")
                axeq.pack(side="left", anchor="w")

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
                self.framehistogram.xbtn.enable()
                self.framehistogram.ybtn.enable()

        except AttributeError:
            endconds, counts = run.inistate.listendconds()
            self.framescatter.xbtn.disable()
            self.framescatter.ybtn.disable()
            self.framescatter.zbtn.disable()
            self.framescatter.cbtn.disable()
            self.framehistogram.xbtn.disable()
            self.framehistogram.ybtn.disable()

        endconds = ["all"] + endconds

        self.framescatter.endc.setvals(endconds, "all")


class StateCanvas(ttk.Frame):

    def init(self):
        return self
