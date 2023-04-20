import tkinter as tk
from tkinter import ttk

from .components import PlotFrame, NumEntry, DropdownMenu, Tickbox

class DistFrame(ttk.Frame):

    def init(self, gui, canvas):
        self.gui = gui

        class GeneralFrame(ttk.Frame):

            def init(self, canvas):
                f1 = ttk.Frame(self)
                f2 = ttk.Frame(self)

                f1.grid(row=0, column=0, sticky="w")
                f2.grid(row=0, column=1, sticky="ne")

                self.columnconfigure(1, weight=1)

                self.source = DropdownMenu(f1, label="Source: ", width=24,
                                           labelwidth=8, labelanchor="w",
                                           trace=self.setcoords)
                self.xcrd = DropdownMenu(f1, label="x: ", width=5, labelwidth=2,
                                         trace=self.setxcoord)
                self.ycrd = DropdownMenu(f1, label="y: ", width=5, labelwidth=2)
                self.source.grid(row=0, column=0)
                self.xcrd.grid(row=1, column=0, sticky="e")
                self.ycrd.grid(row=2, column=0, sticky="e")

                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                return self

            def setcoords(self, *args):
                dist = self.source.var.get()
                dist = dist[1:-1].split(",")
                coords = []
                for d in dist:
                    coords += [d.strip()]

                self.xcrd.setvals(coords, coords[0])
                self.ycrd.setvals(coords[1:]+["None"], "None")


            def setxcoord(self, *args):
                x = self.xcrd.var.get()
                y = self.ycrd.var.get()

                dist = self.source.var.get()
                dist = dist[1:-1].split(",")
                coords = []
                for d in dist:
                    coords += [d.strip()]

                coords.remove(x)
                if x == y:
                    self.ycrd.setvals(coords, "None")
                else:
                    self.ycrd.setvals(coords, y)


        class IntegratedFrame(ttk.Frame):

            def init(self, canvas):
                f1 = ttk.Frame(self)
                f2 = ttk.Frame(self)

                f1.grid(row=0, column=0, sticky="w")
                f2.grid(row=0, column=1, sticky="ne")

                self.columnconfigure(1, weight=1)

                self.source = DropdownMenu(f1, label="  Source: ", width=24,
                                           labelwidth=8, labelanchor="w")
                self.qnt = DropdownMenu(f1, label="Quantity: ", width=24,
                                        labelwidth=8, labelanchor="w")
                self.source.grid(row=0, column=0)
                self.qnt.grid(row=1, column=0)

                self.plotbutton = tk.Button(f2, text="Plot", width=3)
                self.savebutton = tk.Button(f2, text="Store", width=3)
                self.plotbutton.pack(anchor="nw")
                self.savebutton.pack(anchor="nw")

                self.qnt.setvals(["Density"], "Density")

                return self

        master = ttk.Notebook(self)

        self.framegeneral    = GeneralFrame(master).init(canvas)
        self.frameintegrated = IntegratedFrame(master).init(canvas)
        master.add(self.framegeneral, text="General")
        master.add(self.frameintegrated, text="Integrated")
        master.pack(fill="both", expand=True)

        self.canvas = canvas
        return self


    def display(self):
        run = self.gui.ascot.hdf5.active
        dists = []
        try:
            run.dist5d
            dists += ["{R, phi, z, ppa, ppe}"]
        except AttributeError:
            pass

        try:
            run.distrho5d
            dists += ["{rho, phi, theta, ppa, ppe}"]
        except AttributeError:
            pass

        try:
            run.dist6d
            dists += ["{R, phi, z, pR, pphi, pz}"]
        except AttributeError:
            pass

        try:
            run.distrho6d
            dists += ["{rho, phi, theta, pR, pphi, pz}"]
        except AttributeError:
            pass

        self.framegeneral.source.setvals(dists, dists[0])
        self.frameintegrated.source.setvals(dists, dists[0])


class DistCanvas(ttk.Frame):

    def init(self):
        return self
