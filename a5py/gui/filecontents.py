import os
import tkinter as tk
import tkinter.ttk as ttk
from importlib.resources import files as imfiles

import numpy as np

from a5py.ascot5io.options import Opt
from .components import ContentTab, PlotFrame, ScrollableFrame, DropdownMenu

class Info(ContentTab):

    def __init__(self, frame, canvas, gui, *args, **kwargs):
        super().__init__(frame, canvas, gui, *args, **kwargs)
        emptycanvas   = ttk.Frame(self.canvas)
        optionscanvas = Info.OptionsCanvas(self.canvas)

        self.logo = Info.Logo(self)
        self.logo.grid(column=0, row=0, sticky="WE")

        self.desc = Info.Desc(self)
        self.desc.grid(column=0, row=1, sticky="WE")

        self.opt = Info.Options(self, optionscanvas, gui)
        self.opt.grid(column=0, row=2, sticky="WE")

        self.canvas.slideadd(self.logo, emptycanvas)
        self.canvas.slideadd(self.opt, optionscanvas)

        self.logo.grid_forget()
        self.desc.grid_forget()
        self.opt.grid_forget()

    def selecttab(self):
        tree   = self.gui.groups.tree
        qid    = tree.item(tree.selection(), "text")
        parent = tree.item(tree.parent(tree.selection()), "text")

        self.logo.grid_forget()
        self.desc.grid_forget()
        self.opt.grid_forget()

        if parent == "" and qid == "":
            # Show logo if no group is selected (on startup/opening new file)
            self.logo.grid()
            self.canvas.slideshow(self.logo)
            return

        try:
            int(qid)
            if parent == "results":
                group = self.gui.ascot.data["q"+qid]
            else:
                group = self.gui.ascot.data[parent]["q"+qid]
        except ValueError:
            # For parent groups there is nothing to display.
            self.canvas.slideshow(self.logo)
            return

        # Show group description and related information
        self.desc.grid()
        self.desc.view(group)

        if parent == "options":
            self.canvas.slideshow(self.opt)
            self.opt.grid()
            self.opt.view(group)
            return

    class Logo(tk.Canvas):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            logo = imfiles('a5py.gui').joinpath('logo.png')
            logo = tk.PhotoImage(file=logo)
            self.logo = logo # Otherwise garbage collector eats this
            self.create_image(150, 50, image=logo)

    class Desc(ttk.Frame):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            # At top we always show box with the description and buttons
            # to save (or revert) the text.
            f1 = ttk.Frame(self)
            ttk.Label(f1, text="Description:").pack(side="left")
            f1.grid(row=0, column=0, sticky="ew")

            f2 = ttk.Frame(self)
            button_descsave = ttk.Button(f2, text="Save")
            button_descsave.pack(side="right")

            button_descrevert = ttk.Button(f2, text="Revert")
            button_descrevert.pack(side="right")

            f2.grid(row=0, column=1, sticky="ew")

            self.button_descsave   = button_descsave
            self.button_descrevert = button_descrevert

            descbox = tk.Text(self, height=5, width=50)
            descbox.grid(row=1, column=0, columnspan=2, sticky="nsew",
                         padx=2, pady=2)
            self.descbox = descbox

        def revert(self, group):
            """Revert description.
            """
            text = group.get_desc()
            self.descbox.delete("1.0", "end")
            self.descbox.insert("end", text)

        def save(self, group):
            """Save description.
            """
            group.set_desc(self.descbox.get("1.0", "end-1c"))

        def view(self, group):
            self.revert(group)
            self.button_descsave.configure(
                command=lambda:self.save(group))
            self.button_descrevert.configure(
                command=lambda:self.revert(group))

    class Options(ttk.Frame):

        def __init__(self, frame, canvas, gui, *args, **kwargs):
            super().__init__(frame, *args, **kwargs)
            self.gui    = gui
            self.canvas = canvas

            # On frame we have save/revert buttons
            tk.Label(self, text="Options:").pack(side="left")

            self.savebutton = tk.Button(self, text="Save")
            self.savebutton.pack(side="left")

            self.revertbutton = tk.Button(self, text="Revert")
            self.revertbutton.pack(side="left")

        def view(self, group):
            self.savebutton.configure(
                command=lambda:self.canvas.write(group, self.gui))
            self.revertbutton.configure(
                command=lambda:self.canvas.view(group))
            self.canvas.view(group)

    class OptionsCanvas(ttk.Frame):

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            textbox = tk.Text(self, borderwidth=3, relief="sunken")
            textbox.pack(side="left", fill="both", expand=True)

            scrollbar = ttk.Scrollbar(self, command=textbox.yview)
            textbox["yscrollcommand"] = scrollbar.set
            scrollbar.pack(side="left", fill="y")

            # Styles for highlights
            textbox.tag_configure("HEADER", foreground="#004999")
            textbox.tag_configure("DESC",   foreground="#26abff")
            textbox.tag_configure("VALUE",  foreground="black")

            self.text = textbox
            self.text.delete("1.0", tk.END)

        def highlight(self):
            """Highlight values and descriptions.
            """
            lines = self.text.get('1.0', "end").splitlines()
            for i,line in enumerate(lines):
                if len(line) < 2:
                    pass
                elif line[0:2] == "#*":
                    self.text.tag_add("HEADER",str(i+1)+".0",
                                      str(i+1)+".end")
                elif line[0] == "#":
                    self.text.tag_add("DESC",str(i+1)+".0",
                                      str(i+1)+".end")
                else:
                    self.text.tag_add("VALUE",str(i+1)+".0",
                                      str(i+1)+".end")

        def view(self, group):
            """Read and display options in textbox.
            """
            opttext = group.tostring()
            self.text.delete("1.0", tk.END)
            self.text.insert("end", opttext)
            self.highlight()

        def write(self, group, gui):
            """Turn options string to a dict and write it.
            """
            lines = self.text.get("1.0", "end-1c").splitlines()
            opt = Opt.convert_string(lines)

            desc = group.get_desc()
            new = gui.ascot.data.create_input("opt", **opt, desc=desc,
                                              activate=True)
            gui.groups.add_group("options", gui.ascot.data._get_group(new))
