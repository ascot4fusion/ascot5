"""
Content manager for showing group related information.

This content is shown when user selects a group from the treeview. Therefore,
the content that is shown should be displayed fast. As such, here we show only
things that can be shown without ascotpy.

For all input groups we show editable description.

For parent groups, no input is shown.

Group specific content:
 - options    : show and edit options.
 - marker     : number of markers
 - mhd        : individual modes
 - plasma_1DS : 1D profiles


File: contentgroup.py
"""
import os
import tkinter as tk
import tkinter.ttk as ttk

import numpy as np

from a5py.ascot5io.options     import write_hdf5
from a5py.ascot5io.ascot5tools import call_ascot5file
from .components import PlotFrame

class ContentGroup():

    def __init__(self, gui, settings, canvas):
        """
        Initializes all frames and widgets. Does not display anything yet.
        """
        self.gui = gui

        ## Settings frame ##
        settingsframe = tk.Frame(settings)

        # Load graphics which are shown when no group is selected (on startup)
        logo = os.path.join(os.path.dirname(__file__), "logo.png")
        logo = tk.PhotoImage(file=logo)
        logoframe = tk.Canvas(settingsframe)
        logoframe.create_image(150, 50, image=logo)
        self.logo = logo # Otherwise garbage collector eats this

        # At top we always show box with the description and buttons
        # to save (or revert) the text.
        descframe = tk.Frame(settingsframe)
        f1 = tk.Frame(descframe)
        tk.Label(f1, text="Description:").pack(side="left")
        f1.grid(row=0, column=0, sticky="ew")

        f2 = tk.Frame(descframe)
        button_descsave = tk.Button(f2, text="Save")
        button_descsave.pack(side="right")

        button_descrevert = tk.Button(f2, text="Revert")
        button_descrevert.pack(side="right")

        f2.grid(row=0, column=1, sticky="ew")

        class DescBox(tk.Text):
            """
            Text widget with some helpful methods.
            """

            def revert(self, group):
                """
                Revert description.
                """
                text = group.get_desc()
                self.delete("1.0", "end")
                self.insert("end", text)

            def save(self, group):
                """
                Save description.
                """
                group.set_desc(self.get("1.0", "end-1c"))

        descbox = DescBox(descframe, height=5, width=50)
        descbox.grid(row=1, column=0, columnspan=2, sticky="nsew",
                     padx=2, pady=2)

        ## Options group ##
        optionsframe  = tk.Frame(settingsframe)
        optionscanvas = tk.Frame(canvas)

        # On frame we have save/revert buttons
        tk.Label(optionsframe, text="Options:").pack(side="left")

        button_optionssave = tk.Button(optionsframe, text="Save")
        button_optionssave.pack(side="left")

        button_optionsrevert = tk.Button(optionsframe, text="Revert")
        button_optionsrevert.pack(side="left")

        # On canvas we have scrollable textbox showing options
        class OptionsBox(tk.Text):
            """
            Text box for showing and editing options.
            """

            def highlight(self):
                """
                Highlight values and descriptions.
                """
                lines = self.get('1.0', "end").splitlines()
                for i,line in enumerate(lines):
                    if len(line) < 2:
                        pass
                    elif line[0:2] == "#*":
                        self.tag_add("HEADER",str(i+1)+".0", str(i+1)+".end")
                    elif line[0] == "#":
                        self.tag_add("DESC",str(i+1)+".0", str(i+1)+".end")
                    else:
                        self.tag_add("VALUE",str(i+1)+".0", str(i+1)+".end")

            def view(self, group):
                """
                Read and display options in textbox.
                """
                opttext = group.tostring()
                self.delete("1.0", tk.END)
                self.insert("end", opttext)
                self.highlight()

            def write(self, group):
                """
                Turn options string to a dict and write it.
                """
                lines = self.get("1.0", "end-1c").splitlines()
                opt = {}
                for line in lines:
                    if not line.strip() or line[0] == "#":
                        continue

                    p = line.split("=")
                    if(len(p)==2):
                        opt[p[0].strip()] = np.fromstring(p[1].strip(),sep=",")

                desc = group.get_desc()
                grp = write_hdf5(group._file, opt, desc=desc)
                call_ascot5file(group._file, "set_active", grp)
                gui.files.open_new_file(gui.ascot.h5fn)

        optionsbox = OptionsBox(optionscanvas, borderwidth=3, relief="sunken")
        optionsbox.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        scrollbar_options = ttk.Scrollbar(optionscanvas,
                                          command=optionsbox.yview)
        optionsbox["yscrollcommand"] = scrollbar_options.set
        scrollbar_options.grid(row=0, column=1, sticky='nsew')

        # optionsbox fills the whole canvas
        optionscanvas.grid_rowconfigure(0, weight=1)
        optionscanvas.grid_columnconfigure(0, weight=1)

        mhdframe     = tk.Frame(settingsframe)
        mhdcanvas    = tk.Frame(canvas)
        fig_mhdmodes = PlotFrame(mhdcanvas)
        fig_mhdmodes.place(relheight=0.5, relwidth=0.5, anchor="nw")

        # Styles for highlights
        optionsbox.tag_configure("HEADER", foreground="#004999")
        optionsbox.tag_configure("DESC",   foreground="#26abff")
        optionsbox.tag_configure("VALUE",  foreground="black")

        self.settingsframe        = settingsframe
        self.logoframe            = logoframe
        self.descframe            = descframe
        self.descbox              = descbox
        self.button_descsave      = button_descsave
        self.button_descrevert    = button_descrevert
        self.optionsframe         = optionsframe
        self.optionscanvas        = optionscanvas
        self.optionsbox           = optionsbox
        self.button_optionssave   = button_optionssave
        self.button_optionsrevert = button_optionsrevert
        self.mhdcanvas            = mhdcanvas
        self.mhdframe             = mhdframe
        self.fig_mhdmodes         = fig_mhdmodes


    def display(self, parent, qid):
        """
        Display group information.
        """
        # Clear
        self.logoframe.pack_forget()
        self.descframe.pack_forget()
        self.optionsframe.pack_forget()
        self.optionscanvas.pack_forget()
        self.mhdframe.pack_forget()
        self.mhdcanvas.pack_forget()

        if parent == "" and qid == "":
            # Show logo if no group is selected (on startup/opening new file)
            self.logoframe.pack(fill="both", expand=True)
            self.settingsframe.pack()
            return

        try:
            int(qid)
        except:
            # For parent groups there is nothing to display.
            return

        if parent == "results":
            group = self.gui.ascot.hdf5["q"+qid]
        else:
            group = self.gui.ascot.hdf5[parent]["q"+qid]

        # Button functionality
        self.button_descsave.configure(
            command=lambda:self.descbox.save(group))
        self.button_descrevert.configure(
            command=lambda:self.descbox.revert(group))
        self.descbox.revert(group)
        self.descframe.pack()

        ## GROUP SPECIFIC ##
        if parent == "options":

            self.button_optionssave.configure(
                command=lambda:self.optionsbox.write(group))
            self.button_optionsrevert.configure(
                command=lambda:self.optionsbox.view(group))
            self.optionsbox.view(group)

            self.optionsframe.pack()
            self.optionscanvas.pack(fill="both", expand=True)

        if parent == "mhd":
            self.fig_mhdmodes.clear()
            group.plot_amplitude(ax=self.fig_mhdmodes.axis)
            self.fig_mhdmodes.draw()
            self.mhdcanvas.pack(fill="both", expand=True)
            self.mhdframe.pack()

        # Display
        self.settingsframe.pack()
