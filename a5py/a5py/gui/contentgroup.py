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
import tkinter as tk
import tkinter.ttk as ttk

import numpy as np

from a5py.ascot5io.options     import write_hdf5
from a5py.ascot5io.ascot5tools import call_ascot5file

class ContentGroup():

    def _display_group(self, gui, settingsframe, canvasframe, parent, qid):
        """
        Display group information.
        """
        if qid is None:
            # For parent groups there is nothing to display.
            return

        if parent == "results":
            group = gui.ascot.hdf5["q"+qid]
        else:
            group = gui.ascot.hdf5[parent]["q"+qid]

        # Always show description box on top.
        f1 = tk.Frame(settingsframe)
        tk.Label(f1, text="Description:").pack(side="left")
        f1.grid(row=0, column=0, sticky="ew")

        f2 = tk.Frame(settingsframe)
        save = tk.Button(f2, text="Save")
        save.pack(side="right")

        load = tk.Button(f2, text="Revert")
        load.pack(side="right")

        f2.grid(row=0, column=1, sticky="ew")

        descbox = tk.Text(settingsframe, height=5, width=50)
        descbox.grid(row=1, column=0, columnspan=2, sticky="nsew",
                     padx=2, pady=2)
        descbox.insert("end", group.get_desc())

        # Helper function
        def replacetext(text):
            descbox.delete("1.0", "end")
            descbox.insert("end", text)

        # Button functionality
        save.configure(
            command=lambda:group.set_desc(descbox.get("1.0", "end-1c")))
        load.configure(
            command=lambda:replacetext(group.get_desc()))

        # Add extra frame where widgets can be added
        frame = tk.Frame(settingsframe)
        frame.grid(row=2, column=0, columnspan=2, sticky="nsew")


        ## GROUP SPECIFIC ##

        if parent == "options":
            contentgroup_options(
                frame,canvasframe, group,
                lambda : gui.files.open_new_file(gui.ascot.h5fn),
                lambda : descbox.get("1.0", "end-1c")
            )


def contentgroup_options(frame, canvas, group, filechanged,
                         getdescinbox):
    """
    Display options data for a given group.
    """

    ## Add Widgets ##

    # Save and load buttons
    tk.Label(frame, text="Options:").pack(side="left")

    save = tk.Button(frame, text="Save")
    save.pack(side="left")

    load = tk.Button(frame, text="Revert")
    load.pack(side="left")

    # Make a textbox and display it
    textbox = tk.Text(canvas, borderwidth=3, relief="sunken")
    textbox.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

    scrollb = ttk.Scrollbar(canvas, command=textbox.yview)
    scrollb.grid(row=0, column=1, sticky='nsew')
    textbox["yscrollcommand"] = scrollb.set

    # Expand textbox to fill entire frame
    canvas.grid_rowconfigure(0, weight=1)
    canvas.grid_columnconfigure(0, weight=1)

    # Styles for highlights
    textbox.tag_configure("HEADER", foreground="#004999")
    textbox.tag_configure("DESC", foreground="#26abff")
    textbox.tag_configure("VALUE", foreground="black")

    ## Set functionality ##

    def applyhighlights():
        """
        Highlight values and descriptions.
        """
        lines = textbox.get('1.0', "end").splitlines()
        for i,line in enumerate(lines):
            if len(line) < 2:
                pass
            elif line[0:2] == "#*":
                textbox.tag_add("HEADER",str(i+1)+".0", str(i+1)+".end")
            elif line[0] == "#":
                textbox.tag_add("DESC",str(i+1)+".0", str(i+1)+".end")
            else:
                textbox.tag_add("VALUE",str(i+1)+".0", str(i+1)+".end")

    def viewoptions():
        """
        Read and display options in textbox.
        """
        opttext = group.tostring()
        textbox.delete("1.0", tk.END)
        textbox.insert("end", opttext)
        applyhighlights()

    def writeoptions():
        """
        Helper function that turns options string to a dict and writes it.
        """
        lines = textbox.get("1.0", "end-1c").splitlines()
        opt = {}
        for line in lines:
            if not line.strip() or line[0] == "#":
                continue

            p = line.split("=")
            if(len(p)==2):
                opt[p[0].strip()] = np.fromstring(p[1].strip(),sep=",")

        desc = getdescinbox()
        grp = write_hdf5(self.group._file, opt, desc=desc)
        call_ascot5file(self.group._file, "set_active", grp)
        filechanged()


    load.configure(command=viewoptions)
    save.configure(command=writeoptions)

    # Show options
    viewoptions()
