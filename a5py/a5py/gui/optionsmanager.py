"""
Contains definition of OptionsManager class.

File: optionsmanager.py
"""
import tkinter as tk
import tkinter.ttk as ttk

import numpy as np

from a5py.ascot5io.options     import write_hdf5
from a5py.ascot5io.ascot5tools import call_ascot5file

class OptionsManager():
    """
    Frame for viewing simulation options.
    """

    def display_options(self, frame, canvas, group, filechanged, getdescinbox):
        """
        Display options data for a given group.

        Requires sidepanel (frame), canvas (canvas), the options group (group)
        and function that updates GUI when file has changed (filechanged). The
        last one is required to allow user to save options in GUI.
        """
        self.frame   = frame
        self.canvas  = canvas
        self.group   = group


        def readoptions():
            """
            Helper function that reads options data to a single string
            """
            # Read options data value pairs
            opt, info = self.group.read(info=True)

            # Make sure values are not too accurate
            # or have too many leading zeros.
            def trim(val):
                if np.abs(val) >= 1e4 or (np.abs(val) <=1e-4 and val != 0):
                    return "{0:e}".format(val)
                else:
                    return "{0:g}".format(val)

            # Store options in a single string variable
            opttext = ""
            for i in info:
                if len(i) == 1:
                    # Subtitle
                    opttext = "\n" + opttext + i[0] + "\n\n"
                    continue

                # Extract name, description and value
                name        = i[0]
                description = i[1]
                value       = opt[name]

                opttext = opttext + description
                if(len(value) == 1):
                    # Single value
                    opttext = opttext + name + "=" + trim(value[0]) + "\n\n"
                elif(len(value) > 1):
                    # Multiple values, show as a list
                    opttext = opttext + name + "=" \
                            + ",".join([trim(v) for v in value]) \
                            + "\n\n"

            # Remove empty lines from beginning and end
            opttext = opttext[5:-2]

            return opttext


        ## Add Widgets ##

        # Save and load buttons
        tk.Label(self.frame, text="Options:").pack(side="left")

        save = tk.Button(self.frame, text="Save")
        save.pack(side="left")

        load = tk.Button(self.frame, text="Revert")
        load.pack(side="left")

        # Make a textbox and display it
        textbox = tk.Text(self.canvas, borderwidth=3, relief="sunken")
        textbox.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)

        scrollb = ttk.Scrollbar(self.canvas, command=textbox.yview)
        scrollb.grid(row=0, column=1, sticky='nsew')
        textbox["yscrollcommand"] = scrollb.set

        # Expand textbox to fill entire frame
        self.canvas.grid_rowconfigure(0, weight=1)
        self.canvas.grid_columnconfigure(0, weight=1)

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
            opttext = readoptions()
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
