"""
Contains definition of OptionsManager class.

File: optionsmanager.py
"""
import tkinter as tk

import numpy as np

class OptionsManager():
    """
    Frame for viewing simulation options.
    """

    def activateoptions(self, frame, canvas, group, ascotpy):
        self.frame   = frame
        self.canvas  = canvas
        self.group   = group
        self.ascotpy = ascotpy


        # Read options data value pairs
        opt, info = self.group.read(info=True)

        # Make sure values are not too accurate or have too many leading zeros.
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

            # Extranem name, description and value
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

        # Make an immutable textbox and display it
        textbox = tk.Text(self.canvas, borderwidth=3, relief="sunken")
        textbox.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        textbox.delete("1.0", tk.END)
        textbox.insert("end", opttext)
        textbox.configure(state="disabled")

        # Expand textbox to fill entire frame
        self.canvas.grid_rowconfigure(0, weight=1)
        self.canvas.grid_columnconfigure(0, weight=1)
