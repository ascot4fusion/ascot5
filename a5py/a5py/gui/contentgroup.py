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
from .components import PlotFrame, ScrollableFrame, DropdownMenu

class ContentGroup():
    """
    Settingsframe is scrollable. On the settingsframe, the first item on top
    is a frame showing mutable description which is always shown when a group
    is selected. When no group is selected (on startup or new file), the logo is
    shown.

    Below description we show group specific information. Contents of the canvas
    depends on what group is selected.
    """

    def __init__(self, gui, settings, canvas):
        """
        Initializes all frames and widgets. Does not display anything yet.
        """
        self.gui = gui

        # Make this frame scrollable
        from .gui import GUI
        settings = ScrollableFrame(settings, framewidth=GUI.FILEFRAMEWIDTH-40)
        settings.pack(fill='both', expand=True)
        settings = settings.scrollable_frame

        ## Logo frame
        class LogoFrame(tk.Canvas):

            def init(self):

                # Load graphics which are shown when no group is selected
                # (on startup)
                logo = os.path.join(os.path.dirname(__file__), "logo.png")
                logo = tk.PhotoImage(file=logo)
                self.logo = logo # Otherwise garbage collector eats this
                self.create_image(150, 50, image=logo)
                return self

        ## Description frame
        class DescFrame(ttk.Frame):

            def init(self):
                # At top we always show box with the description and buttons
                # to save (or revert) the text.
                f1 = tk.Frame(self)
                tk.Label(f1, text="Description:").pack(side="left")
                f1.grid(row=0, column=0, sticky="ew")

                f2 = tk.Frame(self)
                button_descsave = tk.Button(f2, text="Save")
                button_descsave.pack(side="right")

                button_descrevert = tk.Button(f2, text="Revert")
                button_descrevert.pack(side="right")

                f2.grid(row=0, column=1, sticky="ew")

                self.button_descsave   = button_descsave
                self.button_descrevert = button_descrevert

                descbox = tk.Text(self, height=5, width=50)
                descbox.grid(row=1, column=0, columnspan=2, sticky="nsew",
                             padx=2, pady=2)
                self.descbox = descbox
                return self


            def revert(self, group):
                """
                Revert description.
                """
                text = group.get_desc()
                self.descbox.delete("1.0", "end")
                self.descbox.insert("end", text)


            def save(self, group):
                """
                Save description.
                """
                group.set_desc(self.descbox.get("1.0", "end-1c"))

            def view(self, group):
                self.revert(group)
                self.button_descsave.configure(
                    command=lambda:self.save(group))
                self.button_descrevert.configure(
                    command=lambda:self.revert(group))


        ## Options frame
        class OptionsFrame(ttk.Frame):

            def init(self, canvas):

                # On frame we have save/revert buttons
                tk.Label(self, text="Options:").pack(side="left")

                button_optionssave = tk.Button(self, text="Save")
                button_optionssave.pack(side="left")

                button_optionsrevert = tk.Button(self, text="Revert")
                button_optionsrevert.pack(side="left")

                self.button_optionssave   = button_optionssave
                self.button_optionsrevert = button_optionsrevert
                self.canvas = canvas
                return self


            def view(self, group):
                self.button_optionssave.configure(
                command=lambda:self.canvas.write(group))
                self.button_optionsrevert.configure(
                command=lambda:self.canvas.view(group))


        class OptionsCanvas(ttk.Frame):

            def init(self):
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
                return self


            def highlight(self):
                """
                Highlight values and descriptions.
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
                """
                Read and display options in textbox.
                """
                opttext = group.tostring()
                self.text.delete("1.0", tk.END)
                self.text.insert("end", opttext)
                self.highlight()

            def write(self, group):
                """
                Turn options string to a dict and write it.
                """
                lines = self.text.get("1.0", "end-1c").splitlines()
                opt = {}
                for line in lines:
                    if not line.strip() or line[0] == "#":
                        continue

                    p = line.split("=")
                    if(len(p)==2):
                        opt[p[0].strip()] = np.fromstring(p[1].strip(),sep=",")

                desc = group.get_desc()
                grp = write_hdf5(group._file, opt, desc=desc)
                grp = grp.split("_")[-1]
                gui.ascot.hdf5.reload()
                gui.ascot.hdf5["options"]["q"+grp].set_as_active()
                gui.groups.add_group("options", group)


        class MhdFrame(ttk.Frame):

            def init(self, canvas):
                f1 = ttk.Frame(self)
                f2 = ttk.Frame(self)
                f1.pack(side="left")
                f2.pack(side="left")

                ttk.Label(f1, anchor="e", width=10, text="Amplitude:",
                          font=("Calibri 10")).pack()
                ttk.Label(f1, anchor="e", width=10, text="Mode:",
                          font=("Calibri 10")).pack()

                amplitude = DropdownMenu(f2, width=10)
                mode      = DropdownMenu(f2, width=10)

                amplitude.pack()
                mode.pack()

                amplitude.setvals(["Magnetic", "Electric"], "Magnetic")

                self.amplitude = amplitude
                self.mode      = mode
                self.canvas    = canvas

                return self

            def view(self, group):
                data = group.read()
                modes = [None] * int(data["nmode"])
                for i in range(len(modes)):
                    modes[i] = "(" + str(data["nmodes"][i]) \
                               + "," + str(data["mmodes"][i]) + ")"

                self.mode.setvals(["All"] + modes, "All")

                def plot(*args):
                    amplitude = self.amplitude.getval()
                    if amplitude == "Magnetic":
                        amplitude = "alpha"
                    else:
                        amplitude = "phi"

                    mode = self.mode.getval()
                    if mode == "All":
                        mode = None
                    else:
                        mode = mode[1:-1].split(",")# Removes brackets and comma
                        mode = (int(mode[0]), int(mode[1]))

                    self.canvas.view(group, amplitude, mode)

                self.amplitude.settrace(plot)
                self.mode.settrace(plot)
                plot()

        class MhdCanvas(ttk.Frame):

            def init(self):
                self.fig_modes = PlotFrame(self)
                self.fig_modes.place(relheight=0.5, relwidth=0.5, anchor="nw")
                return self

            def view(self, group, amplitude, mode):
                self.fig_modes.clear()
                group.plot_amplitude(amplitude=amplitude,
                                     mode=mode,
                                     ax=self.fig_modes.axes)
                self.fig_modes.draw()

        self.canvas        = canvas
        self.logoframe     = LogoFrame(settings, width=GUI.FILEFRAMEWIDTH-50).init()
        self.descframe     = DescFrame(settings).init()
        self.optionscanvas = OptionsCanvas(canvas).init()
        self.optionsframe  = OptionsFrame(settings).init(self.optionscanvas)
        self.mhdcanvas     = MhdCanvas(canvas).init()
        self.mhdframe      = MhdFrame(settings).init(self.mhdcanvas)


    def display(self):
        """
        Display group information.
        """

        tree   = self.gui.groups.tree
        qid    = tree.item(tree.selection(), "text")
        parent = tree.item(tree.parent(tree.selection()), "text")

        showlogo = False
        if parent == "" and qid == "":
            # Show logo if no group is selected (on startup/opening new file)
            showlogo = True

        try:
            int(qid)
            if parent == "results":
                group = self.gui.ascot.hdf5["q"+qid]
            else:
                group = self.gui.ascot.hdf5[parent]["q"+qid]

        except:
            # For parent groups there is nothing to display.
            showlogo = True

        if showlogo:
            self.descframe.pack_forget()
            self.optionsframe.pack_forget()
            self.optionscanvas.pack_forget()
            self.mhdframe.pack_forget()
            self.mhdcanvas.pack_forget()
            self.logoframe.pack(fill="both", expand=True)

        else:
            self.logoframe.pack_forget()
            self.descframe.pack(fill="both", expand=True)

            self.descframe.view(group)

            if parent == "options":
                self.mhdframe.pack_forget()
                self.mhdcanvas.pack_forget()
                self.optionsframe.pack(fill="both", expand=True)
                self.optionscanvas.pack(fill="both", expand=True)

                self.optionscanvas.view(group)
                self.optionsframe.view(group)

            elif parent == "mhd":
                self.optionsframe.pack_forget()
                self.optionscanvas.pack_forget()
                self.mhdframe.pack(fill="both", expand=True)
                self.mhdcanvas.pack(fill="both", expand=True)

                self.mhdframe.view(group)

            else:
                self.optionsframe.pack_forget()
                self.optionscanvas.pack_forget()
                self.mhdframe.pack_forget()
                self.mhdcanvas.pack_forget()
