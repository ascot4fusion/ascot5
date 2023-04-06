import os
import sys
import subprocess
import tkinter as tk
from tkinter import ttk

from tkinter.filedialog import askopenfilename, askdirectory
from tkinter import messagebox

import a5py.ascot5io.ascot5tools as tools
from a5py.ascot5io.ascot5 import Ascot, AscotInitException
from a5py.ascotpy.libascot import AscotpyInitException
from a5py.ascotpy.ascotpy import Ascotpy
from a5py.ascot5io.ascot5file import INPUT_PARENTS

from .contentmanager import ContentManager
from .guifile import Datalink

class GUI(tk.Tk):
    """
    Window where all GUI elements are stored.


    The GUI looks like following:
     ____________________________________________
    |  FileFrame   |                             |
    |______________|                             |
    |  GroupFrame  ^                             |
    | (Scrollable) |        CanvasFrame          |
    | (Resizable)  |        (Resizable)          |
    |              |                             |
    |              |                             |
    |              |                             |
    |______________v                             |
    |              ^                             |
    | SettingsFrame|                             |
    | (Scrollable) |                             |
    | (Resizable)  |                             |
    |______________v_____________________________|

    FileFrame has a pop-up window for choosing a file to be opened.

    GroupFrame also has an pop-up window for operations related to
    HDF5 file contents (like activating or removing a group).

    FileFrame and GroupFrame are always displayed. Contents of
    SettingsFrame and CanvasFrame depends on what is clicked on the
    FileFrame/GroupFrame/SettingsFrame and those are managed by
    the ContentManager class which this GUI inherits.

    Access to the data is handled by a separate class.
    """

    # Minimum GUI size in pixels
    MINWIDTH  = 800
    MINHEIGHT = 600

    # FileFrame size (fixed) in pixels. Width is also used in other frames in
    # the same column.
    FILEFRAMEWIDTH  = 350
    FILEFRAMEHEIGHT = 100

    # Other frame heights
    INPUTFRAMEHEIGHT    = 250
    SETTINGSFRAMEHEIGHT = 250

    # Screenwidth method won't regonize dual monitors but instead returns their
    # combined width. If the screenwidth is above this limit (in pixels), assume
    # dual monitors are in use.
    DUALMONITORWIDTH = 2000

    BORDERWIDTH = 2


    def __init__(self, filename=None):
        tk.Tk.__init__(self)

        ## Add decorations: title and icons ##
        self.title("ASCOT5 GUI")
        icon = os.path.join(os.path.dirname(__file__), "icon.png")
        self.tk.call("wm", "iconphoto", self._w, tk.PhotoImage(file=icon))

        ## Set window size and minimum size ##
        sw = self.winfo_screenwidth()
        sh = self.winfo_screenheight()

        # If the screen width is very large it might be due to dual
        # monitors. Divide the screen width by 2 to make GUI look nice
        # with dual monitors.
        if sw > GUI.DUALMONITORWIDTH:
            sw /= 2

        # GUI dimensions and locations
        w = sw*3/4
        h = sh*3/4
        x = (sw/2) - (w/2)
        y = (sh/2) - (h/2)

        # Make sure dimensions are below minimum values
        if w < GUI.MINWIDTH:
            w = GUI.MINWIDTH
        if h < GUI.MINHEIGHT:
            h = GUI.MINHEIGHT

        # Set window size and set minimum size
        self.geometry('%dx%d+%d+%d' % (w, h, x, y))
        self.minsize(GUI.MINWIDTH, GUI.MINHEIGHT)

        # Make sure everything is closed when the window is closed
        def close():
            """
            Close the gui and terminate the program.
            """
            self.destroy()
            exit()

        self.protocol("WM_DELETE_WINDOW", close)


        # Initialize content frames
        files    = FileFrame(self,     width=GUI.FILEFRAMEWIDTH,
                                       height=GUI.FILEFRAMEHEIGHT,
                                       borderwidth=GUI.BORDERWIDTH)
        groups   = GroupFrame(self,    width=GUI.FILEFRAMEWIDTH,
                                       height=GUI.INPUTFRAMEHEIGHT,
                                       borderwidth=GUI.BORDERWIDTH)
        settings = SettingsFrame(self, width=GUI.FILEFRAMEWIDTH,
                                       height=GUI.SETTINGSFRAMEHEIGHT,
                                       borderwidth=GUI.BORDERWIDTH)
        canvas   = CanvasFrame(self,   width=GUI.MINWIDTH-GUI.FILEFRAMEWIDTH,
                                       height=GUI.MINHEIGHT-GUI.FILEFRAMEHEIGHT,
                               bg="white")

        files.grid(   row=0, column=0, sticky="NWES")
        groups.grid(  row=1, column=0, sticky="NWES")
        settings.grid(row=2, column=0, sticky="NWES")
        canvas.grid(  row=0, column=1, sticky="NWSE", rowspan=3)

        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_columnconfigure(1, weight=1)

        self.files    = files
        self.groups   = groups

        # Set up content manager
        self.contentmanager = ContentManager(self, settings.get_frame(), canvas)

        # Read file and show its contents
        self.ascot = Ascotpy()
        path = None if filename is None else os.path.abspath(filename)
        self.files.open_new_file(path)


class SettingsFrame(ttk.Frame):
    """
    Scrollable frame which can contain buttons and text fields.

    Contents of this frame changes depending on what is being plotted or what is
    selected. This frame can be thought of as a settings panel for CanvasFrame.
    """

    SCROLLBARWIDTH = 36

    def __init__(self, container, *args, **kwargs):
        """
        Initializes a frame with a scrollbar.

        Use contentframe as a root for child widgets.

        See: https://blog.tecladocode.com/tkinter-scrollable-frames/
        """
        super().__init__(container, *args, **kwargs)
        canvas = tk.Canvas(self,
                           width=kwargs["width"]-SettingsFrame.SCROLLBARWIDTH,
                           height=kwargs["height"])
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.contentframe = ttk.Frame(canvas)

        self.contentframe.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.contentframe, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")


    def get_frame(self):
        """
        Return the actual contents frame.
        """
        return self.contentframe


class CanvasFrame(tk.Frame):
    """
    Frame that displays plots.
    """

    def __init__(self, container, *args, **kwargs):
        tk.Frame.__init__(self, container, *args, **kwargs)


class FileFrame(tk.Frame):
    """
    Frame for accessing the HDF5 file.
    """

    def __init__(self, gui, *args, **kwargs):
        tk.Frame.__init__(self, gui, *args, **kwargs)
        self.gui = gui

        frame1 = tk.Frame(self)
        frame2 = tk.Frame(self)
        frame3 = tk.Frame(self, width=345, height=100)

        # Label "File: XXX"
        tk.Label(frame1, text="File: ").pack(side="left", anchor="nw")

        # Textbox for showing the current file
        self.filenamefield = tk.Text(frame1, width=30, height=1)
        self.filenamefield.pack(side="left", anchor="nw")

        # Browse button. We need to do this frame trick to set it correct size
        tempframe = tk.Frame(frame1, width=24, height=22)
        tempframe.pack_propagate(False)
        self.browsebutton = tk.Button(tempframe, text="...")
        self.browsebutton.pack(side="left")
        tempframe.pack(side="left", anchor="nw")

        # Label "Size: XXX Mb"
        self.filesizelabel = tk.Label(frame2, text="Size: ")
        self.filesizelabel.pack(side="left", anchor="nw")

        # Sanity check button
        self.sanitybutton = tk.Button(frame3, text="Sanity checks")
        self.sanitybutton.pack(side="left", anchor="nw")

        frame1.pack(side="top", fill="x",    anchor="nw")
        frame2.pack(side="top", fill="x",    anchor="nw")
        frame3.pack(side="top", fill="both", anchor="nw")


        ## Set functionality ##

        # Browse button opens filename dialog when clicked
        self.browsebutton.configure(command=self.browse_file)

        # Unmutable filename field
        self.filenamefield.configure(state="disabled", wrap="none")


    def browse_file(self):
        """
        Open dialog for choosing HDF5 file and open it.

        If a valid file is chosen, deliver the filename to GUI and refresh.
        """
        fn = askopenfilename( title="Select ASCOT5 HDF5 file",
                              filetypes = [("HDF5 files","*.h5")] )
        if len(fn) == 0:
            # Do nothing
            return

        # Get filename and refresh
        fn = os.path.abspath(fn)
        self.open_new_file(fn)


    def open_new_file(self, filename):
        """
        Update the contents of the fileframe if file has changed.
        """
        if filename is not None:

            try:
                self.gui.ascot.reload(filename)
            except AscotInitException:
                messagebox.showerror(
                    "Error",
                    "Could not open file. The file could be corrupted.\n" +
                    "Try Ascot(fn.h5) in Python or h5ls fn.h5 in terminal.")
                return

            except AscotpyInitException:
                messagebox.showwarning(
                    "Warning",
                    "Could not initialize Ascotpy.\n" +
                    "Check that libascot is compiled and included in " +
                    "LD_LIBRARY_PATH")
                return

        self.gui.groups.build_groups()
        self.display_file(filename)


    def display_file(self, filename):
        """
        Display file name and size.
        """
        if filename is None:
            filename = ""
            size = 0.0
        else:
            size = os.path.getsize(filename) / 1e6 # B -> MB

        # Update filename
        self.filenamefield.configure(state="normal")
        self.filenamefield.delete("1.0", tk.END)
        self.filenamefield.insert("end", filename)
        self.filenamefield.configure(state="disabled")
        self.filenamefield.see(tk.END)

        # Update size
        self.filesizelabel.configure(
            text="Size: " + "{0:.3f}".format(size) + " MB")


class GroupFrame(tk.Frame):
    """
    Frame for displaying groups in a treeview.

    Also handles any operations regarding contents of the file.
    """

    def __init__(self, gui, *args, **kwargs):
        """
        Initializes tree which doesn't contain any data.
        """
        super().__init__(*args, **kwargs)

        # Set white background, header font, and default font for tree entries.
        style = ttk.Style()
        style.configure("Treeview.Heading", font=("Calibri", 10),
                        background="white")
        style.map('Treeview.Heading', background=[('focus', 'white')])
        style.configure("Treeview", font=("Calibri", 9),
                        background="white", fieldbackground="white", indent=1)

        # Tree where only one item can be selected at a time
        tree = ttk.Treeview(self, style="mystyle.Treeview", selectmode="browse")

        # Three columns: Group, Type, and Date
        tree["columns"]=("#1","#2")
        tree.column("#0", width= 90, minwidth= 90, stretch=False)
        tree.column("#1", width= 85, minwidth= 85, stretch=False)
        tree.column("#2", width=130, minwidth=130, stretch=True)
        tree.heading("#0", text="Group")
        tree.heading("#1", text="Type")
        tree.heading("#2", text="Date")

        # Styles for different entries:
        # - Every odd entry has a white background and every even a darker one
        # - Parents have white background with bold and larger font.
        # - If input parent don't exist in the file the font is set to red
        # - Active group has a green bolded medium font.
        tree.tag_configure("even",    background="#EFEFEF")
        tree.tag_configure("odd",     background="white")
        tree.tag_configure("nodata",  foreground="red")
        tree.tag_configure("parent",  background="white",
                           font=("Calibri", 10, "bold"))
        tree.tag_configure("active",  foreground="green",
                           font=("Calibri", 9, "bold"))

        # Prevents user from resizing columns
        tree.bind('<Motion>', 'break')
        tree.pack(side=tk.TOP, expand=True, fill="both")

        # Set actions
        tree.bind("<<TreeviewSelect>>", self.select_group)

        # Initialize popupmenu
        groupmenu = self.GroupMenu(gui, self, tree)

        self.tree = tree
        self.gui  = gui


    def build_groups(self):
        """
        (Re)Build tree from scratch from the file in GUI.
        """

        # Clear tree.
        self.tree.delete(*self.tree.get_children())

        # Check if we have an open HDF5 file and proceed.
        if self.gui.ascot.get_filepath() is None:
            return

        # Construct the tree by creating all parents and their children
        parents = self.gui.ascot.hdf5.get_parents()
        for index, parent in enumerate(list(parents.keys()) + ["results"]):

            if parent == "results":
                qids, types, descs, dates = \
                    self.gui.ascot.hdf5.get_resultsinfo(sortbydate=True)
                if len(qids) > 0:
                    activeqid = self.gui.ascot.hdf5.active.get_qid()
                    item = self.tree.insert("", index, text=parent,
                                            tags=("parent"))
            else:
                if parents[parent] is None:
                    # The parent does not exists
                    item = self.tree.insert("", index, text=parent,
                                            tags=("parent", "nodata"))
                    continue

                item = self.tree.insert("", index, text=parent, tags=("parent"))
                qids, types, descs, dates = \
                    self.gui.ascot.hdf5[parent].get_inputinfo(sortbydate=True)
                if len(qids) > 0:
                    activeqid = self.gui.ascot.hdf5[parent].active.get_qid()

            # Generate children:
            # - Rows are tagged "even" or "odd"
            # - Active groups are tagged "active"
            alternatingrowcolor = 0
            for i in range(len(qids)):
                tags = []
                tags.append("odd" if alternatingrowcolor % 2 else "even")
                alternatingrowcolor += 1
                if qids[i] == activeqid:
                    tags.append("active")

                self.tree.insert(item, "end", text=qids[i], tags=tuple(tags),
                                 values=(types[i], dates[i]))


    def select_group(self, event):
        """
        User selected a group.

        Update tree view accordingly and notify GUI.
        """
        qid    = self.tree.item(self.tree.selection(), "text")
        parent = self.tree.item(self.tree.parent(self.tree.selection()), "text")

        # Test if a valid QID or a parent group.
        try:
            int(qid)
        except:
            parent = qid
            qid    = None

        # If selection is a run, highlight the inputs it has used. If selection
        # is a input, highlight the runs that have used it.
        if qid is not None:
            if parent == "results":
                inputqids = []
                for p in INPUT_PARENTS:
                    if p in self.gui.ascot.hdf5["q"+qid]:
                        inputqids.append(
                            self.gui.ascot.hdf5["q"+qid][p].get_qid())

                self._highlightinputs(inputqids)

            else:
                outputqids = self.gui.ascot.hdf5.get_runsfrominput(qid)
                self._highlightoutputs(outputqids)

            # Show contents for this group
            self.gui.contentmanager.display_group(parent, qid)


    def activate_group(self):
        """
        Set item (group) as active in the tree.
        """
        item       = self.tree.selection()
        parent     = self.tree.parent(item)
        qid        = self.tree.item(item, "text")
        parentname = self.tree.item(parent, "text")

        if parentname == "results":
            self.gui.ascot.hdf5["q"+qid].set_as_active()
        else:
            self.gui.ascot.hdf5[parentname]["q"+qid].set_as_active()
        self._activate(item, parent)


    def remove_group(self):
        """
        Remove item (group) from tree.
        """
        item       = self.tree.selection()
        parent     = self.tree.parent(item)
        parentname = self.tree.item(parent, "text")
        qid        = self.tree.item(item, "text")

        # Really remove?
        mbox = messagebox.askquestion (
            "Remove group",
            "Are you sure you want to remove " + qid
            + " permanently from the HDF5 file?",
            icon = "warning")

        if mbox != "yes":
            return

        # Remove the group in HDF5
        self.gui.ascot.hdf5._remove_from_file(qid)

        try:
            self.gui.ascot.hdf5[parentname]["q" + qid]
            messagebox.showerror(
                  "Error",
                  "Could not remove the input group as it has been used in\n"
                + "a simulation. Remove the simulation group first.")
            return
        except:
            pass

        try:
            self.gui.ascot.hdf5[qid]
            messagebox.showerror(
                  "Error",
                  "Could not remove the input group as it has been used in\n"
                + "a simulation. Remove the simulation group first.")
            return
        except:
            pass

        # Update treeview

        # 1. Whole group was removed / ceased to exist
        if qid in INPUT_PARENTS or qid == "results":
            for c in self.tree.get_children(item):
                self.tree.delete(c)

            if qid == "results":
                self.tree.delete(item)
                self._highlightinputs([])
            else:
                self.tree.item(item, tags=("parent", "nodata"))
            return

        # 2. Ordinary group: remove it and set the next one as active
        if parentname == "results":
            activeqid = self.gui.ascot.hdf5.active.get_qid()
            self._highlightinputs([])
        else:
            activeqid = self.gui.ascot.hdf5[parentname].active.get_qid()

        self.tree.delete(item)
        for group in self.tree.get_children(parent):
            if activeqid == self.tree.item(group, "text"):
                self._activate(group, parent)
                break

    def export_group(self):
        """
        Export given group to another HDF5 file.

        Open a prompt for querying destination file.

        This operation requires no GUI re-initialization.
        """
        tree  = self.tree
        item  = tree.selection()
        group = tree.item(tree.selection(), "text")
        fnout = askopenfilename( title="Select ASCOT5 HDF5 file",
                                 filetypes = [("HDF5 files","*.h5")] )
        if len(fnout) == 0:
            # Do nothing
            return

        # Check that the target file is valid.
        try:
            Ascot(fnout)
        except:
            messagebox.showerror(
                    "Error",
                    "Could not open the file.")
            return

        try:
            tools.copygroup(self.gui.ascot.hdf5._hdf5fn, fnout, group)

        except Exception as err:
            messagebox.showerror(
                    "Error",
                    "Could not copy the group as something went wrong.")
            return


    def add_group(self, dummy=False):
        """
        Add a dummy group to a given parent group.

        The GUI is re-initialized after this operation.
        """
        tree   = self.groups.tree
        item   = tree.selection()
        parent = tree.item(tree.selection(), "text")
        #write_dummy_input(self.filename, parent, desc="Dummy")

        #self.ascot   = Ascot(self.filename)
        #self.ascotpy = Ascotpy(self.filename)
        #self.groups.build_groups(self.filename)
        #self.contentmanager.clear()


    def _highlightinputs(self, inputqids):
        """
        Highlight inputs used in the selected run.
        """
        children = self.tree.get_children()
        for c in children:
            if self.tree.item(c, "text") == "results":
                continue

            groups = self.tree.get_children(c)
            for g in groups:
                qid = self.tree.item(g, "text")
                vals = self.tree.item(g, "values")
                date = vals[1]
                if date[0] == "*":
                    date = date[1:-1]

                if qid in inputqids:
                    date = "*" + date + "*"

                self.tree.item(g, values=(vals[0], date))


    def _highlightoutputs(self, outputqids):
        """
        Highlight inputs used in the selected run.
        """
        children = self.tree.get_children()
        for c in children:
            if self.tree.item(c, "text") != "results":
                continue

            groups = self.tree.get_children(c)
            for g in groups:
                qid = self.tree.item(g, "text")
                vals = self.tree.item(g, "values")
                date = vals[1]
                if date[0] == "*":
                    date = date[1:-1]

                if qid in outputqids:
                    date = "*" + date + "*"

                self.tree.item(g, values=(vals[0], date))

    def _activate(self, item, parent):
        """
        Move active tag from the previous active to the current one.
        """
        for c in self.tree.get_children(parent):
            tags = list(self.tree.item(c, "tags"))

            if "active" in tags:
                tags.remove("active")
                self.tree.item(c, tags=tuple(tags))
            if c == item[0] or c == item:
                tags.append("active")
                self.tree.item(c, tags=tuple(tags))


    class GroupMenu(tk.Menu):
        """
        Popupmenu displayed when groups treeview is right-clicked.
        """

        def __init__(self, gui, groups, tree, *args, **kwargs):
            super().__init__(gui, *args, tearoff=0, **kwargs)

            self.add_command(label="Activate",
                             command=groups.activate_group)
            self.add_command(label="Remove",
                             command=groups.remove_group)
            self.add_command(label="Export",
                             command=groups.export_group)
            self.add_command(label="Add dummy input",
                             command=lambda : groups.add_group(dummy=True))

            # Show menu when tree is right-clicked
            tree.bind("<Button-3>", self.showmenu)

            # Hide menu when focus is lost or user clicks anywhere else
            self.bind("<FocusOut>", self.hidemenu)
            gui.bind("<Button-1>", self.hidemenu)

            self.groups = groups


        def showmenu(self, event):
            """
            Show menu if a group was right-clicked and select the group.
            """
            item       = self.groups.tree.identify("item", event.x, event.y)
            itemname   = self.groups.tree.item(item, "text")
            parent     = self.groups.tree.parent(item)
            parentname = self.groups.tree.item(parent, "text")
            if len(itemname) > 0:
                self.post(event.x_root, event.y_root)
                self.focus_set()

                # Disable "activate" option if selected item is a parent
                self.entryconfigure(0, state="normal")
                if itemname in INPUT_PARENTS + ["results"]:
                    self.entryconfigure(0, state="disabled")

                # Disable "remove" option if parent has no children
                self.entryconfigure(1, state="normal")
                if itemname in INPUT_PARENTS and \
                   len(self.groups.tree.get_children(item)) == 0:
                    self.entryconfigure(1, state="disabled")

                # Disable "export" option for parents and runs
                self.entryconfigure(2, state="normal")
                if itemname in INPUT_PARENTS + ["results"] \
                   or parentname == "results":
                    self.entryconfigure(2, state="disabled")

                # Disable "add dummy input" option for other than input parents
                self.entryconfigure(4, state="normal")
                if itemname not in INPUT_PARENTS:
                    self.entryconfigure(4, state="disabled")

                # This will launch a event that will update the GUI.
                self.groups.tree.selection_set(item)


        def hidemenu(self, event):
            """
            Hide menu.
            """
            self.unpost()
