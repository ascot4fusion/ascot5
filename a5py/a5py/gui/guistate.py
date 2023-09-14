import os
import sys
import subprocess
import tkinter as tk
from tkinter import ttk

from tkinter.filedialog import askopenfilename, askdirectory
from tkinter import messagebox

import a5py.ascot5io.ascot5tools as tools
from a5py.ascot5io.ascot5  import Ascot, write_dummy as write_dummy_input
from a5py.ascotpy          import Ascotpy
from a5py.ascot5io.ascot5file import INPUT_PARENTS

from .contentmanager import ContentManager

class GUI(tk.Tk):
    """
    Window where all GUI elements are stored and which keeps track of its state.


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

    GUI organizes these frames and keeps a track of the content of
    InputFrame, and ResultsFrame (i.e. which groups are selected) so that the
    selection remains when window is refereshed when the HDF5 file is altered.

    GUI also provides access to Ascot and Ascotpy objects.
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

        # Initialize popupmenu
        groupmenu = GroupMenu(self, self.groups.tree)

        # Set up content manager
        self.contentmanager = ContentManager(self, settings.get_frame(), canvas)

        # Read file and show its contents
        self.filename = None if filename is None else os.path.abspath(filename)
        self.ascot    = None
        self.ascotpy  = None
        self.filechanged(self.filename)
        self.files.filechanged(self.filename)


    def filechanged(self, filename):
        """
        Update Ascot, Ascotpy, and groups frame when file has changed.
        """
        if filename is not None:
            try:
                self.ascot   = Ascot(filename)
            except:
                raise AscotInitException

            try:
                self.ascotpy = Ascotpy(filename)
            except:
                raise AscotpyInitException

            self.filename = filename

        self.groups.init(filename)


    def selectionchanged(self, parent, qid):
        """
        Update settings and canvas frames to correspond to current selection.
        """
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
                    if p in self.ascot["q"+qid]:
                        inputqids.append(self.ascot["q"+qid][p].get_qid())

                self.groups.highlightinputs(inputqids)

            else:
                outputqids = self.ascot.get_runsfrominput(qid)
                self.groups.highlightoutputs(outputqids)

        self.contentmanager.selectionchanged(
            parent, qid, self.ascot, self.ascotpy)


    def ascotfile_activate(self):
        """
        Set selected group as active in the HDF5 file.
        """
        tree  = self.groups.tree
        item  = tree.selection()
        group = tree.item(tree.selection(),"text")
        tools.call_ascot5file(self.filename, "set_active", group)
        self.groups.activate_item(item)


    def ascotfile_remove(self):
        """
        Remove given group from the HDF5 file.

        First launch a prompt to make sure user understands what is going on and
        then remove. Launch another prompt if the removal was unsuccessful.
        Freeze GUI until the operation is done (it takes a while to
        hdf5_repack).

        The GUI is re-initialized after this operation.
        """
        tree  = self.groups.tree
        item  = tree.selection()
        group = tree.item(tree.selection(),"text")
        mbox = messagebox.askquestion (
            "Remove group",
            "Are you sure you want to remove " + group
            + " permanently from the HDF5 file?",
            icon = "warning")

        if mbox == "yes":
            try:
                tools.removegroup(self.filename, group)
            except RuntimeError:
                messagebox.showerror(
                    "Error",
                    "Could not remove the group. (You will see this message" +
                    " if you are removing an input which was used in a run.)\n"
                    + "Try using a5removegroup directly")
                return
            try:
                rein  = self.filename
                reout = self.filename + "_repack"
                subprocess.call(["h5repack", rein, reout])
                subprocess.call(["mv", reout, rein])
            except:
                messagebox.showerror(
                    "Error",
                    "Could not repack the HDF5 file. (h5repack failed)\n"
                    + "The group was removed but disk space was not freed.")
                return

        self.ascot   = Ascot(self.filename)
        self.ascotpy = Ascotpy(self.filename)
        self.groups.remove_item(item, self.filename)
        self.files.filechanged(self.filename)
        self.contentmanager.clear()


    def ascotfile_export(self):
        """
        Export given group to another HDF5 file.

        Open a prompt for querying destination file.

        This operation requires no GUI re-initialization.
        """
        tree  = self.groups.tree
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
            tools.copygroup(self.filename, fnout, group)

        except:
            messagebox.showerror(
                    "Error",
                    "Could not copy the group as something went wrong.")
            return


    def ascotfile_adddummy(self):
        """
        Add a dummy group to a given parent group.

        The GUI is re-initialized after this operation.
        """
        tree   = self.groups.tree
        item   = tree.selection()
        parent = tree.item(tree.selection(), "text")
        write_dummy_input(self.filename, parent, desc="Dummy")

        self.ascot   = Ascot(self.filename)
        self.ascotpy = Ascotpy(self.filename)
        self.groups.init(self.filename)
        self.contentmanager.clear()


class FileFrame(tk.Frame):
    """
    Frame for accessing the HDF5 file.
    """

    def __init__(self, container, *args, **kwargs):
        tk.Frame.__init__(self, container, *args, **kwargs)

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
        self.browsebutton.configure(command=self.browsefile)

        # Unmutable filename field
        self.filenamefield.configure(state="disabled", wrap="none")


    def browsefile(self):
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
        try:
            self.master.filechanged(fn)
            self.filechanged(fn)

        except AscotInitException:
            # Prompt corrupted file dialog.
            messagebox.showerror(
                "Error",
                "Could not open file. The file could be corrupted.\n" +
                "Try Ascot(fn.h5) in Python or h5ls fn.h5 in terminal.")

        except AscotpyInitException:
            # Prompt libascot not found dialog.
            messagebox.showwarning(
                "Warning",
                "Could not initialize Ascotpy.\n" +
                "Check that libascot is compiled and included in " +
                "LD_LIBRARY_PATH")


    def filechanged(self, filename):
        """
        Update the contents of the fileframe if file has changed.
        """
        if filename is not None:

            # Update filename
            self.filenamefield.configure(state="normal")
            self.filenamefield.delete("1.0", tk.END)
            self.filenamefield.insert("end", filename)
            self.filenamefield.configure(state="disabled")
            self.filenamefield.see(tk.END)

            # Update size
            size = os.path.getsize(filename) / 1e6 # B -> MB
            self.filesizelabel.configure(
                text="Size: " + "{0:.3f}".format(size) + " MB")


class GroupFrame(tk.Frame):
    """
    Frame for displaying groups in a treeview.
    """

    def __init__(self, container, *args, **kwargs):
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
        tree.column("#0", width=100, minwidth=100, stretch=False)
        tree.column("#1", width=125, minwidth=125, stretch=False)
        tree.column("#2", width=110, minwidth=110, stretch=False)
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


        tree.pack(side=tk.TOP, expand=True, fill="both")
        self.tree = tree

        def selectionchanged(event):
            """
            Notify MainWindow when the selection changes.
            """
            selection = self.tree.selection()
            parent    = self.tree.parent(selection)
            container.selectionchanged(
                self.tree.item(parent,    "text"),
                self.tree.item(selection, "text") )

        self.tree.bind("<<TreeviewSelect>>", selectionchanged)


    def init(self, ascotfn):
        """
        Set frame to its default state.
        """

        # Clear tree.
        self.tree.delete(*self.tree.get_children())

        # Check if we have an open HDF5 file and proceed.
        if ascotfn is None:
            return
        ascot = Ascot(ascotfn)

        # Construct the tree by creating all input parents and their children
        parents = INPUT_PARENTS + ["results"]
        for parent in parents:
            index = parents.index(parent)

            # Fetch all QIDs for this parent children and see which is active
            try:
                qids = tools.call_ascot5file(ascotfn, "get_qids", parent)
            except:
                # The parent does not exists
                if parent != "results":
                    item = self.tree.insert("", index, text=parent,
                                            tags=("parent", "nodata"))
                continue

            item = self.tree.insert("", index, text=parent, tags=("parent"))
            activeqid = tools.call_ascot5file(
                ascotfn, "get_activeqid", parent)

            # Get dates and types and sort with respect to date
            dates = []; types = []
            for qid in qids:
                if parent == "results":
                    dates.append(ascot["q"+qid].get_date())
                    types.append("run")
                else:
                    dates.append(ascot[parent]["q"+qid].get_date())
                    types.append(ascot[parent]["q"+qid].get_type())

            sorted_datetypeqid = sorted(zip(dates, types, qids), reverse=True)
            sorted_datetypeqid = list(zip(*sorted_datetypeqid)) # Unzip
            dates = sorted_datetypeqid[0]
            types = sorted_datetypeqid[1]
            qids  = sorted_datetypeqid[2]

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


    def activate_item(self, item):
        """
        Set item (group) as active in the tree.
        """
        parent = self.tree.parent(item)

        # Move active tag from the previous active to the current one.
        for c in self.tree.get_children(parent):
            tags = list(self.tree.item(c, "tags"))

            if "active" in tags:
                tags.remove("active")
                self.tree.item(c, tags=tuple(tags))
            if c == item[0] or c == item:
                tags.append("active")
                self.tree.item(c, tags=tuple(tags))


    def remove_item(self, item, ascotfn):
        """
        Remove item (group) from tree.
        """
        group = self.tree.item(self.tree.selection(),"text")
        if group in INPUT_PARENTS + ["results"]:
            # A whole group was removed.
            for c in self.tree.get_children(item):
                self.tree.delete(c)

            if group == "results":
                self.tree.delete(item)
            else:
                self.tree.item(item, tags=("parent", "nodata"))
            return

        parentitem = self.tree.parent(item)
        parent     = self.tree.item(parentitem, "text")
        try:
            qids = tools.call_ascot5file(ascotfn, "get_qids", parent)
        except:
            # Removed group was the last one.
            self.tree.delete(item)
            if parent == "results":
                self.tree.delete(parentitem)
            else:
                self.tree.item(parentitem, tags=("parent", "nodata"))
            return

        # Remove group and set the next one as active
        self.tree.delete(item)
        activeqid = tools.call_ascot5file(
            ascotfn, "get_activeqid", parent)

        for group in self.tree.get_children(parentitem):
            if activeqid == self.tree.item(group, "text"):
                self.activate_item(group)
                break


    def highlightinputs(self, inputqids):
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


    def highlightoutputs(self, outputqids):
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


class GroupMenu(tk.Menu):
    """
    Popupmenu displayed when groups treeview is right-clicked.
    """

    def __init__(self, mainwindow, tree, *args, **kwargs):
        super().__init__(mainwindow, *args, tearoff=0, **kwargs)

        self.tree = tree


        self.add_command(label="Activate",
                         command=mainwindow.ascotfile_activate)
        self.add_command(label="Remove",
                         command=mainwindow.ascotfile_remove)
        self.add_command(label="Export",
                         command=mainwindow.ascotfile_export)
        self.add_command(label="Add dummy input",
                         command=mainwindow.ascotfile_adddummy)

        # Show menu when tree is right-clicked
        tree.bind("<Button-3>", self.showmenu)

        # Hide menu when focus is lost or user clicks anywhere else
        self.bind("<FocusOut>", self.hidemenu)
        mainwindow.bind("<Button-1>", self.hidemenu)


    def showmenu(self, e):
        """
        Show menu if a group was right-clicked and select the group.
        """
        item = self.tree.identify("item", e.x, e.y)
        itemname = self.tree.item(item, "text")
        if len(itemname) > 0:
            self.post(e.x_root, e.y_root)
            self.focus_set()

            # Disable "add dummy input" option if selected item is not parent
            self.entryconfigure(4, state="normal")
            if itemname not in INPUT_PARENTS + ["results"]:
                self.entryconfigure(4, state="disabled")

            # Disable "activate" option is selected item is a parent
            self.entryconfigure(0, state="normal")
            if itemname in INPUT_PARENTS + ["results"]:
                self.entryconfigure(0, state="disabled")

            # This will launch a event that will update the GUI.
            self.tree.selection_set(item)


    def hidemenu(self, e):
        """
        Hide menu.
        """
        self.unpost()


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


class AscotInitException(Exception):
    """Exception raised when Ascot object could not be initialized."""
    pass


class AscotpyInitException(Exception):
    """Exception raised when Ascotpy object could not be initialized."""
    pass


if __name__ == "__main__":
    gui = GUI("dummy.h5")
    gui.mainloop()
