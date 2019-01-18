"""
Contains definition of IndexFrame class.

File: indexframe.py
"""
import tkinter
from tkinter.filedialog import askopenfilename
import a5py.ascot5io.ascot5tools as tools


class IndexFrame(tkinter.Frame):
    """
    An opening frame where other frames can be accessed and HDF5 file modified.

    File: indexframe.py
    """

    def __init__(self, gui):
        """
        Initialize index frame.

        Index frame contains text panel showing the path to the current file and
        a button to browse and load a new file. For each input parent a small
        panel is displayed, see function make_inputactivationframe().
        """
        super().__init__(gui._root)

        filenameboxlabel = tkinter.Label(self, text="ASCOT5 file:", anchor="w")
        filenamebox      = tkinter.Text(self)

        def setfn(fn):
            gui._h5fn = fn
            filenamebox.configure(state="normal")
            filenamebox.delete("1.0", tkinter.END)
            filenamebox.insert("end", gui._h5fn)
            filenamebox.configure(state="disabled")

        def browse():
            fn = askopenfilename( title="Select file",
                                  filetypes = [("HDF5 files","*.h5")] )
            gui._h5fn = fn
            gui.reload()

        filebrowsebutton = tkinter.Button(self, text="Browse...",
                                          command=browse)
        setfn(gui._h5fn)

        filenameboxlabel.place(x=10, y=20, width=90, height=20)
        filenamebox.place(x=10, y=40, width=360, height=20)
        filebrowsebutton.place(x=100, y=20, width=90, height=20)

        if hasattr(gui._ascot, "options"):
            self.make_inputactivationframe(0, 60, 440, 140, "options", gui)

        if hasattr(gui._ascot, "bfield"):
            self.make_inputactivationframe(0, 200, 440, 140, "bfield", gui)

        if hasattr(gui._ascot, "efield"):
            self.make_inputactivationframe(0, 340, 440, 140, "efield", gui)

        if hasattr(gui._ascot, "marker"):
            self.make_inputactivationframe(0, 480, 440, 140, "marker", gui)

        if hasattr(gui._ascot, "wall"):
            self.make_inputactivationframe(460, 60, 440, 140, "wall", gui)

        if hasattr(gui._ascot, "plasma"):
            self.make_inputactivationframe(460, 200, 440, 140, "plasma", gui)

        if hasattr(gui._ascot, "neutral"):
            self.make_inputactivationframe(460, 340, 440, 140, "neutral", gui)

        if hasattr(gui._ascot, "active"):
            pass


    def make_inputactivationframe(self, x, y, w, h, name, gui):
        """
        Make a panel for showing and adjusting input data.

        Panel contains menu for choosing different inputs within the parent this
        panel corresponds to. Active group is always top on the menu. The active
        group can be changed from this panel as well as description which is
        also shown. These changes are saved to HDF5 file.
        """
        frame = tkinter.Frame(self)
        frame.place(x=x, y=y, width=w, height=h)

        qids      = tools.call_ascot5file(gui._h5fn, "get_qids", name)
        activeqid = tools.call_ascot5file(gui._h5fn, "get_activeqid", name)

        # Put active QID first.
        qids.remove(activeqid)
        qids = [activeqid] + qids
        groups = []
        for qid in qids:
            groups.append(gui._ascot[name]["q"+qid].get_type() + "-" + qid)

        tkvar = tkinter.StringVar(frame)
        tkvar.set(groups[0])

        desc = gui._ascot[name]["q"+qids[0]].get_desc()
        date = gui._ascot[name]["q"+qids[0]].get_date()

        label        = tkinter.Label( frame, text=name)
        datelabel    = tkinter.Label( frame, text="Created: " + date)
        activebutton = tkinter.Button(frame, text="Set active")
        savebutton   = tkinter.Button(frame, text="Save description")
        viewbutton   = tkinter.Button(frame, text="View")
        descbox      = tkinter.Text(frame)
        popupMenu    = tkinter.OptionMenu(frame, tkvar, *groups)

        descbox.delete("1.0", tkinter.END)
        descbox.insert("end", desc)

        def change_input(*args):
            qid = tkvar.get()[-10:]
            newdesc = gui._ascot[name]["q"+qid].get_desc()
            newdate = gui._ascot[name]["q"+qid].get_date()
            descbox.delete("1.0", tkinter.END)
            descbox.insert("end", newdesc)
            datelabel.config(text="Created: " + newdate)

        def set_active():
            g = tkvar.get()
            tools.call_ascot5file(gui._h5fn, "set_active", g)
            gui.reload()

        def set_desc():
            g = tkvar.get()
            newdesc = descbox.get("1.0","end-1c")
            tools.call_ascot5file(gui._h5fn, "set_desc", g, newdesc)
            gui.reload()

        tkvar.trace('w', change_input)

        activebutton.config(command=set_active)
        savebutton.config(command=set_desc)

        label.place(       x=0,   y=10,   width=60,  height=20)
        popupMenu.place(   x=0,   y=30,   width=200, height=30)
        activebutton.place(x=200, y=32,   width=80,  height=26)
        savebutton.place(  x=280, y=32,   width=120, height=26)
        viewbutton.place(  x=400, y=32,   width=40,  height=26)
        descbox.place(     x=0,   y=60,   width=w,   height=h-80)
        datelabel.place(   x=0,   y=h-20, width=w/2,   height=20)
