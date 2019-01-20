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

        self._inputframes = {}
        if hasattr(gui._ascot, "options"):
            self._panels["options"] = InputInfoFrame(gui, self, "options")
            self._panels["options"].place(0, 60, 440, 140)

        if hasattr(gui._ascot, "bfield"):
            self._panels["bfield"] = InputInfoFrame(gui, self, "bfield")
            self._panels["bfield"].place(0, 200, 440, 140)

        if hasattr(gui._ascot, "efield"):
            self._panels["efield"] = InputInfoFrame(gui, self, "efield")
            self._panels["efield"].place(0, 340, 440, 140)

        if hasattr(gui._ascot, "marker"):
            self._panels["marker"] = InputInfoFrame(gui, self, "marker")
            self._panels["marker"].place(0, 480, 440, 140)

        if hasattr(gui._ascot, "wall"):
            self._panels["wall"] = InputInfoFrame(gui, self, "wall")
            self._panels["wall"].place(460, 60, 440, 140)

        if hasattr(gui._ascot, "plasma"):
            self._panels["plasma"] = InputInfoFrame(gui, self, "plasma")
            self._panels["plasma"].place(460, 200, 440, 140)

        if hasattr(gui._ascot, "neutral"):
            self._panels["neutral"] = InputInfoFrame(gui, self, "neutral")
            self._panels["neutral"].place(460, 340, 440, 140)

        if hasattr(gui._ascot, "active"):
            self._panels["results"] = InputInfoFrame(gui, self)
            self._panels["results"].place(460, 480, 440, 140)

    def select_inputs(options=None, bfield=None, efield=None, plasma=None,
                      marker=None, neutral=None, wall=None)
        """
        Set given inputs as selected in their corresponding panels.
        """
        if options is not None:
            self._inputframes["options"].(options)
        if bfield is not None:
            self._inputframes["bfield"].(bfield)
        if efield is not None:
            self._inputframes["efield"].(efield)
        if plasma is not None:
            self._inputframes["plasma"].(plasma)
        if neutral is not None:
            self._inputframes["neutral"].(neutral)
        if marker is not None:
            self._inputframes["marker"].(neutral)
        if wall is not None:
            self._inputframes["wall"].(wall)


class InputInfoFrame(tkinter.Frame):

    """
    Make a panel for showing and adjusting input data.

    Panel contains menu for choosing different inputs within the parent this
    panel corresponds to. Active group is always top on the menu. The active
    group can be changed from this panel as well as description which is
    also shown. These changes are saved to HDF5 file.
    """
    def __init__(self, gui, indexframe, name)
        super().__init_(indexframe)
        self._indexframe = indexframe
        self._name = name

        self._inputs = tools.call_ascot5file(gui._h5fn, "get_qids", name)
        activeqid = tools.call_ascot5file(gui._h5fn, "get_activeqid", name)

        # Put active QID first.
        qids.remove(activeqid)
        qids = [activeqid] + qids
        groups = []
        for qid in qids:
            groups.append(gui._ascot[name]["q"+qid].get_type() + "-" + qid)

        self._inputselection = tkinter.StringVar(self)
        self._inputselection.set(groups[0])

        desc = gui._ascot[name]["q"+qids[0]].get_desc()
        date = gui._ascot[name]["q"+qids[0]].get_date()

        self._namelabel     = tkinter.Label( self, text=name)
        self._datelabel     = tkinter.Label( self, text="Created: " + date)
        self._activebutton  = tkinter.Button(self, text="Set active")
        self._savebutton    = tkinter.Button(self, text="Save description")
        self._viewbutton    = tkinter.Button(self, text="View")
        self._descbox       = tkinter.Text(self)
        self._selectionmenu = tkinter.OptionMenu(self, self._inputselection, 
                                                 *groups)

        self._descbox.delete("1.0", tkinter.END)
        self._descbox.insert("end", desc)

        self._inputselection.trace('w', self._change_input)

        self._activebutton.config(command=self._set_active)
        self._savebutton.config(command=self._set_desc)

        self._namelabel.place(    x=0,   y=10,   width=60,  height=20)
        self._selectionmenu.place(x=0,   y=30,   width=200, height=30)
        self._activebutton.place( x=200, y=32,   width=80,  height=26)
        self._savebutton.place(   x=280, y=32,   width=120, height=26)
        self._viewbutton.place(   x=400, y=32,   width=40,  height=26)
        self._descbox.place(      x=0,   y=60,   width=w,   height=h-80)
        self._datelabel.place(    x=0,   y=h-20, width=w/2,   height=20)

    def select_input():
        pass

    def change_input(*args):
        group = self._inputselection.get()
        desc = self._gui._ascot[self._name][group].get_desc()
        date = self._gui._ascot[self._name][group].get_date()
        self._descbox.delete("1.0", tkinter.END)
        self._descbox.insert("end", desc)
        datelabel.config(text="Created: " + date)

    def set_active():
        g = tkvar.get()
        tools.call_ascot5file(gui._h5fn, "set_active", g)
        gui.reload()

    def set_desc():
        g = tkvar.get()
        newdesc = descbox.get("1.0","end-1c")
        tools.call_ascot5file(gui._h5fn, "set_desc", g, newdesc)
        gui.reload()

    def view():
        g = tkvar.get()
        g = gui._ascot[name][g]
        inputtype = gui._ascot[name][g].get_type()
        if inputtype == "B_GS":
            gui.displayframe(BGSFrame(gui))
        if inputtype == "B_2DS":
            gui.displayframe(B2DSFrame)
        if inputtype == "B_3DS":
            gui.displayframe(B3DSFrame)
        if inputtype == "P_1D":
            gui.displayframe(Plasma1DFrame)
        if inputtype == "wall_2D":
            gui.displayframe(Wall2DFrame)
        if inputtype == "wall_3D":
            gui.displayframe(Wall3DFrame)

class RunInfoPanel(tkinter.Frame):
    """
        
    """

    def __init__(self, gui, indexframe)
        super().__init__(indexframe)
        self._indexframe = indexframe
        self._gui = gui

        # Obtain QIDs for all runs and put the active one first on the list.
        runqids   = tools.call_ascot5file(gui._h5fn, "get_qids", "results")
        activeqid = tools.call_ascot5file(gui._h5fn, "get_activeqid", "results")
        runqids.remove(activeqid)
        runqids = [activeqid] + runqids

        # Initialize run selection.
        desc = gui._ascot["q"+qids[0]].get_desc()
        date = gui._ascot["q"+qids[0]].get_date()
        self._selectedrun = tkinter.StringVar(self)
        self._selectedrun.trace('w', change_run)

        # Create twingets and set functionality.
        self._label           = tkinter.Label( self, text="results")
        self._datelabel       = tkinter.Label( self, text="Simulated: " + date)
        self._activebutton    = tkinter.Button(self, text="Set active")
        self._inputsbutton    = tkinter.Button(self, text="Show inputs")
        self._savebutton      = tkinter.Button(self, text="Save description")
        self._statebutton     = tkinter.Button(self, text="Ini/Endstate")
        self._dist5dbutton    = tkinter.Button(self, text="Dist 5D")
        self._dist6dbutton    = tkinter.Button(self, text="Dist 6D")
        self._distrho5dbutton = tkinter.Button(self, text="Dist rho5D")
        self._distrho6dbutton = tkinter.Button(self, text="Dist rho6D")
        self._descbox         = tkinter.Text(self)
        self._popupMenu       = tkinter.OptionMenu(self, self.selectedrun, 
                                                   *runqids)

        activebutton.config(   command=self._set_active)
        savebutton.config(     command=self._set_desc)
        inputsbutton.config(   command=self._set_inputs)
        statebutton.config(    command=lambda *args : self._view_plot("state"))
        dist5dbutton.config(   command=self._view_plot("dist5d"))
        dist6dbutton.config(   command=self._view_plot("dist6d"))
        distrho5dbutton.config(command=self._view_plot("distrho5d"))
        distrho6dbutton.config(command=self._view_plot("distrho6d"))
        orbitsbutton.config(   command=self._view_plot("orbits"))

        # Place twinkets.
        self._label.place(          x=0,   y=10,   width=60,  height=20)
        self._popupMenu.place(      x=0,   y=30,   width=200, height=30)
        self._activebutton.place(   x=200, y=32,   width=80,  height=26)
        self._inputsbutton.place(   x=300, y=32,   width=40,  height=26)
        self._savebutton.place(     x=280, y=32,   width=120, height=26)
        self._statebutton.place(    x=0,   y=52,   width=40,  height=26)
        self._dist5dbutton.place(   x=40,  y=52,   width=40,  height=26)
        self._dist6dbutton.place(   x=80,  y=52,   width=40,  height=26)
        self._distrho5dbutton.place(x=120, y=52,   width=40,  height=26)
        self._distrho6dbutton.place(x=160, y=52,   width=40,  height=26)
        self._orbitbutton.place(    x=200, y=52,   width=40,  height=26)
        self._descbox.place(        x=0,   y=80,   width=w,   height=h-100)
        self._datelabel.place(      x=0,   y=h-20, width=w/2, height=20)

        # Deactivate buttons that view output data which is not present.
        self._selectedrun("run-" + qids[0])

    def _activatebutton():
        pass

    def _deactivate_button():
        """
        Deactivate button.
        """
        if "" not in gui._ascot[]
        pass

    def select():
        pass

    def _change_run(*args):
        """
        Change run who's information is displayed.
        """
        run  = self._selectedrun.get()
        desc = self._gui._ascot[run].get_desc()
        date = self._gui._ascot[run].get_date()
        self._descbox.delete("1.0", tkinter.END)
        self._descbox.insert("end", desc)
        self._datelabel.config(text="Simulated: " + date)
        self._deactivate_buttons()

    def _set_active():
        """
        Set currently selected run as active in HDF5 and reload GUI.
        """
        run = self._selectedrun.get()
        tools.call_ascot5file(self._gui._h5fn, "set_active", run)
        gui.reload()

    def _set_desc():
        """
        Set description in HDF5 and reload GUI.
        """
        run  = self._selectedrun.get()
        desc = self._descbox.get("1.0", "end-1c")
        tools.call_ascot5file(self._gui._h5fn, "set_desc", run, desc)
        self._gui.reload()

    def _set_inputs():
        """
        Set inputs of the selected runs as selected in input panels.
        """
        run = self._gui._ascot["run-" + self._selectedrun.get()]
        self._indexframe.select_inputs(options=run.options, 
                                       bfield=run.bfield,
                                       efield=run.efield,
                                       plasma=run.plasma,
                                       neutral=run.neutral,
                                       marker=run.marker,
                                       wall=run.wall)

    def _view_plot(outputtype):
        """
        View plot screen of corresponding output data.
        """
        run = self._gui._ascot["run-" + self._selectedrun.get()]
        if outputtype == "state":
            self._gui.displayframe(StateFrame(self._gui, run.state))
        if outputtype == "dist5d":
            self._gui.displayframe(Dist5DFrame(self._gui, run.dist5D))
        if outputtype == "dist6d":
            self._gui.displayframe(Dist6DFrame(self._gui, run.dist6D))
        if outputtype == "distrho5d":
            self._gui.displayframe(DistRho5DFrame(self._gui, run.distrho5D)
        if outputtype == "distrho6d":
            self._gui.displayframe(DistRho6DFrame(self._gui, run.distrho6D))
        if outputtype == "orbits":
            self._gui.displayframe(OrbitsFrame(self._gui, run.orbits))