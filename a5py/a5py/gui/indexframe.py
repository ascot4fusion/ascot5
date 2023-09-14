"""
Contains definition of IndexFrame class.

Index frame is the frame that greets you when the GUI is opened. From there you
can navigate to input and output frames, choose which groups are active, modify
their descriptions or open a different HDF5 file.

File: indexframe.py
"""
import tkinter
import tkinter.ttk as ttk
import a5py.ascot5io.ascot5tools as tools

from .bfieldframe  import BfieldFrame
from .efieldframe  import EfieldFrame
from .plasmaframe  import PlasmaFrame
from .neutralframe import NeutralFrame
from .wallframe    import WallFrame
from .boozerframe  import BoozerFrame
from .mhdframe     import MhdFrame
from .optionsframe import OptionsFrame
from .markerframe  import MarkerFrame

from .stateframe  import StateFrame
from .orbitframe  import OrbitFrame
from .distframe   import DistFrame

from .sanityframe  import SanityFrame


# Some constants defining widget sizes (in px)
DROPDOWNMENUWIDTH = 25
PANELWIDTH  = 500
PANELHEIGHT = 50

MAXROWS = 5

class IndexFrame(tkinter.Frame):
    """
    An index frame where other frames can be accessed and HDF5 file modified.
    """

    def __init__(self, gui):
        """
        Initialize index frame.

        Index frame contains text panel showing the path to the current file and
        a button to browse and load a new file. For each input parent a small
        panel is displayed, and similar panel is displayed for results as well.
        See InputInfoPanel and RunInfoPanel for their contents.

        We don't worry about redrawing as everytime HDF5 file, it's contents, or
        Ascot5 folder is changed we just dump this frame and construct a new
        one.
        """
        super().__init__(gui._root)
        self._gui = gui

        # Create filename panel and fill panels to edges.
        fnpanel = GeneralInfoPanel(gui, self, PANELWIDTH, PANELHEIGHT)
        fnpanel.grid(row=0, column=1, sticky="NSWE")
        tkinter.Frame(self).grid(row=0,column=0,rowspan=MAXROWS+1,
                                 sticky="NSWE")
        tkinter.Frame(self).grid(row=0,column=3,rowspan=MAXROWS+1,
                                 sticky="NSWE")
        tkinter.Frame(self).grid(row=MAXROWS,column=1,columnspan=2,
                                 sticky="NSWE")

        # The columns featuring panels remain unchanged in size
        self.grid_columnconfigure(1, weight=0)
        self.grid_columnconfigure(2, weight=0)

        # The last column is free to resize
        self.grid_columnconfigure(3, weight=1)

        # All rows can grow in size:
        for i in range(MAXROWS+1):
            self.grid_rowconfigure(i, weight=1)

        # Set weight to 1 to obtain horizontally centered view
        self.grid_columnconfigure(0, weight=0)


        self._panels = {}
        ascot = gui.get_ascotobject()

        # Add input panels.
        rowcol = [1, 1]
        for inp in ["options", "bfield", "efield", "marker", "wall", "plasma",
                    "neutral", "boozer", "mhd"]:
            if hasattr(ascot, inp):
                self._panels[inp] = InputInfoPanel(gui, self, inp,
                                                   PANELWIDTH, PANELHEIGHT)
                self._panels[inp].grid(row=rowcol[0], column=rowcol[1],
                                       padx=0, pady=0, sticky="NSWE")

            rowcol[0] += 1
            if rowcol[0] > MAXROWS-1:
                rowcol[0] = 1
                rowcol[1] += 1

        # Add run panel.
        if hasattr(ascot, "active"):
            self._panels["results"] = RunInfoPanel(gui, self,
                                                   PANELWIDTH, PANELHEIGHT)
            self._panels["results"].grid(row=0, column=2,
                                         padx=0, pady=0, sticky="NSWE")


    def browsefile(self):
        """
        Browse a new HDF5 file.
        """
        self._gui.ask_openascot()


    def sanitycheck(self):
        """
        Open sanity check frame.
        """
        ascotpy = self._gui.get_ascotpy()

        if ascotpy is None \
           or "marker" not in self._panels \
           or "bfield" not in self._panels:
            return

        group = self._panels["bfield"]._inputselection.get()
        qid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
        group  = self._panels["marker"]._inputselection.get()
        marker = self._gui.get_ascotobject()["marker"][group]
        self._gui.displayframe(SanityFrame(self._gui, ascotpy, qid, marker))


    def analyze(self):
        """
        Open analysis frame.
        """
        #self._gui.ask_openascot()
        pass


    def viewinput(self, inputtype):
        """
        Open a new frane for viewing the selected input type.
        """
        ascotpy = self._gui.get_ascotpy()

        if ascotpy is None:
            return

        walldata = None
        if "wall" in self._panels:
            group    = self._panels["wall"]._inputselection.get()
            walldata = self._gui.get_ascotobject()["wall"][group]

        if inputtype in ["B_GS", "B_2DS", "B_3DS", "B_STS", "B_TC"]:
            group = self._panels["bfield"]._inputselection.get()
            qid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
            self._gui.displayframe(BfieldFrame(self._gui, ascotpy, qid,
                                               walldata))

        if inputtype in ["E_TC", "E_1DS"]:
            if "bfield" in self._panels:
                group = self._panels["bfield"]._inputselection.get()
                bqid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
            else:
                raise RuntimeError("bfield input is required.")

            group = self._panels["efield"]._inputselection.get()
            qid   = self._gui.get_ascotobject()["efield"][group].get_qid()
            self._gui.displayframe(EfieldFrame(self._gui, ascotpy, qid, bqid,
                                               walldata))

        if inputtype in ["plasma_1D", "plasma_1DS"]:
            if "bfield" in self._panels:
                group = self._panels["bfield"]._inputselection.get()
                bqid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
            else:
                raise RuntimeError("bfield input is required.")

            group = self._panels["plasma"]._inputselection.get()
            qid   = self._gui.get_ascotobject()["plasma"][group].get_qid()
            self._gui.displayframe(PlasmaFrame(self._gui, ascotpy, qid, bqid,
                                               walldata))

        if inputtype in ["N0_3D"]:
            if "bfield" in self._panels:
                group = self._panels["bfield"]._inputselection.get()
                bqid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
            else:
                raise RuntimeError("bfield input is required.")

            group = self._panels["neutral"]._inputselection.get()
            qid   = self._gui.get_ascotobject()["neutral"][group].get_qid()
            self._gui.displayframe(NeutralFrame(self._gui, ascotpy, qid, bqid,
                                                walldata))

        if inputtype in ["wall_2D", "wall_3D"]:
            self._gui.displayframe(WallFrame(self._gui, walldata))

        if inputtype in ["Boozer"]:
            if "bfield" in self._panels:
                group = self._panels["bfield"]._inputselection.get()
                bqid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
            else:
                raise RuntimeError("bfield input is required.")

            group = self._panels["boozer"]._inputselection.get()
            qid   = self._gui.get_ascotobject()["boozer"][group].get_qid()
            self._gui.displayframe(BoozerFrame(self._gui, ascotpy, qid, bqid,
                                               walldata))

        if inputtype in ["MHD_STAT", "MHD_NONSTAT"]:
            if "bfield" in self._panels:
                group = self._panels["bfield"]._inputselection.get()
                bqid   = self._gui.get_ascotobject()["bfield"][group].get_qid()
            else:
                raise RuntimeError("bfield input is required.")

            if "boozer" in self._panels:
                group   = self._panels["boozer"]._inputselection.get()
                boozqid = self._gui.get_ascotobject()["boozer"][group].get_qid()
            else:
                raise RuntimeError("boozer input is required.")

            group = self._panels["mhd"]._inputselection.get()
            qid   = self._gui.get_ascotobject()["mhd"][group].get_qid()
            self._gui.displayframe(MhdFrame(self._gui, ascotpy, qid,bqid,
                                            boozqid, walldata))

        if inputtype in ["opt"]:
            group    = self._panels["options"]._inputselection.get()
            options = self._gui.get_ascotobject()["options"][group]
            self._gui.displayframe(OptionsFrame(self._gui, options))

        if inputtype in ["prt", "gc", "fl"]:
            group  = self._panels["marker"]._inputselection.get()
            marker = self._gui.get_ascotobject()["marker"][group]
            self._gui.displayframe(MarkerFrame(self._gui, marker))


    def view_output(self, outputtype):
        """
        View plot screen of corresponding output data.
        """
        run = "run_" + self._panels["results"]._runselection.get()
        run = self._gui.get_ascotobject()[run]
        if outputtype == "state":
            if hasattr(run, "endstate"):
                self._gui.displayframe(StateFrame(self._gui, run.inistate,
                                                  run.endstate))
            else:
                self._gui.displayframe(StateFrame(self._gui, run.inistate))
        if outputtype == "dist5d":
            self._gui.displayframe(DistFrame(self._gui, run.dist5d))
        if outputtype == "dist6d":
            self._gui.displayframe(DistFrame(self._gui, run.dist6d))
        if outputtype == "distrho5d":
            self._gui.displayframe(DistFrame(self._gui, run.distrho5d))
        if outputtype == "distrho6d":
            self._gui.displayframe(DistFrame(self._gui, run.distrho6d))
        if outputtype == "orbit":
            self._gui.displayframe(OrbitFrame(self._gui, run.orbit, run.wall))


    def select_inputs(self, options=None, bfield=None, efield=None, plasma=None,
                      marker=None, neutral=None, wall=None, boozer=None, mhd=None):
        """
        Set given inputs as selected in their corresponding panels.
        """
        inputobjects = locals()
        for inputname in inputobjects.keys():
            if inputname == "self":
                continue
            if inputobjects[inputname] is not None:
                self._panels[inputname].select(inputobjects[inputname])


class GeneralInfoPanel(ttk.LabelFrame):
    """
    Panel which shows the general information of the HDF5 file.
    """

    def __init__(self, gui, indexframe, width, height):
        super().__init__(indexframe, width=width, height=height,
                         text="ASCOT5 file")

        topframe = tkinter.Frame(self)
        browse   = tkinter.Button(topframe, text="Browse...",
                                  command=indexframe.browsefile, anchor="e")
        browse.pack(side="left")

        textbox = tkinter.Text(self, height=1, width=1)

        botframe = tkinter.Frame(self)
        sanity   = tkinter.Button(botframe, text="Sanity checks",
                                  command=indexframe.sanitycheck, anchor="e")
        analyze  = tkinter.Button(botframe, text="Analyze results",
                                  command=indexframe.sanitycheck, anchor="e")
        sanity.pack(side="left")
        analyze.pack(side="left")

        textbox.configure(state="normal")
        textbox.delete("1.0", tkinter.END)
        textbox.insert("end", gui.get_ascotfilename())
        textbox.configure(state="disabled")

        topframe.pack(anchor="w")
        textbox.pack(anchor="w", fill="x")
        botframe.pack(anchor="w")


class InputInfoPanel(ttk.LabelFrame):
    """
    A panel for showing and adjusting input data.

    Panel contains menu for choosing different inputs within the parent this
    panel corresponds to. Active group is always top on the menu. The active
    group can be changed from this panel as well as description which is
    also shown. These changes are saved to HDF5 file.
    """

    def __init__(self, gui, indexframe, name, width, height):
        super().__init__(indexframe, text=name, width=width, height=height)
        self._gui = gui
        self._indexframe = indexframe
        self._name = name

        # Obtain QIDs for all runs and put the active one first on the list.
        qids = tools.call_ascot5file(gui.get_ascotfilename(), "get_qids", name)
        activeqid = tools.call_ascot5file(gui.get_ascotfilename(),
                                          "get_activeqid", name)
        qids.remove(activeqid)
        qids = [activeqid] + qids

        # Instead of qids, we store the whole group name
        self._inputs = []
        for qid in qids:
            key = gui.get_ascotobject()[name]["q"+qid].get_name()
            self._inputs.append(key.replace("-","_"))

        # Initialize run selection.
        self._inputselection = tkinter.StringVar(self)
        self._inputselection.set(self._inputs[0])
        self._inputselection.trace('w', self._change_selection)

        # Create widgets and set functionality.
        buttonframe   = tkinter.Frame(self)
        selectionmenu = ttk.Combobox(buttonframe, width=DROPDOWNMENUWIDTH,
                                     textvariable=self._inputselection)
        activebutton  = tkinter.Button(buttonframe, text="Set active",
                                       bg="sky blue")
        savebutton    = tkinter.Button(buttonframe, text="Save description",
                                       bg="sky blue")
        viewbutton    = tkinter.Button(buttonframe, text="View",
                                       bg="sky blue")
        descbox       = tkinter.Text(self, height=1, width=1)
        datelabel     = tkinter.Label(self)

        selectionmenu.bind("<<ComboboxSelected>>", self._change_selection)
        selectionmenu["values"] = self._inputs
        activebutton.config(command=self._set_active)
        savebutton.config(command=self._set_desc)
        viewbutton.config(command=self._view)

        # Place widgets.
        selectionmenu.pack(side="left")
        activebutton.pack(side="left")
        savebutton.pack(side="left")
        viewbutton.pack(side="left")

        buttonframe.pack(anchor="w")
        descbox.pack(anchor="w", fill="both")
        datelabel.pack(anchor="w")

        # Finalize and select active input.
        self._selectionmenu = selectionmenu
        self._viewbutton    = viewbutton
        self._descbox       = descbox
        self._datelabel     = datelabel
        self._change_selection()

        # Resize (comment to disable)
        self.bind("<Configure>", self.onresize)


    def onresize(self, event):
        """
        Adjust descbox to better fit the resized frame.
        """
        self.configure(width=event.width, height=event.height)
        if event.height < 120:
            self._descbox.config(height=1)
        elif event.height < 140:
            self._descbox.config(height=2)
        else:
            self._descbox.config(height=3)


    def select(self, group):
        self._inputselection.set(group)
        self._change_selection()


    def _change_selection(self, *args):
        group = self._inputselection.get()
        group = group.replace("-","_")
        desc = self._gui.get_ascotobject()[self._name][group].get_desc()
        date = self._gui.get_ascotobject()[self._name][group].get_date()
        self._descbox.delete("1.0", tkinter.END)
        self._descbox.insert("end", desc)
        self._datelabel.config(text="Created: " + date)


    def _set_active(self):
        group = self._inputselection.get()
        tools.call_ascot5file(self._gui.get_ascotfilename(), "set_active",
                              group)
        self._gui.reload()


    def _set_desc(self):
        group = self._inputselection.get()
        desc = self._descbox.get("1.0","end-1c")
        tools.call_ascot5file(self._gui.get_ascotfilename(), "set_desc", group,
                              desc)
        self._gui.reload()


    def _view(self):
        group = self._inputselection.get()
        group = self._gui.get_ascotobject()[self._name][group]
        inputtype = group.get_type()
        self._indexframe.viewinput(inputtype)


class RunInfoPanel(ttk.LabelFrame):
    """
    A panel for showing data from different runs.

    This works similarly as the InputInfoPanel. Additional button is shown which
    sets input panels to shown inputs for the currently displayed run. Few more
    buttons are included that display frame for viewing output data. If data is
    not present, corresnponding button is set inactive.
    """

    def __init__(self, gui, indexframe, width, height):
        """
        Initialize a new run frame.
        """
        super().__init__(indexframe, width=width, height=height, text="results")
        self._indexframe = indexframe
        self._gui = gui

        # Obtain a list QIDs for all runs and put the active one first.
        runqids   = tools.call_ascot5file(gui.get_ascotfilename(), "get_qids",
                                          "results")
        activeqid = tools.call_ascot5file(gui.get_ascotfilename(),
                                          "get_activeqid", "results")
        runqids.remove(activeqid)
        runqids = [activeqid] + runqids

        # String variable which keeps track of the currently selected run.
        self._runselection = tkinter.StringVar(self)
        self._runselection.set(runqids[0])
        self._runselection.trace('w', self._change_selection)

        ## Create widgets and place them. ##
        #
        # topbuttomframe [selection][set active][show inputs][save desc]
        # botbuttonframe [ini/end][d5][d6][dr5][dr6][orb]
        # textbox        [description]
        # label          [date]

        topbuttonframe  = tkinter.Frame(self)
        selectionmenu   = ttk.Combobox(topbuttonframe,
                                       width=DROPDOWNMENUWIDTH-5,
                                       textvariable=self._runselection)

        activebutton    = tkinter.Button(topbuttonframe,
                                         text="Set active",
                                         bg="sky blue")
        inputsbutton    = tkinter.Button(topbuttonframe,
                                         text="Show inputs",
                                         bg="sky blue")
        savebutton      = tkinter.Button(topbuttonframe,
                                         text="Save description",
                                         bg="sky blue")

        botbuttonframe = tkinter.Frame(self)
        self._statebutton     = tkinter.Button(botbuttonframe,
                                               text="Ini/Endstate")
        self._dist5dbutton    = tkinter.Button(botbuttonframe,
                                               text="Dist 5D")
        self._dist6dbutton    = tkinter.Button(botbuttonframe,
                                               text="Dist 6D")
        self._distrho5dbutton = tkinter.Button(botbuttonframe,
                                               text="Dist rho5D")
        self._distrho6dbutton = tkinter.Button(botbuttonframe,
                                               text="Dist rho6D")
        self._orbitbutton     = tkinter.Button(botbuttonframe,
                                               text="Orbit")
        self._descbox         = tkinter.Text(self, height=1, width=1)

        self._datelabel = tkinter.Label(self)

        selectionmenu.pack(side="left")
        activebutton.pack(side="left")
        savebutton.pack(side="left")
        inputsbutton.pack(side="left")
        self._statebutton.pack(side="left")
        self._dist5dbutton.pack(side="left")
        self._dist6dbutton.pack(side="left")
        self._distrho5dbutton.pack(side="left")
        self._distrho6dbutton.pack(side="left")
        self._orbitbutton.pack(side="left")

        topbuttonframe.pack(anchor="w")
        botbuttonframe.pack(anchor="w")
        self._descbox.pack(anchor="w", fill="x")
        self._datelabel.pack(anchor="w")

        # Set functionality
        selectionmenu.bind("<<ComboboxSelected>>", self._change_selection)
        selectionmenu["values"] = runqids
        activebutton.config(command=self._set_active)
        savebutton.config(  command=self._set_desc)
        inputsbutton.config(command=self._set_inputs)
        self._statebutton.config(
            command=lambda *args : self._indexframe.view_output("state"))
        self._dist5dbutton.config(
            command=lambda *args : self._indexframe.view_output("dist5d"))
        self._dist6dbutton.config(
            command=lambda *args : self._indexframe.view_output("dist6d"))
        self._distrho5dbutton.config(
            command=lambda *args : self._indexframe.view_output("distrho5d"))
        self._distrho6dbutton.config(
            command=lambda *args : self._indexframe.view_output("distrho6d"))
        self._orbitbutton.config(
            command=lambda *args : self._indexframe.view_output("orbit"))

        # Select active run.
        self._change_selection()

        # Resize (comment to disable)
        self.bind("<Configure>", self.onresize)


    def onresize(self, event):
        """
        Adjust descbox to better fit the resized frame.
        """
        self.configure(width=event.width, height=event.height)
        if event.height < 160:
            self._descbox.config(height=1)
        elif event.height < 180:
            self._descbox.config(height=2)
        else:
            self._descbox.config(height=3)


    def select(self, runqid):
        """
        Select run to be displaeyd.
        """
        self._runselection.set(runqid)
        self._change_selection()


    def _change_selection(self, *args):
        """
        Change run who's information is displayed.
        """
        run = self._gui.get_ascotobject()["q"+self._runselection.get()]
        self._descbox.delete("1.0", tkinter.END)
        self._descbox.insert("end", run.get_desc())
        self._datelabel.config(text="Simulated: " + run.get_date())

        # Enable/disable buttons according to whether corresponding output is
        # present in the selected run.
        self._statebutton.config(state="disable")
        if hasattr(run, "inistate"):
            self._statebutton.config(state="normal")

        self._dist5dbutton.config(state="disable")
        if hasattr(run, "dist5d"):
            self._dist5dbutton.config(state="normal")

        self._dist6dbutton.config(state="disable")
        if hasattr(run, "dist6d"):
            self._dist6dbutton.config(state="normal")

        self._distrho5dbutton.config(state="disable")
        if hasattr(run, "distrho5d"):
            self._distrho5dbutton.config(state="normal")

        self._distrho6dbutton.config(state="disable")
        if hasattr(run, "distrho6d"):
            self._distrho6dbutton.config(state="normal")

        self._orbitbutton.config(state="disable")
        if hasattr(run, "orbit"):
            self._orbitbutton.config(state="normal")


    def _set_active(self):
        """
        Set currently selected run as active in HDF5 and reload GUI.
        """
        run = "q"+self._runselection.get()
        tools.call_ascot5file(self._gui.get_ascotfilename(), "set_active", run)
        self._gui.reload()


    def _set_desc(self):
        """
        Set description in HDF5 and reload GUI.
        """
        run  = "q"+self._runselection.get()
        desc = self._descbox.get("1.0", "end-1c")
        tools.call_ascot5file(self._gui.get_ascotfilename(), "set_desc", run,
                              desc)
        self._gui.reload()


    def _set_inputs(self):
        """
        Set inputs of the selected runs as selected in input panels.
        """
        run = self._gui.get_ascotobject()["run_" + self._runselection.get()]
        self._indexframe.select_inputs(options=run.options.get_name(),
                                       bfield=run.bfield.get_name(),
                                       efield=run.efield.get_name(),
                                       plasma=run.plasma.get_name(),
                                       neutral=run.neutral.get_name(),
                                       marker=run.marker.get_name(),
                                       wall=run.wall.get_name(),
                                       boozer=run.boozer.get_name(),
                                       mhd=run.mhd.get_name())
