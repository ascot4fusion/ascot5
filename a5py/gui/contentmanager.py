"""Defines ContentManager class for showing stuff at Canvas and Settings frames.
"""
import tkinter as tk
from tkinter import ttk

from a5py import AscotInitException

from .filecontents  import Info
from .inputcontents import Field, Preflight
from .resultcontents import Summary, StateScatter, StateHistogram, Orbit, \
    Poincare, LossSummary, Dists, Moments, WallLoad, Wall3D
from .runcontents import Trace
from .components import NestedNotebook

class ContentManager(NestedNotebook):
    """Manages the contents in SettingsFrame and in CanvasFrame.

    Settings frame is a notebook widget where changing the tab also changes
    what is shown on Canvas*. This is done so that we always pack_forget the
    current frame on Canvas and pack the new one (Note that the frames are
    kept as that way we preserve the settings user has given previously, which
    makes smoother user-experience when different frames are toggled).

    *For group contents the contents on Canvas also depends what is shown in
    the group treeview.

    This is fairly simple class as the actual contents are defined in Content*
    classes that are imported an used here. If those also have subcontents and
    implement a notebook to choose what is display, then those classes mimic the
    structure in this class. In theory, everything could be defined here but
    that would make this class huge...
    """

    def __init__(self, gui, settingsframe, canvasframe):
        """Generate all content widgets (which also generate their contents).
        """
        super().__init__(settingsframe)
        self.pack(fill="both", expand=True)
        self.gui = gui
        self.canvas = canvasframe

        def tabselected(tab):
            if self._sleep: return
            if tab == "Inputs":
                try:
                    self.gui.ascot.simulation_free()
                except AscotInitException:
                    pass
                init = {"wall" : False}
                for inp in ["bfield", "efield", "plasma", "neutral", "boozer",
                            "mhd", "asigma", "RF"]:
                    init[inp] = inp in self.gui.ascot.data
                #msg = self.gui.pleasehold("Ascotpy is being initialized..."
                self.gui.ascot.input_init(**init, switch=True)
            elif tab == "Results":
                try:
                    self.gui.ascot.simulation_free()
                except AscotInitException:
                    pass
                
                RF = 'RF' in self.gui.ascot.data
                self.gui.ascot.input_init(
                    run=True, bfield=True, efield=True, plasma=True,
                    neutral=True, boozer=True, mhd=True, switch=True,
                    RF=RF)
            elif tab == "Run":
                self.gui.ascot.input_free()
                self.gui.ascot.simulation_initinputs()

        args = (self, self.canvas, self.gui)

        self.add("File", tab=Info(*args))
        self.add("Inputs", tabselected=lambda : tabselected("Inputs"))
        self.traverse("Inputs").add("Plot (R,z)", tab=Field(*args))
        self.traverse("Inputs").add("Preflight", tab=Preflight(*args))
        self.add("Results", tabselected=lambda : tabselected("Results"))
        self.traverse("Results").add("Summary", tab=Summary(*args))
        self.traverse("Results").add("Ini/Endstate")
        self.traverse("Ini/Endstate").add("Scatter", tab=StateScatter(*args))
        self.traverse("Ini/Endstate").add("Histogram", tab=StateHistogram(*args))
        self.traverse("Results").add("Orbits")
        self.traverse("Orbits").add("Trajectory", tab=Orbit(*args))
        self.traverse("Orbits").add("Poincaré", tab=Poincare(*args))
        self.traverse("Results").add("Dists")
        self.traverse("Dists").add("Distribution", tab=Dists(*args))
        self.traverse("Dists").add("Moments", tab=Moments(*args))
        self.traverse("Results").add("Losses")
        self.traverse("Losses").add("Total", tab=LossSummary(*args))
        self.traverse("Losses").add("Wall loads", tab=WallLoad(*args))
        self.traverse("Losses").add("View 3D", tab=Wall3D(*args))
        self.add("Run", tabselected=lambda : tabselected("Run"))
        self.traverse("Run").add("Trace", tab=Trace(*args))

    def update_content(self):
        """Redraw contents on settings and canvas frames.
        """
        tab, name = self.currenttab()
        if name == "File":
            tab.selecttab()

    def restart(self):
        # Free any used resources
        try:
            self.gui.ascot.simulation_free()
        except:
            self.gui.ascot.input_free()

        # Has results?
        if self.gui.ascot.data.active is None:
            self.tab(2, state="disabled")
        else:
            self.tab(2, state="normal")

            run = self.gui.ascot.data.active
            tab0 = self.traverse("Results")

            # Has orbit data?
            if not hasattr(run, "_orbit"):
                tab0.tab(2, state="disabled")
            else:
                tab0.tab(2, state="normal")
                tab1 = self.traverse("Orbits")

                # Has Poincare data?
                if run.getorbit_poincareplanes() is None:
                    tab1.tab(1, state="disabled")
                else:
                    tab1.tab(1, state="normal")

            # Has distribution data?
            if not ( hasattr(run, "_dist5d") or hasattr(run, "_distrho5d") or
                     hasattr(run, "_dist6d") or hasattr(run, "_distrho6d") or
                     hasattr(run, "_distcom") ):
                tab0.tab(3, state="disabled")
            else:
                tab0.tab(3, state="normal")
                tab1 = self.traverse("Dists")

                # Has 5D distribution data to calculate moments?
                if not( hasattr(run, "_dist5d") or hasattr(run, "_distrho5d") ):
                    tab1.tab(1, state="disabled")
                else:
                    tab1.tab(1, state="normal")

            # Has 3D wall?
            tab1 = self.traverse("Losses")
            if run.wall.get_type() != "wall_3D":
                tab1.tab(1, state="disabled")
                tab1.tab(2, state="disabled")
            else:
                tab1.tab(1, state="normal")
                tab1.tab(2, state="normal")
