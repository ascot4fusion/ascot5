"""Wrapper class for RunNode that adds methods to plot and access results.

The motivation for this class is following. Consider that user wants to plot
a figure that shows i) particle orbits in (R,z), ii) their final (R,z) position,
and iii) the wall tiles. It is not clear whether this plot should be implemented
as a method for State or Orbit class as data from both are needed. Furthermore,
neither has access to the wall data. The only logical place for this method
therefore is in the RunNode that has access to all relevant data.
"""
import numpy as np
import pyvista as pv
import unyt

from a5py.exceptions import AscotNoDataException

import a5py.routines.plotting as a5plt
import a5py.wall as wall
from a5py.ascot5io import State, Orbits

class RunMixin():
    """Class with methods to access and plot orbit and state data.

    This class assumes it is inherited by ResultsNode.
    """

    def _require(self, *args):
        """Check if required data is present and raise exception if not.

        This is a helper function to quickly check that the data is available.

        Parameters
        ----------
        *args : `str`
            Name(s) of the required data.

        Raises
        ------
        AscotNoDataException
            Raised if the required data is not present.
        """
        for arg in args:
            if not hasattr(self, arg):
                raise AscotNoDataException(
                    "Data for \"" +  arg + "\" is required but not present.")

    def getstate_list(self):
        """List quantities that can be evaluated with :meth:`getstate`.
        """
        self._require("_inistate")
        return State.listqnts()

    def getorbit_list(self):
        """List quantities that can be evaluated with :meth:`getorbit`.
        """
        self._require("_orbit")
        qnts = Orbits.listqnts()
        if self.options.read()["ORBITWRITE_MODE"] != 0:
            del qnts["pncrid"]
            del qnts["pncrdir"]
        return qnts

    def getstate(self, *qnt, mode="gc", state="ini", ids=None, endcond=None):
        """Evaluate a marker quantity based on its ini/endstate.

        Inistate is marker's phase-space position right at the start of
        the simulation and endstate is the position at the end of
        the simulation.

        This function not only returns the marker phase space coordinates but
        also other quantities that can be inferred from it and information that
        is stored along with coordinates. For a complete list of available
        quantities, see.

        ASCOT5 stores both particle and guiding center phase-space position in
        all simulations. To differentiate these, quantities with suffix "prt",
        e.g. "xprt", return particle quantities and without suffix the guiding
        center quantity is returned.

        Parameters
        ----------
        *qnt : str
            Names of the quantities.
        state : {"ini", "end"}, optional
            Is the quantity evaluated at the ini- or endstate.
        ids : array_like, optional
            Filter markers by their IDs.
        endcond : str or list [str], optional
            Filter markers by their end conditions.

            See for a list of all possible end conditions or to list end
            conditions that are currently present in the data.

            Markers may have multiple end conditions active simultaneously. If
            just the name of the end condition e.g. "POLMAX" is passed, then all
            markers that have (at least) the ``POLMAX`` end condition are
            returned.

            If the end condition is preceded by "NOT", e.g. "NOT POLMAX", then
            markers that don't have that end condition are returned.

            Passing multiple end conditions in a single string returns markers
            that have all listed end conditions active, e.g. "MAXPOL MAXTOR"
            returns markers that have both ``POLMAX`` and ``TORMAX`` active
            simultaneously.

            Passing end condition strings as separate list items acts as
            a logical OR, e.g. ["POLMAX", "TORMAX"] returns markers that have
            either ``POLMAX`` or ``TORMAX`` active.

        Returns
        -------
        val : array_like
            The evaluated quantity sorted by marker ID.

        If multiple quantities are queried, they are returned as a list in
            the order they were listed in ``*qnt``.

        Raises
        ------
        ValueError
            Raised when the queried quantity could not be interpreted.
        AscotNoDataException
            Raised when data required for the operation is not present.
        AscotInitException
            If evaluating quantity required interpolating an input that
            was not initialized.
        """
        self._require("_inistate")
        if endcond is not None: self._require("_endstate")
        if state not in ["ini", "end"]:
            raise ValueError("Unrecognized state: " + state)
        if state == "end": self._require("_endstate")

        # Get or evaluate the quantity
        data = getattr(self, "_" + state + "state").get(*qnt, mode=mode)

        # Parse by ids and endcond
        idx = np.ones(data[0].shape, dtype=bool)
        if endcond is not None:
            if not isinstance(endcond, list): endcond = [endcond]

            # Go through each unique end cond and mark that end cond valid or
            # not. This can then be used to make udix as boolean mask array.
            uecs, uidx = np.unique(self._endstate.get("endcond"),
                                   return_inverse=True)
            mask = np.zeros(uecs.shape, dtype=bool)
            for i, uec in enumerate(uecs):
                for ec in endcond:
                    accept  = State.endcond_check(uec, ec)
                    mask[i] = mask[i] or accept

            idx = mask[uidx]

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self._inistate.get("ids"), ids))

        for i in range(len(data)):
            data[i] = data[i][idx]
        if "mu" in qnt:
            data[qnt.index("mu")].convert_to_units("eV/T")
        return data if len(data) > 1 else data[0]

    def getorbit(self, *qnt, ids=None, pncrid=None, endcond=None):
        """Return orbit data.

        Returns marker phase space coordinates and derived quantities along
        the recorded orbit, if the orbit recording was enabled.

        Parameters
        ----------
        *qnt : str
            Names of the quantities.
        ids : array_like, optional
            Filter markers by their IDs.
        pncrid : array_like, optional
            Filter data points by the Poincaré plane they correspond to.
        endcond : str or list [str], optional
            Filter markers by their end conditions.

            See :meth:`getstate` for details on how this argument is parsed and
            for a list of end conditions present in the data.

        Returns
        -------
        val : array_like
            The queried quantity sorted first by marker ID and then by mileage.

            If multiple quantities are queried, they are returned as a list in
            the order they were listed in ``*qnt``.

        Raises
        ------
        ValueError
            Raised when the queried quantity could not be interpreted.
        AscotNoDataException
            Raised when data required for the operation is not present.
        AscotInitException
            If evaluating quantity required interpolating an input that
            was not initialized.
        """
        self._require("_orbit", "_inistate", "_endstate")
        data = self._orbit.get(self._inistate, self._endstate, *qnt)
        idarr = self._orbit.get(self._inistate, self._endstate, "ids")[0]
        idx = np.ones(data[0].shape, dtype=bool)
        if endcond is not None:
            eids = self.getstate("ids", endcond=endcond)
            idx = np.logical_and(idx, np.in1d(idarr, eids))

        if pncrid is not None:
            pncridarr = self.orbit.get(self._inistate, self._endstate,
                                       "pncrid")[0]
            idx = np.logical_and(idx, np.in1d(pncridarr, pncrid))

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(idarr, ids))

        for i in range(len(data)):
            data[i] = data[i][idx]
        if "mu" in qnt:
            data[qnt.index("mu")].convert_to_units("eV/T")
        return data if len(data) > 1 else data[0]

    def getstate_markersummary(self):
        """Return a summary of marker end conditions present in the data.

        Returns
        -------
        summary : [`str`]
            Array of strings where each string has a format
            "# of markers : end condition".

            If markers were aborted then also error messages are listed as
            separate strings with the same format.

        Raises
        ------
        AscotNoDataException
            Raised when data required for the operation is not present.
        """
        self._require("_endstate")
        econd, emsg, emod, eline = self._endstate.get(
            "endcond", "errormsg", "errormod", "errorline")
        errors = np.unique(np.array([emsg, eline, emod]), axis=1).T

        ec, counts = np.unique(econd, return_counts=True)
        msg = ["End conditions:"]
        for i, e in enumerate(ec):
            econd = State.endcond_tostring(e)
            if econd == "none": econd = "aborted"
            msg.append(str(counts[i]) + " : " + econd)

        # It would be better to import these via libascot, but then this
        # function would require libascot and we don't want that.
        modules = [
            "mccc_wiener.c", "mccc_push.c", "mccc.c", "step_fo_vpa.c",
            "step_gc_cashkarp.c", "step_gc_rk4", "N0_3D.c", "N0_ST.c",
            "B_2DS.c", "B_2DS.c", "B_STS.c", "B_GS.c", "plasma_1D.c",
            "plasma_1DS.c", "plasma.c", "E_field.c", "neutral.c",
            "E_1DS.c", "B_field.c", "particle.c", "boozer.c", "mhd.c"
        ]
        messages = [
            "Input evaluation failed", "Unknown input type",
            "Unphysical quantity when evaluating input",
            "Unphysical marker quantity", "Time step too small/zero/NaN",
            "Wiener array is full or inconsistent",
            "Unphysical result when integrating marker coordinates"
        ]

        for e in errors:
            if np.sum(e) > 0:
                msg.append("Error: " + messages[e[0]-1] \
                           + ". At line " + str(e[1])  \
                           + " in " + modules[e[2]-1])

        return msg

    def getstate_losssummary(self):
        """Return a summary of lost markers.
        """
        self._require("endstate")

        wmrks = self.endstate.get("weight")
        wloss = self.endstate.get("weight", endcond="wall")
        emrks = self.endstate.get("energy")
        eloss = self.endstate.get("energy", endcond="wall")

        markers_lost   = wloss.size
        markers_frac   = wloss.size / wmrks.size
        particles_lost = np.sum(wloss)
        particles_frac = np.sum(wloss) / np.sum(wmrks)
        power_lost     = np.sum(wloss * eloss)
        power_frac     = np.sum(wloss * eloss) / np.sum(wmrks * emrks)

        msg = []
        msg += ["Markers lost: " + str(markers_lost) + " ("
                + str(np.around(markers_frac*100, decimals=1)) + "% of total)"]
        msg += ["Particles lost: " + str(particles_lost) + " ("
                + str(np.around(particles_frac*100, decimals=1)) + "% of total)"]
        msg += ["Energy lost: " + str(power_lost) + " ("
                + str(np.around(power_frac*100, decimals=1)) + "% of total)"]
        return msg

    def getstate_pointcloud(self, endcond=None):
        """Return marker endstate (x,y,z) coordinates in single array.

        Parameters
        ----------
            endcond : str, optional <br>
                Only return markers that have given end condition.
        """
        self._require("_endstate")
        return np.array([self._endstate.get("x", endcond=endcond),
                         self._endstate.get("y", endcond=endcond),
                         self._endstate.get("z", endcond=endcond)]).T

    def getorbit_poincareplanes(self):
        """Return a list of Poincaré planes that were present in the simulation.

        The list consists of "pol", "tor", and "rad" depending on the type of
        the plane.
        """
        opt  = self.options.read()
        npol = opt["ORBITWRITE_POLOIDALANGLES"]
        ntor = opt["ORBITWRITE_TOROIDALANGLES"]
        nrad = opt["ORBITWRITE_RADIALDISTANCES"]
        npol = 0 if npol[0] < 0 else npol.size
        ntor = 0 if ntor[0] < 0 else ntor.size
        nrad = 0 if nrad[0] < 0 else nrad.size

        plane = [None] * (npol + ntor + nrad)
        plane[:npol]                    = ["pol"]
        plane[npol:npol+ntor]           = ["tor"]
        plane[npol+ntor:npol+ntor+nrad] = ["rad"]

        return plane

    def getwall_figuresofmerit(self):
        """
        """
        self._require("endstate")
        ids    = self.endstate.get("walltile", endcond="wall")
        energy = self.endstate.get("energy", endcond="wall")
        weight = self.endstate.get("weight", endcond="wall")
        area   = self.wall.area()

        wetted_total, energy_peak = wall.figuresofmerit(ids, energy, weight, area)
        unit = "J"
        if energy_peak > 1e9:
            energy_peak /= 1e9
            unit = "GJ"
        elif energy_peak > 1e6:
            energy_peak /= 1e6
            unit = "MJ"
        elif energy_peak > 1e3:
            energy_peak /= 1e3
            unit = "KJ"

        msg = []
        msg += ["Total wetted area: " + str(np.around(wetted_total, decimals=2)) + r" $m^2$"]
        msg += ["Peak load: " + str(np.around(energy_peak, decimals=2)) + " " + str(unit) + r"$/m^2$"]
        return msg

    def getwall_loads(self):
        """Get 3D wall loads and associated quantities.

        This method does not return loads on all wall elements (as usually most
        of them receive no loads) but only those that are affected and their
        IDs.

        Returns:
            ids
            edepo
            eload
            pdepo
            pload
            mdepo
            iangle
        """
        self._require("endstate")
        ids    = self.endstate.get("walltile", endcond="wall")
        energy = self.endstate.get("energy", endcond="wall")
        weight = self.endstate.get("weight", endcond="wall")
        area   = self.wall.area()
        return wall.loads(ids, energy, weight, area)


    def getwall_3dmesh(self):
        """Return 3D mesh representation of 3D wall and associated loads.

        Returns
        -------
        wallmesh : Polydata
            Mesh representing the wall.

            The mesh cell data has fields:

            - "pload" particle load in units of prt/m^2 or prt/m^2s,
            - "eload" power/energy load in units of W/m^2 or J/m^2
        """
        wallmesh = pv.PolyData( *self.wall.noderepresentation() )
        ids, _, eload, _, pload, _, iangle = self.getwall_loads()
        ids = ids - 1 # Convert IDs to indices
        ntriangle = wallmesh.n_faces
        wallmesh.cell_data["pload"]       = np.zeros((ntriangle,)) + np.nan
        wallmesh.cell_data["pload"][ids]  = pload
        wallmesh.cell_data["eload"]       = np.zeros((ntriangle,)) + np.nan
        wallmesh.cell_data["eload"][ids]  = pload

        #- "iangle" angle of incidence in units of deg
        #wallmesh.cell_data["iangle"]      = np.zeros((ntriangle,)) + np.nan
        #wallmesh.cell_data["iangle"][ids] = iangle
        return wallmesh

    def plotstate_scatter(self, x, y, z=None, c=None, endcond=None, ids=None,
                          cint=9, cmap=None, axesequal=False, axes=None,
                          cax=None):
        """Make a scatter plot of marker state coordinates.

        The plot is either 2D+1D or 3D+1D, where the extra coordinate is color,
        depending on the number of queried coordinates. The color scale is
        discrete, not continuous.

        The quantity ("qnt") must have one of the following formats:

        - "ini qnt" plots inistate of quantity qnt.
        - "end qnt" plots endstate of quantity qnt.
        - "diff qnt" plots the difference between ini and endstate.
        - "reldiff qnt" plots the relative difference (x1/x0 - 1).
        - "log ini/end/diff/reldiff qnt" plots the logarithmic value.

        Parameters
        ----------
        x : str
            Name of the quantity on x-axis.
        y : str
            Name of the quantity on y-axis.
        z : str, optional
            Name of the quantity on z-axis.

            If not None, the plot will be in 3D.
        c : str, optional
            Name of the quantity shown with color scale or name of a color
            to plot all markers with same color.
        endcond : str, array_like, optional
            Endcond of those markers which are plotted.
        ids : array_like
            IDs of the markers to be plotted.
        cint : int or array_like optional
            Number of colors to be used or bin edges.

            All markers that have ``cint[i] < c <= cint[i+1]`` are plotted with
            same color. If ``cint`` array is not given explicitly,
            cint = np.linspace(min(c), max(c), cint).
        cmap : str, optional
            Name of the colormap where colors are picked.
        axesequal : bool, optional
            Flag whether to set aspect ratio for x and y (and z) axes equal.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.

        Raises
        ------
        AscotNoDataException
            Raised when data required for the opreation is not present.
        ValueError
            If argument could not be parsed.
        """
        def parsearg(arg):
            """Parse arguments of from "log ini qnt" to (qnt, value, islog).
            """
            arg = arg.lower()
            log = "linear"
            if "log" in arg:
                arg = arg.replace("log", "")
                log = "log"

            if "ini" in arg:
                arg = arg.replace("ini", "").strip()
                val = self.getstate(arg, state="ini", endcond=endcond, ids=ids)
                arg = "Initial " + arg
            elif "end" in arg:
                arg = arg.replace("end", "").strip()
                val = self.getstate(arg, state="end", endcond=endcond, ids=ids)
                arg = "Final " + arg
            elif "reldiff" in arg:
                arg = arg.replace("reldiff", "").strip()
                val1 = self.getstate(arg, state="ini", endcond=endcond, ids=ids)
                val2 = self.getstate(arg, state="end", endcond=endcond, ids=ids)
                val = (val2 -val1) / val1
                arg = r"$\Delta x/x_0$ " + arg
            elif "diff" in arg:
                arg = arg.replace("diff", "").strip()
                val1 = self.getstate(arg, state="ini", endcond=endcond, ids=ids)
                val2 = self.getstate(arg, state="end", endcond=endcond, ids=ids)
                val = val2 - val1
                arg = r"$\Delta$ " + arg
            else:
                raise ValueError(
                    "Unclear if a quantity is evaluated from ini or endstate.\n"
                    + "Use \"ini/end/diff/reldiff %s\"" % (arg)
                )
            if any(val < 0) and log == "log":
                log = "symlog"

            # Make sure val is an unyt array
            try:
                val.units
            except AttributeError:
                val *= unyt.dimensionless
            return (arg, val, log)

        x, xc, xlog = parsearg(x)
        y, yc, ylog = parsearg(y)
        x = x + " [" + str(xc.units) + "]"
        y = y + " [" + str(yc.units) + "]"

        cc = None; clog = False;
        if c is not None:
            if len(c.split()) == 1:
                # c is a color string, not quantity
                cc = c
                c = None
            else:
                c, cc, clog = parsearg(c)
                c = c + " [" + str(cc.units) + "]"

        if z is None:
            a5plt.scatter2d(xc, yc, c=cc, xlog=xlog, ylog=ylog, clog=clog,
                            xlabel=x, ylabel=y, clabel=c, cint=cint, cmap=cmap,
                            axesequal=axesequal, axes=axes, cax=cax)
        else:
            z, zc, zlog = parsearg(z)
            z = z + " [" + str(zc.units) + "]"
            a5plt.scatter3d(xc, yc, zc, c=cc, xlog=xlog, ylog=ylog,
                            zlog=zlog, clog=clog, xlabel=x, ylabel=y, zlabel=z,
                            clabel=c, cint=cint, cmap=cmap,
                            axesequal=axesequal, axes=axes, cax=cax)

    def plotstate_histogram(
            self, x, y=None, xbins=None, ybins=None, endcond=None, weight=False,
            iniend=["i", "i"], log=[False, False], logscale=False,
            cmap="viridis", axesequal=False, axes=None, cax=None):
        """Make a histogram plot of marker state coordinates.

        The histogram is either 1D or 2D depending on if the y coordinate is
        provided. In the 1D histogram the markers with different endstate are
        separated by color if the endcond argument is None.

        Parameters
        ----------
        x : str
            Name of the quantity on x-axis.
        y : str, optional
            Name of the quantity on y-axis. Makes the histogram 2D.
        cmap : str, optional
            Name of the colormap used in the 2D histogram.
        endcond : str, array_like, optional
            Endcond of those  markers which are plotted. Separated by color
            in 1D plot otherwise.
        iniend : str, array_like, optional
            Flag whether a corresponding [x, y] coordinate is taken from
            the ini ("i") or endstate ("e").
        log : bool, array_like, optional
            Flag whether the corresponding axis [x, y] is made logarithmic.
            If that is the case, absolute value is taken before the quantity
            is passed to log10.
        axesequal : bool, optional
            Flag whether to set aspect ratio for x and y axes equal in 2D.
        axes : Axes, optional
            The Axes object to draw on. If None, a new figure is displayed.
        cax : Axes, optional
            The Axes object for the color data in 2D, otherwise taken from
            axes.

        Raises
        ------
        AscotNoDataException
            Raised when data required for the opreation is not present.
        """
        if "i" in iniend: self._require("_inistate")
        if "e" in iniend: self._require("_endstate")
        if endcond is not None: self._require("_endstate")

        # Decipher whether quantity is evaluated from ini or endstate
        state = [None] * 2
        for i in range(2):
            if iniend[i] == "i":
                state[i] = self.inistate
            elif iniend[i] == "e":
                state[i] = self.endstate
            else:
                raise ValueError("Only 'i' and 'e' are accepted for iniend.")

        if y is None:
            # 1D plot

            xc = []
            weights = []

            ecs, count = self._endstate.listendconds()
            for ec in ecs:
                if endcond is not None and ec not in [endcond]:
                    pass

                xc.append(state[0].get(x, endcond=ec))
                weights.append(state[0].get("weight", endcond=ec))

            idx = np.argsort([len(i) for i in xc])[::-1]
            xc  = [xc[i] for i in idx]
            ecs = [ecs[i] for i in idx]
            weights = [weights[i] for i in idx]
            if not weight:
                weights = None

            xlabel = x + " [" + str(xc[0].units) + "]"
            a5plt.hist1d(x=xc, xbins=xbins, weights=weights,
                         log=[log[0], logscale],
                         xlabel=xlabel, axes=axes, legend=ecs)

        else:
            # 2D plot
            xc = state[0].get(x, endcond=endcond)
            yc = state[1].get(y, endcond=endcond)
            weights = state[0].get("weight", endcond=endcond)
            if not weight:
                weights = None

            xlabel = x + " [" + str(xc.units) + "]"
            ylabel = y + " [" + str(yc.units) + "]"
            a5plt.hist2d(xc, yc, xbins=xbins, ybins=ybins, weights=weights,
                         log=[log[0], log[1], logscale], xlabel=xlabel,
                         ylabel=ylabel, axesequal=axesequal, axes=axes, cax=cax)

    def plotorbit_trajectory(self, x, y, z=None, c=None, endcond=None, ids=None,
                             cmap=None, axesequal=False, axes=None, cax=None):
        """Plot orbit trajectories in arbitrary coordinates.

        The plot is either 2D+1D or 3D+1D, where the extra coordinate is color,
        depending on the number of queried coordinates. The color scale is
        discrete, not continuous.

        The quantity ("qnt") must have one of the following formats:

        - "diff qnt" plots the difference between inistate and current value.
        - "reldiff qnt" plots the relative difference (x1/x0 - 1).
        - "log <diff/reldiff> qnt" plots the logarithmic value.

        Parameters
        ----------
        x : str
            Name of the quantity on x-axis.
        y : str
            Name of the quantity on y-axis.
        z : str, optional
            Name of the quantity on z-axis.

            If not None, the plot will be in 3D.
        c : str, optional
            The color used to plot the markers or name of the quantity shown
            with color.

            If None, markers are plotted with different colors.
        endcond : str, array_like, optional
            Endcond of those  markers which are plotted.
        ids : array_like
            IDs of the markers to be plotted.
        cmap : str, optional
            Colormap.
        axesequal : bool, optional
            Flag whether to set aspect ratio for x and y (and z) axes equal.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        """
        def parsearg(arg):
            """Parse arguments of from string to (qnt, label, islog).
            """
            arg = arg.lower()
            log = "linear"
            label = arg
            if "log" in arg:
                arg = arg.replace("log", "")
                log = "log"
            if "reldiff" in arg:
                arg = arg.replace("reldiff", "")
                label = r"$\Delta x/x_0$ " + arg
            elif "diff" in arg:
                arg = arg.replace("diff", "")
                label = r"$\Delta$ " + arg

            arg = arg.strip()
            return (arg, label, log)

        cc = c; clabel = None; clog = "linear"# Default values passed to plotter
        x, xlabel, xlog = parsearg(x)
        y, ylabel, ylog = parsearg(y)
        if c is not None: c, clabel, clog = parsearg(c)
        if z is not None: z, zlabel, zlog = parsearg(z)

        if c is not None and not c in self.getorbit_list():
            # c is a color string, not quantity
            c = None

        # Orbit evaluation can be expensive so we try to get all coordinates
        # with a single call
        if z is None and c is None:
            idarr, xc, yc = self.getorbit(
                "ids", x, y, endcond=endcond, ids=ids)
        elif z is not None and c is None:
            idarr, xc, yc, zc = self.getorbit(
                "ids", x, y, z, endcond=endcond, ids=ids)
        elif z is None and c is not None:
            idarr, xc, yc, cc = self.getorbit(
                "ids", x, y, c, endcond=endcond, ids=ids)
        elif z is not None and c is not None:
            idarr, xc, yc, zc, cc = self.getorbit(
                "ids", x, y, z, c, endcond=endcond, ids=ids)

        # Find indices to map values from inistate to orbit array
        _, idarr = np.unique(idarr, return_inverse=True)
        idx = np.where(idarr[:-1] != idarr[1:])[0] + 1
        def parsevals(val, log, label, qnt):
            """Compute values and split them by orbit
            """
            if "x/x_0" in label:
                val0 = self.getstate(qnt, state="ini", endcond=endcond, ids=ids)
                val = ( val / val0[idarr] ) - 1
            elif "Delta" in label:
                val0 = self.getstate(qnt, state="ini", endcond=endcond, ids=ids)
                val = val - val0[idarr]

            if any(val < 0) and log == "log":
                log = "symlog"

            # Make sure val is an unyt array
            try:
                val.units
            except AttributeError:
                val *= unyt.dimensionless
            vmin = np.amin(val)
            vmax = np.amax(val)
            val = np.split(val, idx)
            return (val, log, vmin, vmax)

        xc, xlog, xmin, xmax = parsevals(xc, xlog, xlabel, x)
        yc, ylog, ymin, ymax = parsevals(yc, ylog, ylabel, y)
        xlabel = xlabel + " [" + str(xc[0].units) + "]"
        ylabel = ylabel + " [" + str(yc[0].units) + "]"
        bbox = [xmin, xmax, ymin, ymax, None, None]
        if z is not None:
            zc, zlog, zmin, zmax = parsevals(zc, zlog, zlabel, z)
            zlabel + " [" + str(zc[0].units) + "]"
            bbox = [xmin, xmax, ymin, ymax, zmin, zmax, None, None]
        if c is not None:
            cc, clog, bbox[-2], bbox[-1] = parsevals(cc, clog, clabel, c)
            clabel + " [" + str(cc[0].units) + "]"

        if z is None:
            a5plt.line2d(xc, yc, c=cc, xlog=xlog, ylog=ylog, clog=clog,
                         xlabel=xlabel, ylabel=ylabel, clabel=clabel, bbox=bbox,
                         cmap=cmap, axesequal=axesequal, axes=axes, cax=cax)
        else:
            a5plt.line3d(xc, yc, zc, c=cc, xlog=xlog, ylog=ylog, zlog=zlog,
                         clog=clog, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                         clabel=clabel, bbox=bbox, cmap=cmap,
                         axesequal=axesequal, axes=axes, cax=cax)

    def plotorbit_poincare(self, plane, conlen=True, axes=None, cax=None):
        """Create a Poincaré plot where the color separates different markers
        or shows the connection length.

        Parameters
        ----------
        plane : int
            Index number of the plane to be plotted in the list given by
            getorbit_poincareplanes.
        conlen : bool, optional
            If true, trajectories of lost markers are colored (in blue
            shades) where the color shows the connection length at that
            position. Confined (or all if conlen=False) markers are shown
            with shades of red where color separates subsequent
            trajectories.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.

        Raises
        ------
        AscotNoDataException
            Raised when data required for the opreation is not present.
        """
        if not self.has_orbit:
            raise LackingDataError("orbit")
        try:
            self.orbit["pncrid"]
        except:
            raise LackingDataError("orbit(poincare)")
        if not self.has_endstate:
            raise LackingDataError("endstate")

        # Set quantities corresponding to the planetype
        planes = self.getorbit_poincareplanes()
        if len(planes) <= plane:
            raise ValueError("Unknown plane.")

        ptype = planes[plane]
        if ptype == "pol":
            x = "r"
            y = "z"
            xlabel = "R [m]"
            ylabel = "z [m]"
            xlim = None # TBD
            ylim = None # TBD
            axesequal = True
        elif ptype == "tor":
            x = "rho"
            y = "phimod"
            xlabel = "Normalized poloidal flux"
            ylabel = "Toroidal angle [deg]"
            xlim = [0, 1.1]
            ylim = [0, 360]
            axesequal = False
        elif ptype == "rad":
            x = "thetamod"
            y = "phimod"
            xlabel = "Poloidal angle [deg]"
            ylabel = "Toroidal angle [deg]"
            xlim = [0, 360]
            ylim = [0, 360]
            axesequal = True

        ids = self.orbit.get("id", pncrid=plane)
        x   = self.orbit.get(x,    pncrid=plane)
        y   = self.orbit.get(y,    pncrid=plane)

        if ptype == "pol":
            # Determine (R,z) limits by making sure the data fits nicely and
            # then rounding to nearest decimal.
            xlim = [np.floor(np.amin(x)*9) / 10, np.ceil(np.amax(x)*11) / 10]
            ylim = [np.floor(np.amin(y)*11) / 10, np.ceil(np.amax(y)*11) / 10]

        if conlen:
            # Calculate connection length using endstate (only markers that pass
            # this plane are considered)
            total = self.endstate.get("mileage", ids=np.unique(ids))

            # Find how many data points each marker has in orbit data arrays.
            # idx are indices of last data point for given id.
            idx = np.argwhere(ids[1:] - ids[:-1]).ravel()
            idx = np.append(idx, ids.size-1)
            # Getting the number of data points now is just a simple operation:
            idx[1:] -= idx[:-1]
            idx[0] += 1

            # Now we can repeat the total mileage by the number of data points
            # each marker has, so we can calculate the connection length at each
            # point as: conlen = total - mileage
            total  = np.repeat(total, idx)
            conlen = total - self.orbit.get("mileage", pncrid=plane)

            # Now set confined markers as having negative connection length
            lost1 = self.endstate.get("id", endcond="rhomax")
            lost2 = self.endstate.get("id", endcond="wall")

            idx = ~np.logical_or(np.in1d(ids, lost1), np.in1d(ids, lost2))
            conlen[idx] *= -1
            clabel = "Connection length [" + str(conlen.units) + "]"
        else:
            conlen = None
            clabel = None

        a5plt.poincare(x, y, ids, conlen=conlen, xlim=xlim, ylim=ylim,
                       xlabel=xlabel, ylabel=ylabel, clabel=clabel,
                       axesequal=axesequal, axes=axes, cax=cax)


    def plotstate_summary(self, axes_inirho=None, axes_endrho=None,
                          axes_mileage=None, axes_energy=None,
                          axes_rz=None, axes_rhophi=None):
        """Plot several graphs that summarize the simulation.

        Following graphs are plotted:
          - inirho: Initial rho histogram with colors marking the endcond.
          - endrho: Initial rho histogram with colors marking the endcond.
          - mileage: Final mileage histogram with colors marking the endcond.
          - energy: Final energy histogram with colors marking the endcond.
          - Rz: Final R-z scatterplot.
          - rhophi: Final rho-phi scatterplot.

        Parameters
        ----------
            axes_inirho: {Axes, bool}, optional <br>
                Axes where inirho is plotted else a new figure is created. If
                False, then this plot is omitted.
            axes_endrho: {Axes, bool}, optional <br>
                Axes where inirho is plotted else a new figure is created. If
                False, then this plot is omitted.
            axes_mileage: {Axes, bool}, optional <br>
                Axes where inirho is plotted else a new figure is created. If
                False, then this plot is omitted.
            axes_energy: {Axes, bool}, optional <br>
                Axes where inirho is plotted else a new figure is created. If
                False, then this plot is omitted.
            axes_rz: {Axes, bool}, optional <br>
                Axes where inirho is plotted else a new figure is created. If
                False, then this plot is omitted.
            axes_rhophi: {Axes, bool}, optional <br>
                Axes where inirho is plotted else a new figure is created. If
                False, then this plot is omitted.

        Raises
        ------
        AscotNoDataException
            Raised when data required for the opreation is not present.
        """
        self._require("_inistate", "_endstate")
        # Initial rho histogram with colors marking endcond
        fig = None
        if axes_inirho is None:
            fig = plt.figure()
            axes_inirho = fig.add_subplot(1,1,1)
        if axes_inirho != False:
            axes_inirho.set_xlim(0,1.1)
            axes_inirho.set_title("Initial radial position")
            self.plotstate_histogram(
                "rho", xbins=np.linspace(0,1.1,55), weight=True,
                iniend=["i", "i"], axes=axes_inirho)
            if fig is not None: plt.show()

        # Final rho histogram with colors marking endcond
        fig = None
        if axes_endrho is None:
            fig = plt.figure()
            axes_endrho = fig.add_subplot(1,1,1)
        if axes_endrho != False:
            axes_endrho.set_xlim([0,1.1])
            axes_endrho.set_title("Final radial position")
            self.plotstate_histogram(
                "rho", xbins=np.linspace(0,1.1,55), weight=True,
                iniend=["e", "i"], axes=axes_endrho)
            if fig is not None: plt.show()

        # Mileage histogram with colors marking endcond
        fig = None
        if axes_mileage is None:
            fig = plt.figure()
            axes_mileage = fig.add_subplot(1,1,1)
        if axes_mileage != False:
            axes_mileage.set_title("Final mileage")
            self.plotstate_histogram(
                "mileage", xbins=55, weight=True,
                iniend=["e", "i"], log=[True, False], axes=axes_mileage)
            if fig is not None: plt.show()

        # Final energy histogram with colors marking endcond
        fig = None
        if axes_energy is None:
            fig = plt.figure()
            axes_energy = fig.add_subplot(1,1,1)
        if axes_energy != False:
            axes_energy.set_title("Final energy")
            self.plotstate_histogram(
                "energy", xbins=55, weight=True,
                iniend=["e", "i"], log=[True, False], axes=axes_energy)
            if fig is not None: plt.show()

        # Final Rz scatter positions
        fig = None
        if axes_rz is None:
            fig = plt.figure()
            axes_rz = fig.add_subplot(1,1,1)
        if axes_rz != False:
            axes_rz.set_title("Final R-z positions")
            self.plotstate_scatter(
                "R", "z", color="C0", endcond=None,
                iniend=["e", "e", "i", "i"], axesequal=True, axes=axes_rz)
            if fig is not None: plt.show()

        # Final rho-phi scatter positions
        fig = None
        if axes_rhophi is None:
            fig = plt.figure()
            axes_rhophi = fig.add_subplot(1,1,1)
        if axes_rhophi != False:
            axes_rhophi.set_xlim(left=0)
            axes_rhophi.set_ylim([0,360])
            axes_rhophi.set_xticks([0, 0.5, 1.0])
            axes_rhophi.set_yticks([0, 180, 360])
            axes_rhophi.set_title("Final rho-phi positions")
            self.plotstate_scatter(
                "rho", "phimod", color="C0", endcond=None,
                iniend=["e", "e", "i", "i"], axesequal=False, axes=axes_rhophi)
            if fig is not None: plt.show()

    def plotwall_loadvsarea(self, axes=None):
        ids, _, eload, _, _, _, _ = self.getwall_loads()
        area = self.wall.area()[ids-1]
        a5plt.loadvsarea(area, eload, axes=axes)

    def plotwall_3dstill(self, wallmesh=None, points=None, data=None, log=False,
                         cpos=None, cfoc=None, cang=None, axes=None, cax=None):
        """Take a still shot of the mesh and display it using matplotlib
        backend.

        The rendering is done using vtk but the vtk (interactive) window is not
        displayed. It is recommended to use the interactive plot to find desired
        camera position and produce the actual plot using this method. The plot
        is shown using imshow and acts as a regular matplotlib plot.

        Parameters
        ----------
        wallmesh : Polydata
            Mesh representing the wall. If not given then construct the mesh
            from the wall data.
        points : {array_like, bool}, optional
            Array Npoint x 3 defining points (markers) to be shown. For each
            point [x, y, z] coordinates are given. If boolean True is given,
            then markers are read from the endstate.
        cpos : array_like, optional
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional
            Camera angle [azimuth, elevation, roll].
        axes : Axes, optional
            The Axes object to draw on.
        cax : Axes, optional
            The Axes object for the color data (if c contains data), otherwise
            taken from axes.
        """
        if wallmesh is None:
            wallmesh = self.getwall_3dmesh()
        if isinstance(points, bool) and points == True:
            points = self.getstate_pointcloud(endcond="wall")

        (cpos0, cfoc0, cang0) = a5plt.defaultcamera(wallmesh)
        if cpos is None: cpos = cpos0
        if cfoc is None: cfoc = cfoc0
        if cang is None: cang = cang0

        a5plt.still(wallmesh, points=points, data=data, log=log, cpos=cpos,
                    cfoc=cfoc, cang=cang, axes=axes, cax=cax)


    def plotwall_3dinteractive(self, wallmesh=None, *args, points=None,
                               data=None, log=False, cpos=None, cfoc=None,
                               cang=None):
        """Open vtk window to display interactive view of the wall mesh.

        Parameters
        ----------
        wallmesh, optional : Polydata
            Mesh representing the wall. If not given then construct the
            mesh from the wall data.
        *args : tuple (str, method), optional
            Key (str) method pairs. When key is pressed when the plot is
            displayed, the associated method is called. The method should
            take Plotter instance as an argument.
        points : array_like, optional
            Array Npoint x 3 defining points (markers) to be shown. For
            each point [x, y, z] coordinates are given. If boolean True is
            given, then markers are read from the endstate.
        cpos : array_like, optional
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional
            Camera angle [azimuth, elevation, roll].
        """
        if wallmesh is None:
            wallmesh = self.getwall_3dmesh()
        if isinstance(points, bool) and points == True:
            points = self.getstate_pointcloud(endcond="wall")

        (cpos0, cfoc0, cang0) = a5plt.defaultcamera(wallmesh)
        if cpos is None: cpos = cpos0
        if cfoc is None: cfoc = cfoc0
        if cang is None: cang = cang0

        a5plt.interactive(wallmesh, *args, points=points, data=data, log=log,
                          cpos=cpos, cfoc=cfoc, cang=cang)
