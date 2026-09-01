from typing import Optional

import unyt
import numpy as np

from a5py.exceptions import AscotNoDataException
from a5py.data.access import OutputVariant, Leaf
from a5py.data.marker.state import MarkerState
from a5py.data import Orbit, Hist, SimulationOptions

@Leaf.register
class Run(OutputVariant):

    _options: SimulationOptions

    @property
    def options(self):
        return self._options

    def _load(self, file):
        super()._load(file)
        for name in self._file.children:
            if name == "endstate":
                self._diagnostics[name] = MarkerState()
                self._diagnostics[name]._file = self._file.access_data("endstate")
            elif name == "orbit":
                self._diagnostics[name] = Orbit()
                self._diagnostics[name]._file = self._file.access_data("orbit")
            elif "hist" in name:
                self._diagnostics[name] = Hist()
                self._diagnostics[name]._file = self._file.access_data(name)
            elif name == "options":
                file = self._file.access_data("options")
                self._options = SimulationOptions.from_hdf5(file)
                self.options._freeze()

        for i, hist in enumerate(self.options.histograms):
            self._diagnostics[f"hist_{i}"]._params = hist

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
            if not hasattr(self, arg) and arg not in self._diagnostics:
                raise AscotNoDataException(
                    "Data for \"" +  arg + "\" is required but not present.")

    def _save_data(self):
        file = self._file.access_data("options")
        self.options._write_hdf5(file)
        for name, diag in self._diagnostics.items():
            if name == "endstate":
                file = self._file.access_data("endstate")
                diag.save(file)
            if name == "orbit":
                file = self._file.access_data("orbit")
                diag.save(file)
            if "hist" in name:
                file = self._file.access_data(name)
                diag.save(file)

    def getstate(
            self, *qnt: str,
            state: str="ini",
            filter: Optional[list[str | int]]=None,
            mode: Optional[str]=None,
            ) -> np.ndarray | unyt.unyt_array:
        """Evaluate quantities from markers' ini- and endstate.

        Inistate is marker's phase-space position right at the start of the
        simulation and endstate is the position at the end of the simulation.

        This function not only returns the marker phase space coordinates but
        also other quantities that can be inferred from it. For a complete list
        of available quantities, see. The returned arrays are sorted by marker
        ID by default.

        All simulations store both particle and guiding center phase-space
        coordinates. *By default, the returned quantity is in guiding-center
        picture*.

        Parameters
        ----------
        *qnt : str
            Names of the quantities.
        mode : {"gc", "prt"}, *optional*
            Choose whether to return guiding-center or particle quantity.
        state : {"ini", "end"}, *optional*
            Choose whether to return inistate or endstate quantity.
        ids : array_like, *optional*
            Filter markers by their IDs.

            If given, the values in the returned arrays in the same order as
            the given IDs. Cannot be used with ``endcond``.
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
        *val : array_like or tuple(array_like)
            The evaluated quantity sorted by marker ID.

            If multiple quantities are queried, they are returned as a tuple in
            the order they were listed in ``*qnt``.

        Raises
        ------
        ValueError
            Raised when the queried quantity could not be interpreted.
        AscotDataException
            Raised when data required for the operation is not present.
        """
        if state not in ["ini", "end"]:
            raise ValueError("Unrecognized state: " + state)
        self._require("marker")
        #if endcond is not None or state == "end":
        #    self._require("endstate")

        if state == "end":
            state = self._diagnostics["endstate"]
        else:
            state = self.marker

        data = state.evaluate(*qnt, filter=filter, bfield=self["bfield"])
        if len(qnt) == 1:
            data = data[0]
        return data

    def getorbit(
        self, *qnt: str,
        mode: Optional[str]=None,
        filter: Optional[str]=None,
        bfield=None
        ) -> np.ndarray | unyt.unyt_array:

        orbit = self._diagnostics["orbit"]
        state = self["marker"]
        if bfield is None:
            bfield = self["bfield"]
        if mode is None:
            match self.options.simulation.mode:
                case "guiding-center":
                    mode = "guidingcenter"
                case "gyro-orbit":
                    mode = "particle"
                case "field-line":
                    mode = "field-line"

        return orbit.evaluate_quantity(mode, state, bfield, *qnt, filter=filter)
    
    def getdist(self, idx, dimensions):
        return self._diagnostics[f"hist_{idx}"].extract(dimensions)

    def plotorbit_trajectory(
            self,
            x: str,
            y: str,
            z: Optional=None,
            c: Optional=None,
            filter: Optional=None,
            cmap: Optional=None,
            axesequal: bool=False,
            axes: Optional=None,
            cax: Optional=None,
            ) -> None:
        """Plot orbit trajectories in given coordinates.

        The arguments ``x``, ``y``, ``z``, and optionally ```c`` are names of
        the quantities to be plotted. They may also have any of the following
        formats:

        - "diff <quantity>" plots the difference between inistate and current value.
        - "reldiff <quantity>" plots the relative difference (x1/x0 - 1).
        - "log <quantity>" plots the logarithmic value of the
           quantity and "log <diff/reldiff> <quantity>" the logarithm of the
           difference/relative difference.

        The plot is either 2D+1D or 3D+1D, where the extra coordinate is the
        color, depending on the number of plotted quantities.

        Note that the color scale is discrete, not continuous.

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
                "ids", x, y, c, filter=filter)
        elif z is not None and c is not None:
            idarr, xc, yc, zc, cc = self.getorbit(
                "ids", x, y, z, c, endcond=endcond, ids=ids)

        # Find indices to map values from inistate to orbit array
        _, idarr = np.unique(idarr, return_inverse=True)
        idx = np.where(idarr[:-1] != idarr[1:])[0] + 1
        def parsevals(val, log, label, qnt, mode="gc"):
            """Compute values and split them by orbit

            Parameters
            ----------
            val : array_like
                Values to be split.
            log : str
                Logarithmic flag.
            label : str
                Label of the quantity.
            qnt : str
                Name of the quantity.
            mode : {"prt", "gc", "fl"}
                Mode of the quantity.
            """
            if "x/x_0" in label:
                val0 = self.getstate(qnt, state="ini", endcond=endcond, ids=ids, mode=mode)
                val = ( val / val0[idarr] ) - 1
            elif "Delta" in label:
                val0 = self.getstate(qnt, state="ini", endcond=endcond, ids=ids, mode=mode)
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
            zlabel += " [" + str(zc[0].units) + "]"
            bbox = [xmin, xmax, ymin, ymax, zmin, zmax, None, None]
        if c is not None:
            cc, clog, bbox[-2], bbox[-1] = parsevals(cc, clog, clabel, c)
            clabel += " [" + str(cc[0].units) + "]"

        if z is None:
            a5plt.line2d(xc, yc, c=cc, xlog=xlog, ylog=ylog, clog=clog,
                         xlabel=xlabel, ylabel=ylabel, clabel=clabel, bbox=bbox,
                         cmap=cmap, axesequal=axesequal, axes=axes, cax=cax)
        else:
            a5plt.line3d(xc, yc, zc, c=cc, xlog=xlog, ylog=ylog, zlog=zlog,
                         clog=clog, xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
                         clabel=clabel, bbox=bbox, cmap=cmap,
                         axesequal=axesequal, axes=axes, cax=cax)

    def plotorbit_poincare(self, plane, connlen=True, markersize=2,
                           alternative_coords=False, axes=None, cax=None):
        """Create a Poincaré plot where the color separates different markers
        or shows the connection length.

        Parameters
        ----------
        plane : str
            The Poincaré plane to be plotted.

            The argument is expected to have format "pol/tor/rad i" where the
            first part specifies the type of the plane (poloidal, toroidal,
            radial) and the second argument is the plane index. For example,
            ``plane="pol 2"`` plots the second plane with a contant theta angle.
            The order of the planes is same as given by
            :meth:`getorbit_poincareplanes`.
        connlen : bool, optional
            Show connection length and separated lost and confined markers with
            color.

            If true, trajectories of lost markers are colored (in blue
            shades) where the color shows the connection length at that
            position. Confined (or all if conlen=False) markers are shown
            with shades of red where color separates subsequent
            trajectories.
        markersize : int, optional
            Marker size on plot.
        alternative_coords : bool, optional
            Use alternative coordinate axes to visualize the plot.

            The default axes are (R,z) for the poloidal plot and (rho,phi) for
            the toroidal plot. If `alternative_coords`=True, these are changed
            to (rho,theta) and (R,phi), respectively.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.

        Raises
        ------
        ValueError
            Raised when the plane is unknown.
        """
        def choose(default, alternative):
            """Helper function to avoid repetitive if else statements when
            choosing coordinates.
            """
            return alternative if alternative_coords else default

        # Set quantities corresponding to the planetype
        pol, tor, rad = self.getorbit_poincareplanes()
        if "pol" in plane:
            plane = int(plane.split()[1]) - 1
            if plane < len(pol) and plane >= 0:
                plane = pol[plane][1]
                x = choose("rho", "r")
                y = "phimod"
                xlabel = choose("Normalized poloidal flux", "R [m]")
                ylabel = "Toroidal angle [deg]"
                xlim = choose([0, 1.1], None) # None is determined separately
                ylim = [0, 360]
                axesequal = False
            else:
                raise ValueError("Unknown plane.")
        elif "tor" in plane:
            plane = int(plane.split()[1]) - 1
            if plane < len(tor) and plane >= 0:
                plane = tor[plane][1]
                x = choose("r", "rho")
                y = choose("z", "polmod")
                xlabel = choose("R [m]", "Normalized poloidal flux")
                ylabel = choose("z [m]", "Poloidal angle [deg]")
                xlim = choose(None, [0, 1.1]) # None is determined separately
                ylim = choose(None, [0, 360]) # None is determined separately
                axesequal = choose(True, False)
            else:
                raise ValueError("Unknown plane.")
        elif "rad" in plane:
            plane = int(plane.split()[1]) - 1
            if plane < len(rad) and plane >= 0:
                plane = rad[plane][1]
                x = "thetamod"
                y = "phimod"
                xlabel = "Poloidal angle [deg]"
                ylabel = "Toroidal angle [deg]"
                xlim = [0, 360]
                ylim = [0, 360]
                axesequal = True
            else:
                raise ValueError("Unknown plane.")
        else:
            raise ValueError("Unknown plane.")

        plotconnlen = connlen
        ids, x, y, connlen = self.getorbit("ids", x, y, "connlen", pncrid=plane)

        # Determine (R,z) limits for poloidal plane by making sure the data
        # fits nicely and then rounding to nearest decimal.
        if xlim == None:
            xlim = [np.floor(np.amin(x)*9) / 10, np.ceil(np.amax(x)*11) / 10]
        if ylim == None:
            ylim = [np.floor(np.amin(y)*11) / 10, np.ceil(np.amax(y)*11) / 10]

        if plotconnlen:
            # Now set confined markers as having negative connection length
            connlen *= -1
            lost1 = self.getstate("ids", state="end",
                                  endcond=["rhomax", "wall"])

            idx = ~np.in1d(ids, lost1)
            connlen[idx] *= -1
            clabel = "Connection length [" + str(connlen.units) + "]"
        else:
            connlen = None
            clabel = None

        a5plt.poincare(x, y, ids, connlen=connlen, xlim=xlim, ylim=ylim,
                       xlabel=xlabel, ylabel=ylabel, clabel=clabel,
                       markersize=markersize, axesequal=axesequal, axes=axes,
                       cax=cax)

    @classmethod
    def make_fresh(cls, mrk, inputs, params):
        """Initialize a new run that has not been simulated.

        Parameters
        ----------
        mrk : Marker
            Markers to simulate.

            Effectively this argument is the inistate, which is then copied here
            to make an endstate. The endstate markers are simulated.
        inputs
        params

        Returns
        -------
        run : Run
            The new run, ready for simulation.
        """
        run = cls(inputs=inputs)
        run._diagnostics["endstate"] = MarkerState.from_params(mrk._cdata)
        orbit = Orbit.from_params(mrk.n, params.orbit)
        if orbit is not None:
            run._diagnostics["orbit"] = orbit
        for i, h in enumerate(params.histograms):
            run._diagnostics[f"hist_{i}"] = Hist.from_params(h)
        run._options = params
        return run
