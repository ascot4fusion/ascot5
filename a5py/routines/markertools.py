import numpy as np

def filterdata(data, qnt, ids=None, endcond=None):
    # Parse by ids and endcond
    idx = np.ones(data[0].shape, dtype=bool)
    if endcond is not None:
        if not isinstance(endcond, list): endcond = [endcond]

        # Go through each unique end cond and mark that end cond valid or
        # not. This can then be used to make udix as boolean mask array.
        uecs, uidx = np.unique(self._endstate.get("endcond")[0],
                                return_inverse=True)
        mask = np.zeros(uecs.shape, dtype=bool)
        for i, uec in enumerate(uecs):
            for ec in endcond:
                ec = ec.replace(" and ", " ")
                accept  = State.endcond_check(uec, ec)
                mask[i] = mask[i] or accept

        idx = mask[uidx]

    if ids is not None:
        idx = np.logical_and(idx, np.in1d(self._inistate.get("ids"), ids))

    for i in range(len(qnt)):
        data[i] = data[i][idx]
    if "mu" in qnt:
        data[qnt.index("mu")].convert_to_units("eV/T")
    if "psi" in qnt:
        data[qnt.index("psi")].convert_to_units("Wb")
    return data if len(data) > 1 else data[0]


def filterorbit():
    idx = np.ones(data[0].shape, dtype=bool)
    if endcond is not None:
        eids = self.getstate("ids", endcond=endcond)
        idx = np.logical_and(idx, np.in1d(idarr, eids))

    if pncrid is not None:
        pncridarr = self._orbit.get(self._inistate, self._endstate,
                                    "pncrid")[0]
        idx = np.logical_and(idx, np.in1d(pncridarr, pncrid))

    if ids is not None:
        idx = np.logical_and(idx, np.in1d(idarr, ids))

    for i in range(len(data)):
        data[i] = data[i][idx]
    if "mu" in qnt:
        data[qnt.index("mu")].convert_to_units("eV/T")
    if "psi" in qnt:
        data[qnt.index("psi")].convert_to_units("Wb")
    return data if len(data) > 1 else data[0]


def plotstate_scatter(self, x, y, z=None, c=None, xmode="gc", ymode="gc",
                          zmode="gc", cmode="gc", endcond=None, ids=None,
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
    xmode : {"prt", "gc"}, optional
        Evaluate x quantity in particle or guiding-center phase space.
    ymode : {"prt", "gc"}, optional
        Evaluate y quantity in particle or guiding-center phase space.
    zmode : {"prt", "gc"}, optional
        Evaluate z quantity in particle or guiding-center phase space.
    cmode : {"prt", "gc"}, optional
        Evaluate color quantity in particle or guiding-center phase space.
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

    cc = None
    clog = False
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

def plotstate_histogram(self, x, y=None, xbins=10, ybins=10, xmode="gc",
                        ymode="gc", endcond=None, ids=None, weight=False,
                        logscale=False, cmap=None, axesequal=False,
                        axes=None, cax=None):
    """Make a histogram plot of marker state coordinates.

    The histogram is either 1D or 2D depending on if the y coordinate is
    provided. In the 1D histogram the markers with different endstate are
    separated by color if the endcond argument is None.

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
    y : str, optional
        Name of the quantity on y-axis.

        If not None, the histogram will be in 2D.
    xbins : int or array_like, optional
        Bin edges for the x coordinate or the number of bins.
    ybins : int or array_like, optional
        Bin edges for the y coordinate or the number of bins.
    xmode : {"prt", "gc"}, optional
        Evaluate x quantity in particle or guiding-center phase space.
    ymode : {"prt", "gc"}, optional
        Evaluate y quantity in particle or guiding-center phase space.
    endcond : str, array_like, optional
        Endcond of those markers which are plotted.

        If None and the plot is in 1D, the histogram will be stacked with
        different colors indicating different end conditions.
    ids : array_like
        IDs of the markers to be plotted.
    weight : bool
        Whether to weight markers when they are binned to the histogram.

        The weighted histogram represents physical particles whereas
        the unweighted histogram corresponds to markers.
    cmap : str, optional
        Name of the colormap used in the 2D histogram.
    axesequal : bool, optional
        Flag whether to set aspect ratio for x and y axes equal in 2D.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    def parsearg(arg, mode, endcond):
        """Parse arguments of from "log ini qnt" to (qnt, value, islog).
        """
        arg = arg.lower()
        log = "linear"
        if "log" in arg:
            arg = arg.replace("log", "")
            log = "log"

        if "ini" in arg:
            arg = arg.replace("ini", "").strip()
            val = self.getstate(arg, state="ini", endcond=endcond, ids=ids,
                                mode=mode)
            if log == "log": arg = "|" + arg + "|"
            arg = "Initial " + arg
        elif "end" in arg:
            arg = arg.replace("end", "").strip()
            val = self.getstate(arg, state="end", endcond=endcond, ids=ids,
                                mode=mode)
            if log == "log": arg = "|" + arg + "|"
            arg = "Final " + arg
        elif "reldiff" in arg:
            arg = arg.replace("reldiff", "").strip()
            val1 = self.getstate(arg, state="ini", endcond=endcond, ids=ids,
                                    mode=mode)
            val2 = self.getstate(arg, state="end", endcond=endcond, ids=ids,
                                    mode=mode)
            val = (val2 -val1) / val1
            arg = r"$\Delta x/x_0$ " + arg
            if log == "log": arg = "|" + arg + "|"
        elif "diff" in arg:
            arg = arg.replace("diff", "").strip()
            val1 = self.getstate(arg, state="ini", endcond=endcond, ids=ids,
                                    mode=mode)
            val2 = self.getstate(arg, state="end", endcond=endcond, ids=ids,
                                    mode=mode)
            val = val2 - val1
            arg = r"$\Delta$ " + arg
            if log == "log": arg = "|" + arg + "|"
        else:
            raise ValueError(
                "Unclear if a quantity is evaluated from ini or endstate.\n"
                + "Use \"ini/end/diff/reldiff %s\"" % (arg)
            )
        if log == "log": val = np.abs(val)
        # Make sure val is an unyt array
        try:
            val.units
        except AttributeError:
            val *= unyt.dimensionless
        return (arg, val, log)

    if y is None:
        # 1D plot
        xcs = []
        weights = []

        #ecs, _ = self.getstate_markersummary()
        #for ec in ecs:
        #    if endcond is not None and ec[1] not in [endcond]:
        #        continue

        #    ids, w = self.getstate("ids", "weight", endcond=ec[1], ids=ids)
        #    xcs.append(xc[ids-1])
        #    weights.append(w)

        #xc = [xc]

        x0       = x
        xcs      = []
        weights  = []
        endconds = []
        ecs, _ = self.getstate_markersummary()
        for ec in ecs:
            if endcond is None or ec[1] in endcond:
                w = self.getstate("weight", endcond=ec[1], ids=ids)
                x, xc, xlog = parsearg(x0, xmode, ec[1])
                xcs.append(xc)
                weights.append(w)
                endconds.append(ec)

        if len(xcs) == 0: return

        # Sort data so that when the stacked histogram is plotted, the stack
        # with most markers is at the bottom.
        idx = np.argsort([len(i) for i in xcs])[::-1]
        xcs = [xcs[i] for i in idx]
        ecs = [endconds[i][1] + " : %.2e" % endconds[i][0] for i in idx]
        weights = [weights[i] for i in idx]
        if not weight: weights = None
        a5plt.hist1d(x=xcs, xbins=xbins, weights=weights, xlog=xlog,
                        logscale=logscale, xlabel=x, axes=axes, legend=ecs)

    else:
        # 2D plot
        x, xc, xlog = parsearg(x, xmode, endcond)
        y, yc, ylog = parsearg(y, ymode, endcond)
        weights = self.getstate("weight", state="ini", endcond=endcond,
                                ids=ids)
        if not weight: weights = None
        a5plt.hist2d(xc, yc, xbins=xbins, ybins=ybins, weights=weights,
                        xlog=xlog, ylog=ylog, logscale=logscale, xlabel=x,
                        ylabel=y, axesequal=axesequal, axes=axes, cax=cax)

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
        ``plane="pol 2"`` plots the second poloidal plane. The order of
        the planes is same as given by :meth:`getorbit_poincareplanes`.
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
            x = choose("r", "rho")
            y = choose("z", "polmod")
            xlabel = choose("R [m]", "Normalized poloidal flux")
            ylabel = choose("z [m]", "Poloidal angle [deg]")
            xlim = choose(None, [0, 1.1]) # None is determined separately
            ylim = choose(None, [0, 360]) # None is determined separately
            axesequal = choose(True, False)
        else:
            raise ValueError("Unknown plane.")
    elif "tor" in plane:
        plane = int(plane.split()[1]) - 1
        if plane < len(tor) and plane >= 0:
            plane = tor[plane][1]
            x = choose("rho", "r")
            y = "phimod"
            xlabel = choose("Normalized poloidal flux", "R [m]")
            ylabel = "Toroidal angle [deg]"
            xlim = choose([0, 1.1], None) # None is determined separately
            ylim = [0, 360]
            axesequal = False
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