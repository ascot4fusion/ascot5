def getwall_figuresofmerit(self, flags=None):
    """Get peak power load and other 0D quantities related to wall loads.

    Parameters
    ----------
    flags : array_like, optional
        Filter output to include only the elements whose flag is in this
        list (the values can either be integers or strings).

    Returns
    -------
    warea : float
        Total wetted area.
    pload : float
        Peak power load.
    """
    _, area, edepo, _, _ = self.getwall_loads(flags=flags)
    wetted_total = np.sum(area)
    energy_peak  = np.amax(edepo/area)
    return wetted_total, energy_peak

def getwall_loads(self, weights=True, p_ids=None, flags=None):
    """Get wall loads and associated quantities.

    This method does not return loads on all wall elements (as usually most
    of them receive no loads) but only those that are affected and their
    IDs.

    For 2D walls, iangle, angle of incidence is not calculated.

    Parameters
    ----------
    weights : bool, optional
        Include marker weights to get physical results (otherwise particle
        deposition would be just the number of markers that hit the tile).

        Dropping weights is useful to check how many markers hit a tile
        which tells us how good the statistics are.
    p_ids : array_like, optional
        Calculate wall loads only for the particles with the given IDs.
    flags : str, int, or array_like, optional
        Filter output to include only the elements whose flag is in this
        list (the values can either be integers or strings).

    Returns
    -------
    ids : array_like
        Indices of loaded wall tiles.
    area : array_like
        Area of each tile.
    edepo : array_like
        Energy/power deposition per tile.
    pdepo : array_like
        Particle/flux deposition per tile.
    iangle : array_like
        Angle of incidence i.e. the angle between particle velocity and
        the surface normal.

        When weights are included, the angle of incidence is the average of
        all markers that hit the tile weighted by the marker energy and
        weight. Otherwise it is the mean value of all markers without any
        weights.
    """
    self._require("_endstate")
    ids, energy, weight, pr, pphi, pz, pnorm, phi = self.getstate(
        "walltile", "ekin", "weight", "pr", "pphi", "pz", "pnorm", "phi",
        state="end", mode="prt", endcond="wall", ids=p_ids)
    if flags is not None:
        w = self.wall.read()
        flag = w["flag"]
        labels = w["labels"]
        if not isinstance(flags, list):
            flags = [flags]
        mask = np.array([
            labels[f] if f in labels else f for f in flags
        ])
        idx = np.where(np.isin(flag[ids-1], mask))[0]
        ids, energy, weight, pr, pphi, pz, pnorm, phi = (
            ids[idx], energy[idx], weight[idx], pr[idx], pphi[idx],
            pz[idx], pnorm[idx], phi[idx]
        )

    energy.convert_to_units("J")
    eunit = (energy.units * weight.units)
    try:
        eunit = (1 * eunit).to("W")
    except unyt.exceptions.UnitConversionError:
        pass
    wetted = np.unique(ids)
    edepo  = np.zeros(wetted.shape)
    pdepo  = np.zeros(wetted.shape)
    iangle = np.zeros(wetted.shape)
    if weights:
        edepo = edepo * eunit
        pdepo = pdepo * weight.units
    else:
        edepo = edepo * energy.units

    area, nvec = self.wall.area(normal=True)
    area = area[wetted-1]
    punit = np.array([pr, pphi, pz]) / pnorm

    # Sort markers by wall tile ID to make processing faster
    idx    = np.argsort(ids)
    ids    = ids[idx]
    energy = energy[idx]
    weight = weight[idx]
    punit  = punit[:,idx].v
    idx = np.append(np.argwhere(ids[1:] - ids[:-1]).ravel(), ids.size-1)

    # Convert unit vectors to cartesian basis
    punit = np.array([
        punit[0,:] * np.cos(phi) - punit[1,:] * np.sin(phi),
        punit[0,:] * np.sin(phi) + punit[1,:] * np.cos(phi),
        punit[2,:]])

    i0 = 0
    for i in range(wetted.size):
        i1 = idx[i] + 1
        vec = nvec[:,wetted[i]-1]
        dotprod = punit[0,i0:i1] * vec[0] + punit[1,i0:i1] * vec[1] \
                            + punit[2,i0:i1] * vec[2]
        angles = np.arccos(dotprod) * (180 / np.pi) * unyt.deg
        angles[angles.v > 90] = 180 * unyt.deg - angles[angles.v > 90]

        if weights: #include weights in the deposition
            pdepo[i]  = np.sum(weight[i0:i1])
            edepo[i]  = np.sum(energy[i0:i1]*weight[i0:i1])
            iangle[i] = np.sum(angles * energy[i0:i1] * weight[i0:i1]) \
                / edepo[i]
        else:
            pdepo[i]  = i1 - i0
            edepo[i]  = np.sum(energy[i0:i1])
            iangle[i] = np.mean(angles)
        i0 = i1

    return wetted, area, edepo, pdepo, iangle

def getwall_3dmesh(self, w_indices=None, p_ids=None):
    """Return 3D mesh representation of 3D wall and associated loads

    Parameters
    ----------
    w_indices : array_like, optional
        List of triangle indices for which the 3D mesh is made.
    p_ids : array_like, optional
        List of particle ids for which the wall loads are calculated.

    Returns
    -------
    wallmesh : Polydata
        Mesh representing the wall.

        The mesh cell data has fields:

        - "pload" particle load in units of prt/m^2 or prt/m^2s,
        - "eload" power/energy load in units of W/m^2 or J/m^2
        - "mload" marker load in units of markers
        - "iangle" angle of incidence (the angle between power flux and
            the surface normal) in deg
        - "label" flag of the wall element
    """
    wallmesh = pv.PolyData(
        *self.wall.noderepresentation(w_indices=w_indices))
    ids, area, eload, pload, iangle = self.getwall_loads(p_ids=p_ids)
    # get the marker load at each tile
    _, _, _, mload, _ = self.getwall_loads(weights=False, p_ids=p_ids)
    ids = ids - 1 # Convert IDs to indices

    w = self.wall.read()
    n_tri, flag = w["nelements"], w["flag"]
    wall_f = np.arange(n_tri, dtype=int)
    if w_indices is not None:
        wall_f = wall_f[w_indices]

    eload_tot = np.zeros((n_tri, )) + np.nan
    eload_tot[ids] = eload/area
    pload_tot = np.zeros((n_tri, )) + np.nan
    pload_tot[ids] = pload/area
    mload_tot = np.zeros((n_tri, )) + np.nan
    mload_tot[ids] = mload
    iangle_tot = np.zeros((n_tri, )) + np.nan
    iangle_tot[ids] = iangle

    wallmesh.cell_data["label"] = flag[wall_f]
    wallmesh.cell_data["pload"] = pload_tot[wall_f]
    wallmesh.cell_data["eload"] = eload_tot[wall_f]
    wallmesh.cell_data["mload"] = mload_tot[wall_f]
    wallmesh.cell_data["iangle"] = iangle_tot[wall_f]
    return wallmesh


def input_plotwallcontour(self, phi=0, axes=None):
    """Plot intersection of the wall and the poloidal plane at the given
    toroidal angle.

    Parameters
    ----------
    phi : float
        Toroidal angle of the poloidal plane.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    """
    ls = self.data.wall.active.getwallcontour(phi=phi)
    line2d(ls[:,:,0], ls[:,:,1], c="black", axesequal=True, axes=axes,
            xlabel="R [m]", ylabel="z [m]")


@a5plt.openfigureifnoaxes(projection=None)
def plotwall_convergence(self, qnt, nmin=1000, nsample=10, axes=None,
                            flags=None):
    """Plot convergence of losses by subsampling the results.

    This function works by picking randomly a subset of
    n (< Total number of markers) from the results, and reweighting those so
    that the total weight of the subset is equal to the original. The subset
    is then used to calculate the figure of merit. Several samples are taken
    with n = np.logspace(nmin, ntotal, nsample) to show how the results
    converge.

    Parameters
    ----------
    qnt : {"lostpower", "peakload"}
        The name of the quantity for which the converge is plotted.
    nmin : int, optional
        Number of markers in the smallest sample.
    nsample : int, optional
        Number of samples.
    flags : array_like, optional
        Filter output to include only the elements whose flag is in this
        list (the values can either be integers or strings).
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    """
    self._require("_endstate")
    ids, weight, ekin = self.getstate("ids", "weight", "ekin", state="end")
    lost = self.getstate("ids", endcond='WALL')
    lost = np.in1d(ids, lost)
    ntotal = ids.size
    nsubset = np.logspace(np.log10(nmin), np.log10(ntotal), nsample)
    nsubset = nsubset.astype('i8')
    wtotal = np.sum(weight)

    val = np.zeros((nsample,))
    rng = np.random.default_rng()
    for i, n in enumerate(nsubset):
        idx = rng.choice(ntotal, replace=False, size=n)
        ids0  = ids[idx]
        lost0 = lost[idx]
        if qnt == 'lostpower':
            if flags is None:
                val[i] = ( np.sum((lost*weight*ekin.to('J'))[idx])
                    * wtotal / np.sum(weight[idx])
                )
            else:
                _, area, loads, _, _ = self.getwall_loads(
                    weights=True, p_ids=ids0[lost0], flags=flags,
                )
                val[i] = (
                    np.sum( loads * area ) * wtotal / np.sum(weight[idx])
                )
            axes.set_ylabel('Lost power [W]')
        elif qnt == 'peakload':
            _, area, loads, _, _ = self.getwall_loads(
                weights=True, p_ids=ids0[lost0], flags=flags,
                )
            val[i] = np.amax(loads/area) * wtotal / np.sum(weight[idx])
            axes.set_ylabel(r'Peak load [W/m$^2$]')
        else:
            raise ValueError(
                f"Unrecognized quantity: {qnt}. Recognized quantities are:"
                f" lostpower, peakload."
                )

    axes.set_xscale('log')
    axes.set_xlabel('Number of markers')
    axes.plot(nsubset, val)

def plotwall_loadvsarea(self, flags=None, axes=None):
    """Plot histogram showing area affected by at least a given load.

    Parameters
    ----------
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    """
    _, area, eload, _, _ = self.getwall_loads(flags=flags)
    a5plt.loadvsarea(area, eload/area, axes=axes)


@a5plt.openfigureifnoaxes(projection=None)
def plotwall_poloidalslice(self, r, z, data=None, axes=None):
    """Project points onto the closest segments and sum the values for each
    segment.

    Parameters
    ----------
    r : array-like
        R-coordinates of the curve.
    z : array-like
        z-coordinates of the curve.
    """
    w = self.wall.read()
    ri = np.sqrt(  np.mean(w["x1x2x3"], axis=1)**2
                    + np.mean(w["y1y2y3"], axis=1)**2 )
    zi = np.mean(w["z1z2z3"], axis=1)
    ri = ri[data[0]-1]
    zi = zi[data[0]-1]
    ci = data[1]
    n_segments = len(r) - 1
    c = np.zeros(n_segments)

    # Midpoints of each segment
    midpoints = np.column_stack([(r[:-1] + r[1:]) / 2, (z[:-1] + z[1:]) / 2])

    # Build KDTree for efficient nearest neighbor search
    tree = cKDTree(midpoints)
    # Find closest segment midpoint for each point
    distances, indices = tree.query(np.column_stack([ri, zi]))

    # Sum ci values to the corresponding segment
    np.add.at(c, indices, ci)
    a5plt.line2d([r], [z], c=[c],
                    bbox=[np.amin(r), np.amax(r), np.amin(z), np.amax(z),
                        0, np.amax(c)], axes=axes)
    print(c)


@a5plt.openfigureifnoaxes(projection=None)
def plotwall_2D_parametrized(self, axes=None, particle_load=False,
                                ref_indices=None, colors=None, xlabel=None,
                                ylabel=None, ref_cmap='jet'):
    """Plot heat loads as a function of the parametrized wall contour.

    Parameters
    ----------
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    particle_load : bool, optional
        Plot particle fluxes instead of powerloads.
    ref_indices : {array_like, int}, optional
        Plot 2D wall points given by ref_indices.
    colors : arraylike, optional
        Colors of all drawn indices. Must be same length as ref_indices.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    ref_cmap : str, optional
        Name of the colormap for ref_indices.
    """

    wetted, area, edepo, pdepo, iangle = self.getwall_loads()
    xdif = self.wall.getnormedline().flatten()
    x_orig = np.cumsum(xdif)
    x = x_orig - xdif/2 #initialize loads at half point
    #x = x.flatten()
    y = np.zeros(x.shape)
    if particle_load:
        y[wetted] = pdepo / area
    else:
        y[wetted] = edepo / area


    # wrap x and y values
    x = np.r_[-xdif[-1]/2, x, x_orig[-1]+xdif[0]/2]
    y = np.r_[y[-1],       y, y[0]]

    if xlabel is None:
        xlabel = 's ['+str(xdif.units)+']'
    if ylabel is None:
        if particle_load:
            ylabel = 'Particleload ['+str(pdepo.units/area.units)+']'
        else:
            ylabel = 'Heatload ['+str(edepo.units/area.units)+']'
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)

    axes.plot(x, y)

    #set axes limits
    xmax = x_orig[-1] #for limiting axes
    xmin = 0
    ymax = np.max(y)
    ymin = 0

    axes.set_xlim(xmin, xmax)
    axes.set_ylim(ymin, ymax*1.1) #10% padding before top


    if ref_indices is not None:

        x_ori = np.r_[0, x_orig.flatten()]
        dx    = np.r_[xdif.flatten()/2, xdif.flatten()[-1]/2]
        dx2   = x[1:]-x[:-1]
        dy   = y[1:]-y[:-1]
        y_ori = y[1:]- dy*(dx/dx2)

        xp = x_ori[ref_indices]
        yp = y_ori[ref_indices]

        x_max = x_ori[-1]

        if colors is None:
            N_ref = len(ref_indices)
            import matplotlib as mp
            colors = mp.cm.get_cmap(ref_cmap)(np.linspace(0,1,N_ref))
        else:
            if colors.shape[0] != len(ref_indices):
                raise ValueError('ref_indices and colors len do not match')

        if 0 in ref_indices:
            ind = ref_indices.index(0)
            xp     = np.append(xp, xmax)
            yp     = np.append(yp, yp[ind])
            colors = np.append(colors, colors[ind, :][np.newaxis, :], axis=0)

        axes.scatter(xp, yp, s=70, c=colors)

@a5plt.openfigureifnoaxes(projection=None)
def plotwall_2D_contour(self, axes=None, particle_load=False,
                        ref_indices=None, colors=None, clog="linear",
                        cmap=None, ref_cmap='jet', xlabel=None, ylabel=None,
                        clabel=None):
    """Plot heat loads as in the 2D wall contour with a colorbar.

    Parameters
    ----------
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    particle_load : bool, optional
        Plot particle fluxes instead of powerloads.
    ref_indices : {array_like, int}, optional
        Plot 2D wall points given by ref_indices.
    colors : arraylike, optional
        Colors of all drawn indices. Must be same length as ref_indices.
    clog : {"linear", "log", "symlog"}, optional
        color-axis scaling.
    cmap : str, optional
        Name of the colormap.
    ref_cmap : str, optional
        Name of the colormap for ref_indices.
    xlabel : str, optional
        Label for the x-axis.
    ylabel : str, optional
        Label for the y-axis.
    clabel : str, optional
        Label for the color axis.
    """
    wetted, area, edepo, pdepo, iangle = self.getwall_loads()
    lines = self.wall.getwallcontour()
    y = np.zeros((lines.shape[0], ))
    if particle_load:
        y[wetted] = pdepo / area
    else:
        y[wetted] = edepo / area

    if xlabel is None:
        xlabel = 'R [m]'
    if ylabel is None:
        ylabel = 'z [m]'
    if clabel is None:
        if particle_load:
            clabel = 'Particleload ['+str(pdepo.units/area.units)+']'
        else:
            clabel = 'Heatload ['+str(edepo.units/area.units)+']'
    bbox = [np.min(lines[:, :, 0]), np.max(lines[:, :, 0]),\
            np.min(lines[:, :, 1]), np.max(lines[:, :, 1]),\
            np.min(y),              np.max(y)]

    import matplotlib as mp
    if cmap is None:
        cmap ='viridis'
    cols = mp.cm.get_cmap(cmap)((y - bbox[-2])/(bbox[-1] - bbox[-2]))
    a5plt.line2d(x=lines[:, :, 0], y=lines[:, :, 1], c=cols, clog=clog,\
                xlabel=xlabel, ylabel=ylabel, clabel=clabel, bbox=bbox,\
                cmap=cmap, axesequal=True, axes=axes)

    if ref_indices is not None:

        w =self.wall.read()
        x_ori = w['r']
        y_ori = w['z']

        xp = x_ori[ref_indices]
        yp = y_ori[ref_indices]

        if colors is None:
            N_ref = len(ref_indices)
            colors = mp.cm.get_cmap(ref_cmap)(np.linspace(0,1,N_ref))
        else:
            if colors.shape[0] != len(ref_indices):
                raise ValueError('ref_indices and colors len do not match')

        axes.scatter(xp, yp, s=70, c=colors)

@a5plt.openfigureifnoaxes(projection=None)
def plotwall_torpol(self, qnt='eload', getaxis=None,log=True, clim=None,
                    cmap=None, axes=None, cax=None):
    """Plot the toroidal-poloidal projection of the wall loads.

    Note that the resulting doesn't present the areas accurately and even
    the locations are approximate. The purpose of the plot is to roughly
    illustrate how the losses are distributed.

    Parameters
    ----------
    qnt : {'eload', 'pload', 'iangle', 'label'} or tuple [int, float],
    optional
        Quantity to plot.

        One can also provide a tuple instead to plot a quantity that is
        not predefined (e.g. sputtering yield). In this case the tuple
        should consists of wallids and the associated value of the quantity
        for that tile.
    getaxis : (float, float) or callable
        Location of the magnetic or geometrical axis which is used to
        determine the poloidal angle.

        If this parameter is None (default), the magnetic axis evaluated via
        libascot is used (requires that bfield is initialized).

        For tokamaks, a tuple of two floats corresponding to the axis (R,z)
        coordinates can be given explicitly.

        For stellarators, one can provide a function that takes one
        positional argument ``phi`` (`float`), which is the toroidal angle
        in degress, and returns a tuple of axis (R,z) coordinates at that
        toroidal angle.

    log : bool, optional
        Make the color scale logarithmic.
    clim : [float, float], optional
        Minimum and maximum values in the color scale.
    cmap : str, optional
        Colormap to be used.

        Defaults to 'Reds' for 'eload' and 'pload' and 'viridis' for
        'iangle'.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    """
    d = self.wall.read()
    if qnt == 'label':
        nelement = d["nelements"]
        x1x2x3 = d["x1x2x3"]
        y1y2y3 = d["y1y2y3"]
        z1z2z3 = d["z1z2z3"]
        color = d["flag"].ravel()
        clabel = r"Label"
        if cmap is None: cmap = 'viridis'
    elif isinstance(qnt, str):
        wetted, area, edepo, pdepo, iangle = self.getwall_loads()
        nelement = wetted.size
        x1x2x3 = d["x1x2x3"][wetted-1]
        y1y2y3 = d["y1y2y3"][wetted-1]
        z1z2z3 = d["z1z2z3"][wetted-1]
        if qnt == 'eload':
            color = edepo/area
            clabel = r"Wall load [W/m$^2$]"
            if cmap is None: cmap = 'Reds'
        elif qnt == 'pload':
            color = pdepo/area
            clabel = r"Particle flux [prt s$^{-1}$m$^{-2}$]"
            if cmap is None: cmap = 'Reds'
        elif qnt == 'iangle':
            color = iangle
            clabel = r"Angle of incidence [deg]"
            if cmap is None: cmap = 'viridis'
    else:
        wetted = qnt[0]
        nelement = wetted.size
        x1x2x3 = d["x1x2x3"][wetted-1]
        y1y2y3 = d["y1y2y3"][wetted-1]
        z1z2z3 = d["z1z2z3"][wetted-1]
        color = np.array(qnt[1])
        clabel = r""
        if cmap is None: cmap = 'Reds'

    # Toroidal angle for each vertex
    tor = np.rad2deg( np.arctan2( y1y2y3, x1x2x3 ))

    # Poloidal angle for each vertex
    if getaxis is None:
        out = self._root._ascot._eval_bfield(
            1*unyt.m, tor*unyt.deg, 1*unyt.m, 0*unyt.s, evalaxis=True)
        axisr  = out["axisr"].v
        axisz  = out["axisz"].v
    elif getaxis is tuple and len(getaxis) == 2:
        axisr = getaxis[0]
        axisz = getaxis[1]
    elif getaxis is callable:
        axisr, axisz = getaxis(tor*unyt.deg)

    rmajor = np.sqrt( x1x2x3**2 + y1y2y3**2 )
    pol    = np.rad2deg( np.arctan2(z1z2z3 - axisz, rmajor - axisr) )

    # If triangle would otherwise be split on different ends of [0, 360]
    # range in tor or pol, move it to top/right.
    dx = (np.amax(tor, axis=1) - np.amin(tor, axis=1)) > 180.0
    tor[dx, :] = np.mod(tor[dx, :], 360)
    dy = (np.amax(pol, axis=1) - np.amin(pol, axis=1)) > 180.0
    pol[dy, :] = np.mod(pol[dy, :], 360)

    # Move the horizontal axis to 0...360
    tor[np.amin(tor, axis=1) < 0, :] += 360.0

    # Sort the patches so that triangles closest to the axis are drawn
    # last (so that they are on top of others)
    rminor = np.sqrt( (rmajor - axisr)**2 + (z1z2z3 - axisz)**2 )
    idx    = np.flipud( np.argsort(np.amin(rminor, axis=1), axis=0) )
    tor    = tor[idx, :]
    pol    = pol[idx, :]
    color  = color[idx]

    # Convert the data into format [vert0, vert1, ...] where verti contains
    # vertices of triangle i in (3,2) array.
    patches = np.zeros((nelement,3,2))
    patches[:,:,0] = tor[:,:]
    patches[:,:,1] = pol[:,:]
    a5plt.triangularpatch(
        patches, color, xlim=[-2, 362], ylim=[-182, 182], clim=clim,
        log=log, axes=axes, cax=cax, clabel=clabel,
        cmap=cmap,
        xlabel="Toroidal angle [deg]", ylabel="Poloidal angle [deg]")
    axes.set_xticks([0, 90, 180, 270, 360])
    axes.set_yticks([-180, -90, 0, 90, 180])

def plotwall_3dstill(
        self, wallmesh=None, points=None, orbit=None, data=None, log=False,
        clim=None, cpos=None, cfoc=None, cang=None, p_ids=None,
        w_indices=None, axes=None, cax=None, **kwargs):
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
    orbit : int, optional
        ID of a marker whose orbit is plotted.
    data : str, optional
        Name of the cell data in the wall mesh that is shown in color.
    log : bool, optional
        Color range is logarithmic if True.
    clim : [float, float], optional
        Color [min, max] limits.
    cpos : array_like, optional
        Camera position coordinates [x, y, z].
    cfoc : array_like, optional
        Camera focal point coordinates [x, y, z].
    cang : array_like, optional
        Camera angle [azimuth, elevation, roll].
    p_ids : array_like, optional
        List of ids of the particles for which the heat load is shown.
    w_indices : array_like, optional
        List of wall indices which are included in the wall mesh.
    axes : :obj:`~matplotlib.axes.Axes`, optional
        The axes where figure is plotted or otherwise new figure is created.
    cax : :obj:`~matplotlib.axes.Axes`, optional
        The color bar axes or otherwise taken from the main axes.
    **kwargs
        Keyword arguments passed to :obj:`~pyvista.Plotter`.
    """
    if wallmesh is None:
        wallmesh = self.getwall_3dmesh(p_ids=p_ids, w_indices=w_indices)
    if isinstance(points, bool) and points == True:
        points = self.getstate_pointcloud(endcond="wall")
    if orbit is not None:
        x,y,z = self.getorbit("x", "y", "z", ids=orbit)
        orbit = np.array([x,y,z]).T

    (cpos0, cfoc0, cang0) = a5plt.defaultcamera(wallmesh)
    if cpos is None: cpos = cpos0
    if cfoc is None: cfoc = cfoc0
    if cang is None: cang = cang0

    a5plt.still(wallmesh, points=points, data=data, orbit=orbit, log=log,
                clim=clim, cpos=cpos, cfoc=cfoc, cang=cang, axes=axes,
                cax=cax, **kwargs)

def plotwall_3dinteractive(self, wallmesh=None, *args, points=None,
                            orbit=None, data=None, log=False, clim=None,
                            cpos=None, cfoc=None, cang=None,
                            p_ids=None, w_indices=None, **kwargs):
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
    orbit : int, optional
        ID of a marker whose orbit is plotted.
    data : str, optional
        Name of the cell data in the wall mesh that is shown in color.
    log : bool, optional
        Color range is logarithmic if True.
    clim : [float, float], optional
        Color [min, max] limits.
    cpos : array_like, optional
        Camera position coordinates [x, y, z].
    cfoc : array_like, optional
        Camera focal point coordinates [x, y, z].
    cang : array_like, optional
        Camera angle [azimuth, elevation, roll].
    p_ids : array_like, optional
        List of ids of the particles for which the heat load is shown.
    w_indices : array_like, optional
        List of wall indices which are included in the wall mesh.
    **kwargs
        Keyword arguments passed to :obj:`~pyvista.Plotter`.
    """
    if wallmesh is None:
        wallmesh = self.getwall_3dmesh(p_ids=p_ids, w_indices=w_indices)
    if isinstance(points, bool) and points == True:
        points = self.getstate_pointcloud(endcond="wall")
    if orbit is not None:
        x,y,z = self.getorbit("x", "y", "z", ids=orbit)
        orbit = np.array([x,y,z]).T

    (cpos0, cfoc0, cang0) = a5plt.defaultcamera(wallmesh)
    if cpos is None: cpos = cpos0
    if cfoc is None: cfoc = cfoc0
    if cang is None: cang = cang0

    a5plt.interactive(wallmesh, *args, points=points, data=data,
                        orbit=orbit, log=log, clim=clim,
                        cpos=cpos, cfoc=cfoc, cang=cang, **kwargs)

