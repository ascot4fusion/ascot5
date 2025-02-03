"""Class for BBNBI post-processing.
"""
import numpy as np
import unyt
import a5py.routines.plotting as a5plt
from a5py.routines.distmixin import DistMixin

from a5py.data import Marker, State
from a5py.exceptions import AscotNoDataException

class BBNBIMixin(DistMixin):
    """Mixin class with post-processing results related to AFSI.
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

    def getstate(self, *qnt, mode="gc", ids=None, endcond=None):
        """Evaluate a marker quantity based on its state.

        Parameters
        ----------
        *qnt : str
            Names of the quantities.
        ids : array_like, optional
            Filter markers by their IDs.
        endcond : str or list [str], optional
            Filter markers by their end conditions.

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
        self._require("_state")

        # Get or evaluate the quantity
        data = getattr(self, "_state").get(*qnt, mode=mode)

        # Parse by ids and endcond
        idx = np.ones(data[0].shape, dtype=bool)
        if endcond is not None:
            if not isinstance(endcond, list): endcond = [endcond]

            # Go through each unique end cond and mark that end cond valid or
            # not. This can then be used to make udix as boolean mask array.
            uecs, uidx = np.unique(self._state.get("endcond")[0],
                                   return_inverse=True)
            mask = np.zeros(uecs.shape, dtype=bool)
            for i, uec in enumerate(uecs):
                for ec in endcond:
                    ec = ec.replace(" and ", " ")
                    accept  = State.endcond_check(uec, ec)
                    mask[i] = mask[i] or accept

            idx = mask[uidx]

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self._state.get("ids"), ids))

        for i in range(len(data)):
            data[i] = data[i][idx]
        if "mu" in qnt:
            data[qnt.index("mu")].convert_to_units("eV/T")
        if "psi" in qnt:
            data[qnt.index("psi")].convert_to_units("Wb")
        return data if len(data) > 1 else data[0]

    def getdist(self, dist, exi=False, ekin_edges=None, pitch_edges=None,
                plotexi=False):
        """Return distribution function.

        Parameters
        ----------
        dist : {"5d", "6d", "rho5d", "rho6d", "com"}
            Name of the distribution.

            If ``dist`` is ``None``, returns a list of available distributions.
        exi : bool, optional
            Convert the momentum space to energy-pitch.

            The distribution is normalized to conserve the particle number.

            Not applicable for ``dist`` = {"6d", "rho6d", "com"}.
        ekin_edges : int or array_like, optional
            Number of bins or bin edges in the energy abscissa if ``exi=True``.
        pitch_edges : int or array_like, optional
            Number of bins or bin edges in the pitch abscissa if ``exi=True``.
        plotexi : bool, optional
            Visualize the transformation from ppar-perp to energy-pitch if
            if ``exi=True``.

            Use this option to adjust energy and pitch abscissae to your liking.

        Returns
        -------
        data : :class:`DistData` or list [str]
            The distribution data object or a list of available distributions if
            ``dist`` is ``None``.
        """
        if dist is None:
            dists = ["5d", "6d", "rho5d", "rho6d", "com"]
            for d in dists:
                try:
                    self._require("_dist" + d)
                except AscotNoDataException:
                    dists.remove(d)
            return dists

        if dist == "5d":
            self._require("_dist5d")
            distout = self._dist5d.get()
        elif dist == "6d":
            self._require("_dist6d")
            distout = self._dist6d.get()
        elif dist == "rho5d":
            self._require("_distrho5d")
            distout = self._distrho5d.get()
        elif dist == "rho6d":
            self._require("_distrho6d")
            distout = self._distrho6d.get()
        elif dist == "com":
            self._require("_distcom")
            distout = self._distcom.get()
        else:
            raise ValueError("Unknown distribution")

        mass = np.mean(self.getstate("mass"))
        return self._getdist(distout, mass, exi=exi, ekin_edges=ekin_edges,
                             pitch_edges=pitch_edges, plotexi=plotexi)

    def getstate_markersummary(self):
        """Return a summary of marker end conditions and errors present in
        the data.

        Returns
        -------
        econds : list [(int, str)]
            List of present end conditions.

            Each list member is a tuple, where first item is the number of
            markers with the end condition specified in the second item.

        emsg : list [(str, int, str)]
            List of present errors.

            Each list member is a tuple, where first item is the error
            message, second the line number where error was raised, and
            third is the name of the file.

        Raises
        ------
        AscotNoDataException
            Raised when data required for the operation is not present.
        """
        self._require("_state")
        econd, emsg, emod, eline = self._state.get(
            "endcond", "errormsg", "errormod", "errorline")
        errors = np.unique(np.array([emsg, eline, emod]), axis=1).T

        ec, counts = np.unique(econd, return_counts=True)
        econds = []
        for i, e in enumerate(ec):
            econd = State.endcond_tostring(e)
            if econd == "none": econd = "aborted"
            econds.append( (counts[i], econd) )

        # It would be better to import these via libascot, but then this
        # function would require libascot and we don't want that.
        modules = [
            "mccc_wiener.c", "mccc_push.c", "mccc_coefs", "mccc.c",
            "step_fo_vpa.c", "step_gc_cashkarp.c", "step_gc_rk4", "N0_1D.c",
            "N0_3D.c", "N0_ST.c", "B_3DS.c", "B_2DS.c", "B_STS.c", "B_GS.c",
            "plasma_1D.c", "plasma_1DS.c", "plasma.c", "E_field.c", "neutral.c",
            "E_1DS.c", "B_field.c", "particle.c", "boozer.c", "mhd.c",
            "atomic.c", "asigma.c", "asigma_loc.c",
        ]
        messages = [
            "Input evaluation failed", "Unknown input type",
            "Unphysical quantity when evaluating input",
            "Unphysical marker quantity", "Time step too small/zero/NaN",
            "Wiener array is full or inconsistent",
            "Unphysical result when integrating marker coordinates"
        ]

        emsg = []
        for e in errors:
            if np.sum(e) > 0:
                emsg.append( (messages[e[0]-1], e[1], modules[e[2]-1]) )

        return econds, emsg

    def getstate_markers(self, mrktype, ids=None):
        """Convert state to marker input.

        Parameters
        ----------
        mrktype : {"prt", "gc", "fl"}
            Type of marker input to be created.
        ids : array_like, optional
            Select only these markers for the output.

        Returns
        -------
        mrk : dict
            Markers parameters that can be supplied to :meth:`Prt.write_hdf5`,
            :meth:`GC.write_hdf5` or :meth:`FL.write_hdf5` depending on
            ``mrktype`` value.
        """
        ids = self.getstate("ids", ids=ids)
        mrk = Marker.generate(mrktype, ids.size)
        mrk["ids"] = ids
        if mrktype == "prt":
            qnt = ["r", "phi", "z", "weight", "time", "vr", "vphi", "vz",
                   "mass", "charge", "anum", "znum"]
            state = self.getstate(*qnt, mode="prt", ids=ids)
            for i, q in enumerate(qnt):
                mrk[q] = state[i]
        elif mrktype == "gc":
            qnt = ["r", "phi", "z", "weight", "time", "ekin", "pitch", "zeta",
                   "mass", "charge", "anum", "znum"]
            state = self.getstate(*qnt, mode="gc", ids=ids)
            for i, q in enumerate(qnt):
                if q == "ekin":
                    mrk["energy"] = state[i]
                else:
                    mrk[q] = state[i]
        elif mrktype == "fl":
            qnt = ["r", "phi", "z", "weight", "time", "pitch"]
            state = self.getstate(*qnt, mode="gc", ids=ids)
            for i in qnt:
                mrk[i] = qnt[i]
        return mrk

    def getwall_figuresofmerit(self):
        """Get peak power load and other 0D quantities related to wall loads.

        Returns
        -------
        warea : float
            Total wetted area.
        pload : float
            Peak power load.
        """
        wetted, area, edepo, pdepo, iangle = self.getwall_loads()
        wetted_total = np.sum(area)
        energy_peak  = np.amax(edepo/area)
        return wetted_total, energy_peak

    def getwall_loads(self, weights=True):
        """Get wall loads and associated quantities.

        This method does not return loads on all wall elements (as usually most
        of them receive no loads) but only those that are affected and their
        IDs.

        Parameters
        ----------
        weights : bool, optional
            Include marker weights to get physical results (otherwise particle
            deposition would be just the number of markers that hit the tile).

            Dropping weights is useful to check how many markers hit a tile
            which tells us how good the statistics are.

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
            Angle of incidence (TODO).
        """
        self._require("_state")
        ids, energy, weight = self.getstate("walltile", "ekin", "weight",
                                            endcond="wall")
        energy.convert_to_units("J")
        eunit = (energy.units * weight.units)
        try:
            eunit = (1 * eunit).to("W")
        except unyt.exceptions.UnitConversionError:
            pass
        wetted = np.unique(ids)
        area   = self.wall.area()[wetted-1]
        edepo  = np.zeros(wetted.shape) * eunit
        pdepo  = np.zeros(wetted.shape) * weight.units
        iangle = np.zeros(wetted.shape)

        # Sort markers by wall tile ID to make processing faster
        idx    = np.argsort(ids)
        ids    = ids[idx]
        energy = energy[idx]
        weight = weight[idx]

        idx = np.append(np.argwhere(ids[1:] - ids[:-1]).ravel(), ids.size-1)

        i0 = 0
        for i in range(wetted.size):
            i1 = idx[i] + 1
            pdepo[i]  = np.sum(weight[i0:i1])
            edepo[i]  = np.sum(energy[i0:i1]*weight[i0:i1])
            iangle[i] = 0 #TBD
            i0 = i1

        return wetted, area, edepo, pdepo, iangle

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
        import pyvista as pv
        wallmesh = pv.PolyData( *self.wall.noderepresentation() )
        ids, area, eload, pload, iangle = self.getwall_loads()
        ids = ids - 1 # Convert IDs to indices
        ntriangle = wallmesh.n_faces
        wallmesh.cell_data["pload"]       = np.zeros((ntriangle,)) + np.nan
        wallmesh.cell_data["pload"][ids]  = pload / area
        wallmesh.cell_data["eload"]       = np.zeros((ntriangle,)) + np.nan
        wallmesh.cell_data["eload"][ids]  = eload / area

        #- "iangle" angle of incidence in units of deg
        #wallmesh.cell_data["iangle"]      = np.zeros((ntriangle,)) + np.nan
        #wallmesh.cell_data["iangle"][ids] = iangle
        return wallmesh

    def plotwall_3dstill(self, wallmesh=None, points=None,
                         data=None, log=False, cpos=None, cfoc=None, cang=None,
                         axes=None, cax=None, **kwargs):
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
        data : str, optional
            Name of the cell data in the wall mesh that is shown in color.
        cpos : array_like, optional
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional
            Camera angle [azimuth, elevation, roll].
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        **kwargs
            Keyword arguments passed to :obj:`~pyvista.Plotter`.
        """
        if wallmesh is None:
            wallmesh = self.getwall_3dmesh()
        if isinstance(points, bool) and points == True:
            points = self.getstate_pointcloud(endcond="wall")

        (cpos0, cfoc0, cang0) = a5plt.defaultcamera(wallmesh)
        if cpos is None: cpos = cpos0
        if cfoc is None: cfoc = cfoc0
        if cang is None: cang = cang0

        a5plt.still(wallmesh, points=points, data=data, log=log,
                    cpos=cpos, cfoc=cfoc, cang=cang, axes=axes, cax=cax,
                    **kwargs)

    def plotwall_3dinteractive(self, wallmesh=None, *args, points=None,
                               orbit=None, data=None, log=False,
                               cpos=None, cfoc=None, cang=None, **kwargs):
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
        data : str, optional
            Name of the cell data in the wall mesh that is shown in color.
        cpos : array_like, optional
            Camera position coordinates [x, y, z].
        cfoc : array_like, optional
            Camera focal point coordinates [x, y, z].
        cang : array_like, optional
            Camera angle [azimuth, elevation, roll].
        **kwargs
            Keyword arguments passed to :obj:`~pyvista.Plotter`.
        """
        if wallmesh is None:
            wallmesh = self.getwall_3dmesh()
        if isinstance(points, bool) and points == True:
            points = self.getstate_pointcloud(endcond="wall")

        (cpos0, cfoc0, cang0) = a5plt.defaultcamera(wallmesh)
        if cpos is None: cpos = cpos0
        if cfoc is None: cfoc = cfoc0
        if cang is None: cang = cang0

        a5plt.interactive(wallmesh, *args, points=points, data=data,
                          log=log, cpos=cpos, cfoc=cfoc, cang=cang, **kwargs)

    @a5plt.openfigureifnoaxes(projection=None)
    def plotwall_torpol(self, log=True, clim=None, axes=None, cax=None):
        """Plot the toroidal-poloidal projection of the wall loads.

        Note that the resulting doesn't present the areas accurately and even
        the locations are approximate. The purpose of the plot is to roughly
        illustrate how the losses are distributed.

        Parameters
        ----------
        log : bool, optional
            Make the color scale logarithmic.
        clim : [float, float], optional
            Minimum and maximum values in the color scale.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        """
        wetted, area, edepo, pdepo, iangle = self.getwall_loads()
        d = self.wall.read()
        nelement = wetted.size
        color    = edepo/area
        x1x2x3   = d["x1x2x3"][wetted-1]
        y1y2y3   = d["y1y2y3"][wetted-1]
        z1z2z3   = d["z1z2z3"][wetted-1]

        # Toroidal angle for each vertex
        tor = np.rad2deg( np.arctan2( y1y2y3, x1x2x3 ))

        # Poloidal angle for each vertex
        out = self._root._ascot._eval_bfield(
            1*unyt.m, tor*unyt.deg, 1*unyt.m, 0*unyt.s, evalaxis=True)
        axisr  = out["axisr"].v
        axisz  = out["axisz"].v
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
            log=log, axes=axes, cax=cax, clabel=r"Wall load [W/m$^2$]",
            cmap="Reds",
            xlabel="Toroidal angle [deg]", ylabel="Poloidal angle [deg]")
        axes.set_xticks([0, 90, 180, 270, 360])
        axes.set_yticks([-180, -90, 0, 90, 180])
