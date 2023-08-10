"""Wrapper class for RunNode that adds methods to plot and access results.

The motivation for this class is following. Consider that user wants to plot
a figure that shows i) particle orbits in (R,z), ii) their final (R,z) position,
and iii) the wall tiles. It is not clear whether this plot should be implemented
as a method for State or Orbit class as data from both are needed. Furthermore,
neither has access to the wall data. The only logical place for this method
therefore is in the RunNode that has access to all relevant data.
"""
import numpy as np
import warnings
try:
    import pyvista as pv
except ImportError:
    warnings.warn("Could not import pyvista. 3D wall plotting disabled.")
import unyt

from a5py.exceptions import AscotNoDataException

import a5py.routines.plotting as a5plt
import a5py.wall as wall
from a5py.ascot5io import Marker, State, Orbits, Dist
from a5py.ascot5io.dist import DistMoment
import a5py.physlib as physlib

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
            pncridarr = self._orbit.get(self._inistate, self._endstate,
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
        self._require("_endstate")
        econd, emsg, emod, eline = self._endstate.get(
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

        emsg = []
        for e in errors:
            if np.sum(e) > 0:
                emsg.append( (messages[e[0]-1], e[1], modules[e[2]-1]) )

        return econds, emsg

    def getstate_losssummary(self):
        """Return a summary of lost markers.
        """
        self._require("_endstate")

        wmrks, emrks = self.getstate("weight", "ekin", state="end")
        wloss, eloss = self.getstate("weight", "ekin", state="end",
                                     endcond="wall")

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
        endcond : str, optional
            Only return markers that have given end condition.
        """
        self._require("_endstate")
        return np.array([self._endstate.get("x", endcond=endcond),
                         self._endstate.get("y", endcond=endcond),
                         self._endstate.get("z", endcond=endcond)]).T

    def getstate_markers(self, mrktype, ids=None):
        """Convert endstate to marker input.

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
        ids = self.getstate("ids", state="ini", ids=ids)
        mrk = Marker.generate(mrktype, ids.size)
        mrk["ids"] = ids
        if mrktype == "prt":
            qnt = ["r", "phi", "z", "weight", "time", "vr", "vphi", "vz",
                   "mass", "charge", "anum", "znum"]
            state = self.getstate(*qnt, mode="prt", state="end", ids=ids)
            for i in qnt:
                mrk[i] = qnt[i]
        elif mrktype == "gc":
            qnt = ["r", "phi", "z", "weight", "time", "ekin", "pitch", "zeta",
                   "mass", "charge", "anum", "znum"]
            state = self.getstate(*qnt, mode="gc", state="end", ids=ids)
            for i, q in enumerate(qnt):
                if q == "ekin":
                    mrk["energy"] = state[i]
                else:
                    mrk[q] = state[i]
        elif mrktype == "fl":
            qnt = ["r", "phi", "z", "weight", "time", "pitch"]
            state = self.getstate(*qnt, mode="gc", state="end", ids=ids)
            for i in qnt:
                mrk[i] = qnt[i]
        return mrk

    def getorbit_poincareplanes(self):
        """Return a list of Poincaré planes that were present in the simulation
        and the corresponding pncrid.

        Returns
        -------
        pol : list [(float, int)]
            List of tuples with a toroidal angle and pncrid for each poloidal
            plane or empty list if there are none.
        tor : list [float]
            List of tuples with a poloidal angle and pncrid for each toroidal
            plane or empty list if there are none.
        rad : list [float]
            List of tuples with radius and pncrid for radial surfaces or empty
            list if there are none.
        """
        opt  = self.options.read()
        hasdata = opt["ENABLE_ORBITWRITE"] == 1 and opt["ORBITWRITE_MODE"] == 0
        if not hasdata: return None

        polval = opt["ORBITWRITE_POLOIDALANGLES"]
        torval = opt["ORBITWRITE_TOROIDALANGLES"]
        radval = opt["ORBITWRITE_RADIALDISTANCES"]

        pncrid = 0
        def opt2list(val):
            """Convert the option string to a list with tuples (val, pncrid)
            """
            nonlocal pncrid
            tuples = []
            if val[0] < 0: val = []
            for i in range(len(val)):
                tuples.append((val[i], pncrid))
                pncrid += 1
            return tuples

        return opt2list(polval), opt2list(torval), opt2list(radval)

    def getorbit_average(self, qnt, ids):
        """Calculate average of a quantity during a single poloidal transit.

        Parameters
        ----------
        qnt : str
            Name of the averaged quantity.
        
        """
        qnt, mileage, r, z, p, pitch, pol = \
            self.getorbit(qnt, "mileage", "r", "z", "phi", "pitch", "theta",
                          ids=ids)
        if any(pitch < 0) and any(pitch > 0):
            if pitch[0] < 0 or pitch[1] < 0:
                i1 = np.argmax(pitch > 0)
                i2 = np.argmax(pitch[i1:] < 0) + i1
                i2 = np.argmax(pitch[i2:] > 0) + i2
            else:
                i1 = np.argmax(pitch < 0)
                i2 = np.argmax(pitch[i1:] > 0) + i1
                i2 = np.argmax(pitch[i2:] < 0) + i2
        else:
            i1 = 1
            if pol[-1] > pol[0]:
                i2 = np.argmax(pol > (pol[0] + 360*unyt.deg))
            else:
                i2 = np.argmax(pol < (pol[0] - 360*unyt.deg))
        qnt = qnt[i1:i2]
        mileage = mileage[i1-1:i2]
        r = r[i1-1:i2]
        z = z[i1-1:i2]
        p = p[i1-1:i2]
        x, y, _ = physlib.pol2cart(r, p)
        ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2 + np.diff(z)**2)
        qnt = np.sum(qnt*ds) / np.sum(ds)
        return mileage[1:], r, z, qnt

    def getwall_figuresofmerit(self):
        """
        """
        self._require("_endstate")
        ids, energy, weight = self.getstate("walltile", "ekin", "weight",
                                            state="end", endcond="wall")
        area   = self.wall.area()

        wetted_total, energy_peak = wall.figuresofmerit(ids, energy, weight,
                                                        area)
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
        msg += ["Total wetted area: " + str(np.around(wetted_total, decimals=2))
                + r" $m^2$"]
        msg += ["Peak load: " + str(np.around(energy_peak, decimals=2)) + " "
                + str(unit) + r"$/m^2$"]
        return msg

    def getwall_loads(self):
        """Get 3D wall loads and associated quantities.

        This method does not return loads on all wall elements (as usually most
        of them receive no loads) but only those that are affected and their
        IDs.

        Returns
        -------
        ids : array_like
            a
        edepo : array_like
            a
        eload : array_like
            a
        pdepo : array_like
            a
        pload : array_like
            a
        mdepo : array_like
            a
        iangle : array_like
            a
        """
        self._require("_endstate")
        ids, energy, weight = self.getstate("walltile", "ekin", "weight",
                                            state="end", endcond="wall")
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

    def getdist(self, dist, exi=False, ekin_edges=None, pitch_edges=None):
        """Return distribution function.

        Parameters
        ----------
        dist : {"5d", "6d", "rho5d", "rho6d", "com"}
            Name of the distribution.

            If ``dist`` is ``None``, returns a list of available distributions.
        exi : bool, optional
            Convert the momentum space to energy-pitch.

            Not applicable for ``dist`` = {"6d", "rho6d", "com"}.

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

        mass = np.mean(self.getstate("mass"))
        if dist == "5d":
            self._require("_dist5d")
            dist = self._dist5d.get()
            if exi:
                if ekin_edges is None:
                    ekin_edges  = int(dist.abscissa("ppar").size / 2)
                if pitch_edges is None:
                    pitch_edges = dist.abscissa("pperp").size
                dist = Dist.ppappe2ekinpitch(
                    dist, mass, ekin_edges=ekin_edges, pitch_edges=pitch_edges)
        elif dist == "6d":
            self._require("_dist6d")
            dist = self._dist6d.get()
        elif dist == "rho5d":
            self._require("_distrho5d")
            dist = self._distrho5d.get()
            if exi:
                if ekin_edges is None:
                    ekin_edges  = int(dist.abscissa("ppar").size / 2)
                if pitch_edges is None:
                    pitch_edges = dist.abscissa("pperp").size
                dist = Dist.ppappe2ekinpitch(
                    dist, mass, ekin_edges=ekin_edges, pitch_edges=pitch_edges)
        elif dist == "rho6d":
            self._require("_distrho6d")
            dist = self._distrho6d.get()
        elif dist == "com":
            self._require("_distcom")
            dist = self._distcom.get()
        else:
            raise ValueError("Unknown distribution")

        return dist

    def getdist_moments(self, dist, *moments):
        """Calculate moments of distribution.

        Parameters
        ----------
        dist : str or :class:`DistData`
            Distribution from which moments are calculated.
        *moments : str
            Moments to be calculated.

        Returns
        -------
        out : :class:`DistMoment`
            Distribution object containing moments as ordinates.
        """
        if "rho" in dist.abscissae:
            rhodist = True
            nrho   = dist.abscissa_edges("rho").size
            ntheta = dist.abscissa_edges("theta").size
            nphi   = dist.abscissa_edges("phi").size
            volume, area, r, phi, z = self._root._ascot.input_rhovolume(
                method="prism", nrho=nrho, ntheta=ntheta, nphi=nphi,
                return_area=True, return_coords=True)
            out = DistMoment(
                dist.abscissa_edges("rho"), dist.abscissa_edges("theta"),
                dist.abscissa_edges("phi"), r, phi, z, area, volume, rhodist)
        elif all([x in dist.abscissae for x in ["r", "phi", "z"]]):
            rhodist = False
            r = dist.abscissa_edges("r")
            p = dist.abscissa_edges("phi").to("rad").v
            z = dist.abscissa_edges("z")
            v1, v2 = np.meshgrid( r[1:] - r[:-1], z[1:] - z[:-1])
            area = v1 * v2
            v1, v2, v3 = np.meshgrid(
                p[1:] - p[:-1], r[1:]**2 - r[:-1]**2, z[1:] - z[:-1])
            volume = 0.5 * v1 * v2 * v3
            r, phi, z = np.meshgrid(
                dist.abscissa("r"), dist.abscissa("phi"), dist.abscissa("z"))
            out = DistMoment(
                dist.abscissa_edges("r"), dist.abscissa_edges("phi"),
                dist.abscissa_edges("z"), r, phi, z, area, volume, rhodist)
        else:
            raise ValueError("Distribution has neither rho nor R, phi and z")

        ordinates = {}
        if "density" in moments: 
            Dist.density(dist, out)
        if "chargedensity" in moments:
            Dist.chargedensity(dist, out)
        if "energydensity" in moments: 
            Dist.energydensity(dist, out)
        if "toroidalcurrent" in moments:
            mass = np.mean(self.getstate("mass"))
            Dist.toroidalcurrent(self._root._ascot, mass, dist, out)
        if "parallelcurrent" in moments:
            mass = np.mean(self.getstate("mass"))
            Dist.parallelcurrent(self._root._ascot, mass, dist, out)
        if "pressure" in moments:
            vnorm = self.getstate("vnorm")
            mass = np.mean(self.getstate("mass"))
            Dist.pressure(vnorm, mass, dist, out)
        if "powerdep" in moments:
            mass = np.mean(self.getstate("mass"))
            Dist.powerdep(self._root._ascot, mass, dist, out)
        if "electronpowerdep" in moments:
            mass = np.mean(self.getstate("mass"))
            Dist.electronpowerdep(self._root._ascot, mass, dist, out)
        if "ionpowerdep" in moments:
            mass = np.mean(self.getstate("mass"))
            Dist.ionpowerdep(self._root._ascot, mass, dist, out)
        if "jxBTorque" in moments:
            mass = np.mean(self.getstate("mass"))
            Dist.jxBTorque(self._root._ascot, mass, dist, out)
        if "collTorque" in moments:
            # WIP
            Dist.collTorque(dist, out)
        if "canMomentTorque" in moments:
            # WIP
            Dist.canMomentTorque(dist, out)
        return out

    def plotdist(self, dist, axes=None, cax=None):
        """Plot distribution in 1D or 2D.

        This method assumes that the input distribution has been integrated,
        sliced, and interpolated so that only one or two dimensions have
        a size above one.

        Parameters
        ----------
        dist : :class:`DistData`
            The distribution data object.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        """
        x = None; y = None;
        for key in dist.abscissae:
            val = dist.abscissa_edges(key)
            if val.size > 2:
                if x is None:
                    x = val
                    xlabel = key + " [" + str(x.units) + "]"
                elif y is None:
                    y = val
                    ylabel = key + " [" + str(y.units) + "]"
                else:
                    raise ValueError(
                        "The distribution has more than two dimensions with "
                        + "size greater than one")
        if x is None: raise ValueError("The distribution is zero dimensional")

        ordinate = np.squeeze(dist.distribution())
        if y is None:
            ylabel = "f" + " [" + str(ordinate.units) + "]"
            a5plt.mesh1d(x, ordinate, xlabel=xlabel, ylabel=ylabel, axes=axes)
        else:
            axesequal = x.units == y.units
            clabel = "f" + " [" + str(ordinate.units) + "]"
            a5plt.mesh2d(x, y, ordinate, axesequal=axesequal, xlabel=xlabel,
                         ylabel=ylabel, clabel=clabel, axes=axes, cax=cax)

    def plotdist_moments(self, moment, ordinate, axes=None, cax=None):
        """Plot radial or (R,z) profile of a distribution moment.

        The plotted profile is the average of (theta, phi) or phi depending
        on whether the input is calculated from a rho distribution or not.

        Parameters
        ----------
        moment: class:`DistMoment`
            Moments calculated from the distribution.
        ordinate : str
            Name of the moment to be plotted.
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.
        """
        if moment.rhodist:
            ylabel = ordinate
            ordinate = moment.ordinate(ordinate, toravg=True, polavg=True)
            ylabel += " [" + str(ordinate.units) + "]"
            a5plt.mesh1d(moment.rho, ordinate,
                         xlabel="Normalized poloidal flux",
                         ylabel=ylabel, axes=axes)
        else:
            clabel = ordinate
            ordinate = moment.ordinate(ordinate, toravg=True)
            clabel += " [" + str(ordinate.units) + "]"
            a5plt.mesh2d(moment.r, moment.z, ordinate, axesequal=True,
                         xlabel="R [m]", ylabel="z [m]", clabel=clabel,
                         axes=axes, cax=cax)

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
            xcs = [xcs[i].v for i in idx]
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

    def plotorbit_poincare(self, plane, connlen=True, axes=None, cax=None):
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
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        cax : :obj:`~matplotlib.axes.Axes`, optional
            The color bar axes or otherwise taken from the main axes.

        Raises
        ------
        ValueError
            Raised when the plane is unknown.
        """
        # Set quantities corresponding to the planetype
        pol, tor, rad = self.getorbit_poincareplanes()
        if "pol" in plane:
            plane = int(plane.split()[1]) - 1
            if plane < len(pol) and plane >= 0:
                plane = pol[plane][1]
                x = "r"
                y = "z"
                xlabel = "R [m]"
                ylabel = "z [m]"
                xlim = None # Determined spearately
                ylim = None # Determined separately
                axesequal = True
            else:
                raise ValueError("Unknown plane.")
        elif "tor" in plane:
            plane = int(plane.split()[1]) - 1
            if plane < len(tor) and plane >= 0:
                plane = tor[plane][1]
                x = "rho"
                y = "phimod"
                xlabel = "Normalized poloidal flux"
                ylabel = "Toroidal angle [deg]"
                xlim = [0, 1.1]
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

        if xlim == None:
            # Determine (R,z) limits for poloidal plane by making sure the data
            # fits nicely and then rounding to nearest decimal.
            xlim = [np.floor(np.amin(x)*9) / 10, np.ceil(np.amax(x)*11) / 10]
            ylim = [np.floor(np.amin(y)*11) / 10, np.ceil(np.amax(y)*11) / 10]

        if plotconnlen:
            # Now set confined markers as having negative connection length
            connlen *= -1
            lost1 = self.getstate("ids", state="end", endcond="rhomax wall")

            idx = ~np.in1d(ids, lost1)
            connlen[idx] *= -1
            clabel = "Connection length [" + str(connlen.units) + "]"
        else:
            connlen = None
            clabel = None

        a5plt.poincare(x, y, ids, connlen=connlen, xlim=xlim, ylim=ylim,
                       xlabel=xlabel, ylabel=ylabel, clabel=clabel,
                       axesequal=axesequal, axes=axes, cax=cax)

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
