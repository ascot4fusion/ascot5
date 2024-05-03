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
from a5py.ascot5io import Marker, State, Orbits, Dist
from a5py.ascot5io.dist import DistMoment
import a5py.physlib as physlib
from a5py.routines.distmixin import DistMixin

class RunMixin(DistMixin):
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
                    ec = ec.replace(" and ", " ")
                    accept  = State.endcond_check(uec, ec)
                    mask[i] = mask[i] or accept

            idx = mask[uidx]

        if ids is not None:
            idx = np.logical_and(idx, np.in1d(self._inistate.get("ids"), ids))

        for i in range(len(data)):
            data[i] = data[i][idx]
        if "mu" in qnt:
            data[qnt.index("mu")].convert_to_units("eV/T")
        if "psi" in qnt:
            data[qnt.index("psi")].convert_to_units("Wb")
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
        if "psi" in qnt:
            data[qnt.index("psi")].convert_to_units("Wb")
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

    def getstate_losssummary(self):
        """Return a summary of lost markers.
        """
        self._require("_endstate")

        wmrks, emrks = self.getstate("weight", "ekin", state="ini")
        wloss, eloss = self.getstate("weight", "ekin", state="end",
                                     endcond="wall")

        markers_lost   = wloss.size
        markers_frac   = wloss.size / wmrks.size
        particles_lost = np.sum(wloss)
        particles_frac = np.sum(wloss) / np.sum(wmrks)
        power_lost     = np.sum(wloss * eloss.to("J"))
        power_frac     = np.sum(wloss * eloss.to("J")) / np.sum(wmrks * emrks)

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
            for i, q in enumerate(qnt):
                mrk[q] = state[i]
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
            for i, q in enumerate(qnt):
                mrk[q] = state[i]
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
            if not isinstance(val, list): val = [val]
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

    def getwall_loads(self, weights=True, p_ids=None):
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

        p_ids : array_like, optional
            Calculate wall loads only for the particles with the given indices.

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
                List of triangle indecies for which the 3D mesh is made.

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
        """

        wallmesh = pv.PolyData(
            *self.wall.noderepresentation(w_indices=w_indices))
        ids, area, eload, pload, iangle = self.getwall_loads(p_ids=p_ids)
        # get the marker load at each tile
        _, _, _, mload, _ = self.getwall_loads(weights=False, p_ids=p_ids)
        ids = ids - 1 # Convert IDs to indices

        #here possible area filter to ids, TBI
        #filt_area = area > threshold
        n_tri = self.wall.getNumberOfElements()

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

        ntriangle = wallmesh.n_faces
        wallmesh.cell_data["pload"]  = pload_tot[wall_f]
        wallmesh.cell_data["eload"]  = eload_tot[wall_f]
        wallmesh.cell_data["mload"]  = mload_tot[wall_f]
        wallmesh.cell_data["iangle"] = iangle_tot[wall_f]
        return wallmesh

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

    def getdist_moments(self, dist, *moments, volmethod="prism"):
        """Calculate moments of distribution.

        Parameters
        ----------
        dist : str or :class:`DistData`
            Distribution from which moments are calculated.
        *moments : str
            Moments to be calculated.
        volmethod : {"mc", "prism"}, optional
            Method used to calculate the volume.

            It is a good idea to verify that same results are obtained with both
            methods.

        Returns
        -------
        out : :class:`DistMoment`
            Distribution object containing moments as ordinates.
        """
        # Initialize the moment object and evaluate the volume.
        if all([x in dist.abscissae for x in ["rho", "theta", "phi"]]):
            rhodist = True
            nrho   = dist.abscissa_edges("rho").size
            ntheta = dist.abscissa_edges("theta").size
            nphi   = dist.abscissa_edges("phi").size
            volume, area, r, phi, z = self._root._ascot.input_rhovolume(
                method=volmethod, tol=1e-1, nrho=nrho, ntheta=ntheta, nphi=nphi,
                return_area=True, return_coords=True)
            volume[volume == 0] = 1e-8 # To avoid division by zero
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
            raise ValueError(
                "Distribution has neither (rho, theta, phi) nor (R, phi, z)")

        ordinates = {}
        mass = np.mean(self.getstate("mass"))
        if "density" in moments:
            Dist.density(dist, out)
        if "chargedensity" in moments:
            Dist.chargedensity(dist, out)
        if "energydensity" in moments:
            Dist.energydensity(mass, dist, out)
        if "toroidalcurrent" in moments:
            Dist.toroidalcurrent(self._root._ascot, mass, dist, out)
        if "parallelcurrent" in moments:
            Dist.parallelcurrent(self._root._ascot, mass, dist, out)
        if "pressure" in moments:
            Dist.pressure(mass, dist, out)
        if "powerdep" in moments:
            Dist.powerdep(self._root._ascot, mass, dist, out)
        if "electronpowerdep" in moments:
            Dist.electronpowerdep(self._root._ascot, mass, dist, out)
        if "ionpowerdep" in moments:
            Dist.ionpowerdep(self._root._ascot, mass, dist, out)
        #if "jxbtorque" in moments:
        #    Dist.jxbtorque(self._root._ascot, mass, dist, out)
        #if "colltorque" in moments:
        #    Dist.collTorque(self._root._ascot, mass, dist, out)
        #if "canmomtorque" in moments:
        #    Dist.canMomentTorque(dist, out)
        return out

    def getdist_list(self, show=True):
        """List all available distributions and moments.

        Parameters
        ----------
        show : bool, optional
            Print the output.

        Returns
        -------
        dists : list [str]
            List of available distributions.
        moms : list [(str, str)]
            List of distribution moments that can be evaluated and their
            meaning.
        """
        dists = []
        if hasattr(self, "_dist5d"):    dists.append("5d")
        if hasattr(self, "_distrho5d"): dists.append("rho5d")
        if hasattr(self, "_dist6d"):    dists.append("6d")
        if hasattr(self, "_distrho6d"): dists.append("rho6d")
        if hasattr(self, "_distcom"):   dists.append("com")

        moms = []
        if "5d" not in dists and "rho5d" not in dists:
            if show:
                print("Available distributions:")
                for d in dists:
                    print(d)
                print("\nMoments available only for 5d/rho5d dists.")
            return dists, moms

        moms = [
            ("density", "Particle number density"),
            ("chargedensity", "Charge density"),
            ("energydensity", "Energy density"),
            ("toroidalcurrent", "Toroidal current"),
            ("parallelcurrent", "Parallel current"),
            ("pressure", "Pressure"),
            ("powerdep", "Total deposited power"),
            ("ionpowerdep", "Power deposited to ions"),
            ("electronpowerdep", "Power deposited to electrons"),
            #("jxbtorque", "j_rad x B_pol torque"),
            #("colltorque", "Torque from collisions"),
            #("canmomtorque", "Torque from change in can. tor. ang. momentum"),
        ]
        if show:
            print("Available distributions:")
            for d in dists:
                print(d)
            print("Available Moments:")
            for d in moms:
                print(d[0] + " : " + d[1])
        return dists, moms

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

    @a5plt.openfigureifnoaxes(projection=None)
    def plotwall_convergence(self, qnt, nmin=1000, nsample=10, axes=None):
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

            if qnt == 'lostpower':
                val[i] = np.sum((lost*weight*ekin.to('J'))[idx]) \
                    * wtotal / np.sum(weight[idx])
                axes.set_ylabel('Lost power [W]')
            elif qnt == 'peakload':
                ids0  = ids[idx]
                lost0 = lost[idx]
                _, area, loads, _, _ = self.getwall_loads(weights=True, p_ids=ids0[lost0])
                val[i] = np.amax(loads/area) * wtotal / np.sum(weight[idx])
                axes.set_ylabel(r'Peak load [W/m$^2$]')

        axes.set_xscale('log')
        axes.set_xlabel('Number of markers')
        axes.plot(nsubset, val)

    def plotwall_loadvsarea(self, axes=None):
        """Plot histogram showing area affected by at least a given load.

        Parameters
        ----------
        axes : :obj:`~matplotlib.axes.Axes`, optional
            The axes where figure is plotted or otherwise new figure is created.
        """
        _, area, eload, _, _ = self.getwall_loads()
        a5plt.loadvsarea(area, eload/area, axes=axes)

    @a5plt.openfigureifnoaxes(projection=None)
    def plotwall_torpol(self, qnt='eload', getaxis=None,log=True, clim=None,
                        cmap=None, axes=None, cax=None):
        """Plot the toroidal-poloidal projection of the wall loads.

        Note that the resulting doesn't present the areas accurately and even
        the locations are approximate. The purpose of the plot is to roughly
        illustrate how the losses are distributed.

        Parameters
        ----------
        qnt : {'eload', 'pload', 'iangle'}, optional
            Quantity to plot.
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
        wetted, area, edepo, pdepo, iangle = self.getwall_loads()
        d = self.wall.read()
        nelement = wetted.size
        if qnt == 'eload':
            color = edepo/area
            clabel = r"Wall load [W/m$^2$]"
            if cmap is None: cmap = 'Reds'
        if qnt == 'pload':
            color = pdepo/area
            clabel = r"Particle flux [prt s$^{-1}$m$^{-2}$]"
            if cmap is None: cmap = 'Reds'
        elif qnt == 'iangle':
            color = iangle
            clabel = r"Angle of incidence [deg]"
            if cmap is None: cmap = 'viridis'

        x1x2x3   = d["x1x2x3"][wetted-1]
        y1y2y3   = d["y1y2y3"][wetted-1]
        z1z2z3   = d["z1z2z3"][wetted-1]

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
