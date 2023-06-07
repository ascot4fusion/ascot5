"""
Methods to init, run, and process interactive simulations.

Interactive simulations are simulations that are run via libascot.so.
They are equivalent to normal ascot5_main simulations, except that the
output data is not written to the HDF5 file but the C structs are accessed
directly from Python.
"""
import unyt
import ctypes
import numpy as np
from a5py.ascotpy import ascot2py

from numpy.ctypeslib import ndpointer

from a5py.physlib.gamma import momentum_velocity

class LibSimulate():

    def simulation_initoptions(self, opt):
        """
        Set simulation options for the interactive simulation.

        Args:
            opt : dict
                Options dictionary with the desired options.
        """

        # Simulation mode options
        self._sim.sim_mode    = int(opt["SIM_MODE"]);
        self._sim.enable_ada  = int(opt["ENABLE_ADAPTIVE"])
        self._sim.record_mode = int(opt["RECORD_MODE"])

        # Time step
        self._sim.fix_usrdef_use    = int(opt["FIXEDSTEP_USE_USERDEFINED"])
        self._sim.fix_usrdef_val    = opt["FIXEDSTEP_USERDEFINED"]
        self._sim.fix_gyrodef_nstep = int(opt["FIXEDSTEP_GYRODEFINED"])
        self._sim.ada_tol_orbfol    = opt["ADAPTIVE_TOL_ORBIT"]
        self._sim.ada_tol_clmbcol   = opt["ADAPTIVE_TOL_CCOL"]
        self._sim.ada_max_drho      = opt["ADAPTIVE_MAX_DRHO"]
        self._sim.ada_max_dphi      = opt["ADAPTIVE_MAX_DPHI"]

        # Physics
        self._sim.enable_orbfol       = int(opt["ENABLE_ORBIT_FOLLOWING"])
        self._sim.enable_clmbcol      = int(opt["ENABLE_COULOMB_COLLISIONS"])
        self._sim.enable_mhd          = int(opt["ENABLE_MHD"])
        self._sim.disable_gctransform = int(opt["DISABLE_FIRSTORDER_GCTRANS"])
        self._sim.disable_energyccoll = int(opt["DISABLE_ENERGY_CCOLL"])
        self._sim.disable_pitchccoll  = int(opt["DISABLE_PITCH_CCOLL"])
        self._sim.disable_gcdiffccoll = int(opt["DISABLE_GCDIFF_CCOLL"])

        # Which end conditions are active
        self._sim.endcond_active = 0;
        if opt["ENDCOND_SIMTIMELIM"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_tmax
        if opt["ENDCOND_CPUTIMELIM"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_cpumax
        if opt["ENDCOND_RHOLIM"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_rhomin
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_rhomax
        if opt["ENDCOND_ENERGYLIM"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_emin
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_therm
        if opt["ENDCOND_WALLHIT"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_wall
        if opt["ENDCOND_MAXORBS"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_polmax
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_tormax
        self._sim.endcond_torandpol  = 0 + 1 * int(opt["ENDCOND_MAXORBS"] == 2)

        # End condition parameters
        self._sim.endcond_max_simtime = opt["ENDCOND_MAX_SIMTIME"]
        self._sim.endcond_max_mileage = opt["ENDCOND_MAX_MILEAGE"]
        self._sim.endcond_max_cputime = opt["ENDCOND_MAX_CPUTIME"]
        self._sim.endcond_max_rho     = opt["ENDCOND_MAX_RHO"]
        self._sim.endcond_min_rho     = opt["ENDCOND_MIN_RHO"]
        self._sim.endcond_min_thermal = opt["ENDCOND_MIN_THERMAL"]
        self._sim.endcond_min_ekin    = opt["ENDCOND_MIN_ENERGY"]*unyt.e.base_value
        self._sim.endcond_max_polorb  = 2*np.pi * opt["ENDCOND_MAX_POLOIDALORBS"]
        self._sim.endcond_max_tororb  = 2*np.pi * opt["ENDCOND_MAX_TOROIDALORBS"]

        diag = self._sim.diag_offload_data
        diag.dist5D_collect    = int(opt["ENABLE_DIST_5D"]) * 0    # Not implemented
        diag.dist6D_collect    = int(opt["ENABLE_DIST_6D"]) * 0    # Not implemented
        diag.distrho5D_collect = int(opt["ENABLE_DIST_RHO5D"]) * 0 # Not implemented
        diag.distrho6D_collect = int(opt["ENABLE_DIST_RHO6D"]) * 0 # Not implemented
        diag.diagtrcof_collect = int(opt["ENABLE_TRANSCOEF"]) * 0  # Not implemented
        diag.diagorb_collect   = int(opt["ENABLE_ORBITWRITE"])

        diagorb = diag.diagorb
        diagorb.mode          = int(opt["ORBITWRITE_MODE"])
        diagorb.Npnt          = int(opt["ORBITWRITE_NPOINT"])
        diagorb.writeInterval = opt["ORBITWRITE_INTERVAL"]

        diagorb.record_mode = self._sim.sim_mode
        if self._sim.record_mode and \
           (self._sim.sim_mode == ascot2py.simulate_mode_fo or
            self._sim.sim_mode == ascot2py.simulate_mode_hybrid):
            diagorb.record_mode = ascot2py.simulate_mode_gc


        torangs = opt["ORBITWRITE_TOROIDALANGLES"]
        torangs = torangs if isinstance(torangs, list) else [torangs]
        polangs = opt["ORBITWRITE_POLOIDALANGLES"]
        polangs = polangs if isinstance(polangs, list) else [polangs]
        radials = opt["ORBITWRITE_RADIALDISTANCES"]
        radials = radials if isinstance(radials, list) else [radials]

        diagorb.ntoroidalplots = len(torangs)
        if torangs: diagorb.ntoroidalplots = 0
        for i in range(diagorb.ntoroidalplots):
            diagorb.toroidalangles[i] = torangs[i] * np.pi / 180

        diagorb.npoloidalplots = len(polangs)
        if polangs: diagorb.npoloidalplots = 0
        for i in range(diagorb.npoloidalplots):
            diagorb.poloidalangles[i] = polangs[i] * np.pi / 180

        diagorb.nradialplots = len(radials)
        if radials: diagorb.nradialplots = 0
        for i in range(diagorb.nradialplots):
            diagorb.radialdistances[i] = radials[i]


    def simulation_initinput(self):
        """
        Prepare input fields for the interactive simulation.

        All necessary inputs are initialized (if they haven't been already) and
        they are packed. Packing means that the data arrays are stored as
        a single array which is required in order to offload the data in ASCOT5.
        In practice this means that the initialized inputs cannot be changed
        before the array is unpacked.

        This method must be called before running the interactive simulation as
        the input array must be packed for that. This method does not concern
        options and markers, so you can modify those even after the input has
        been packed.
        """
        self.input_init(bfield=True, efield=True, plasma=True, neutral=True,
                        wall=True, boozer=True, mhd=True, switch=True)
        self._pack()


    def simulation_initmarkers(self, mrk):
        """
        Create markers for the interactive simulations.

        Any existing markers are deallocated when new ones are created.
        Note that you can use getstate also after creating the markers,
        but before running the simulation. In that case the result
        correspond to marker inistate.

        Args:
            mrk : dict
                Marker dictionary which can be of any marker type.
        """
        if self._nmrk.value > 0:
            ascot2py.libascot_deallocate(self._markers)
            self._nmrk.value = 0

        nmrk = int(mrk["n"])
        pin = ascot2py.libascot_allocate_input_particles(nmrk)
        prttypes = ascot2py.input_particle_type__enumvalues

        # particle
        if "v_r" in mrk:
            for i in range(nmrk):
                pin[i].type = ascot2py.input_particle_type_p
                p = pin[i].c__SA_input_particle_0.p

                vvec = np.array([mrk["v_r"][i], mrk["v_phi"][i], mrk["v_z"][i]])
                pvec = momentum_velocity(mrk["mass"][i], vvec)

                p.r       = mrk["r"][i]
                p.phi     = mrk["phi"][i] * np.pi / 180
                p.z       = mrk["z"][i]
                p.p_r     = pvec[0]
                p.p_phi   = pvec[1]
                p.p_z     = pvec[2]
                p.mass    = mrk["mass"][i]
                p.charge  = mrk["charge"][i]
                p.anum    = mrk["anum"][i]
                p.znum    = mrk["znum"][i]
                p.weight  = mrk["weight"][i]
                p.time    = mrk["time"][i]
                p.mileage = mrk["mileage"][i]
                p.id      = mrk["ids"][i]

        # particle gc
        elif "energy" in mrk:
            for i in range(nmrk):
                pin[i].type = ascot2py.input_particle_type_gc
                p = pin[i].c__SA_input_particle_0.p_gc

                p.r       = mrk["r"][i]
                p.phi     = mrk["phi"][i] * np.pi / 180
                p.z       = mrk["z"][i]
                p.energy  = mrk["energy"][i] * unyt.elementary_charge.value
                p.pitch   = mrk["pitch"][i]
                p.zeta    = mrk["zeta"][i]
                p.mass    = mrk["mass"][i] * unyt.atomic_mass_unit.value
                p.charge  = mrk["charge"][i] * unyt.elementary_charge.value
                p.anum    = mrk["anum"][i]
                p.znum    = mrk["znum"][i]
                p.weight  = mrk["weight"][i]
                p.time    = mrk["time"][i]
                p.id      = mrk["ids"][i]

        # particle fl
        else:
            for i in range(nmrk):
                pin[i].type = ascot2py.input_particle_type_ml
                p = pin[i].c__SA_input_particle_0.p_ml

                p.r       = mrk["r"][i]
                p.phi     = mrk["phi"][i] * np.pi / 180
                p.z       = mrk["z"][i]
                p.pitch   = mrk["pitch"][i]
                p.weight  = mrk["weight"][i]
                p.time    = mrk["time"][i]
                p.id      = mrk["ids"][i]

        pout = ctypes.POINTER(ascot2py.struct_c__SA_input_particle)()
        ascot2py.prepare_markers(
            ctypes.byref(self._sim), self._mpi_size, self._mpi_rank, nmrk,
            ctypes.byref(pin), ctypes.byref(self._markers),
            ctypes.byref(self._nmrk), self._bfield_offload_array)


    def simulation_run(self, printsummary=True):
        """
        Run the interactive simulation.

        This method runs the interactive simulation using the options set with
        simulation_initoptions, markers created with simulation_initmarkers and
        input fields that were initialized in simulation_initinputs.

        After running the simulation, the marker endstates can be accessed with
        simulation_getstate and orbits (if they were enabled in options) with
        simulation_getorbit. Output must be deallocated with simulation_free
        before rerunning the simulation

        Args:
            printsummary : bool, optional
                If True, summary of marker endstates is printed after
                the simulation completes.
        """
        if not self._offload_ready:
            raise AscotInitException(
                "Simulation can't be run before input is initialized")
        if self._nmrk.value == 0:
            raise AscotInitException(
                "Simulation can't be run before markers are initialized")
        if self._diag_occupied:
            raise AscotInitException(
                "Simulation can't be run while previous result is not free'd")

        pout = ctypes.POINTER(ascot2py.struct_c__SA_particle_state)()

        ascot2py.offload_and_simulate(
            ctypes.byref(self._sim), self._mpi_size, self._mpi_rank,
            self._mpi_root, self._nmrk, self._nmrk, self._markers,
            ctypes.byref(self._offload_data), self._offload_array,
            self._int_offload_array, ctypes.byref(pout),
            ctypes.byref(self._diag_offload_array))

        self._markers = pout
        self._diag_occupied = True
        if printsummary:
            ascot2py.print_marker_summary(self._markers, self._nmrk)


    def simulation_free(self, inputs=False, markers=False, diagnostics=False):
        """Free resources used by the interactive simulation.

        Args:
          inputs : bool, optional
            If True, inputs are unpacked (but not free'd).
          markers : bool, optional
            If True, markers are free'd.
          diagnostics : bool, optional
            If True, diagnostics data is free'd.
        """
        if inputs:
            self._unpack()

        if markers:
            if self._nmrk.value == 0:
                raise AscotpyInitException("No markers exist")

            self._nmrk.value = 0
            ascot2py.libascot_deallocate(self._markers)
        if diagnostics:
            if not self._diag_occupied:
                raise AscotpyInitException("No result exists")

            self._diag_occupied = False
            ascot2py.libascot_deallocate(self._diag_offload_array)


    def simulation_getstate(self, quantity):
        """
        Access marker state data in interactive simulations.

        This method returns quantities from marker endstate (or inistate if
        called before running the simulation).

        Args:
            quantity : str
                Name of the quantity.
        Returns:
          np.array
            Array with the quantity for all markers.
        """
        if self._nmrk.value == 0:
            raise AscotpyInitException("No markers exist")

        arr = [ getattr(self._markers[i], quantity)
                for i in range(self._nmrk.value) ]
        return np.array(arr)


    def simulation_getorbit(self, q, ids=None):
        """
        Access orbit data in interactive simulations.

        This method can be called once simulation has been run but before
        the resources are free'd. The simulation must also collect the orbit
        data for it to be accessible.

        Args:
          quantity : str
            Name of the quantity.
          ids : array_like, int, optional
            ID(s) of the marker(s) whose orbit data is returned. If None, then
            data for all markers is returned.
        Returns:
          np.array
            Requested quantity sorted by marker ID and mileage.
        """
        if not self._diag_occupied:
            raise AscotInitException("No result exists")

        if self._sim.diag_offload_data.diagorb_collect == 0:
            raise AscotInitException("No orbit diagnostics")

        arraylength = self._sim.diag_offload_data.diagorb.Nmrk \
            * self._sim.diag_offload_data.diagorb.Npnt

        mrk = ids
        ids = np.array(self._diag_offload_array[arraylength*0:arraylength*1])
        mil = np.array(self._diag_offload_array[arraylength*1:arraylength*2])
        idx = np.lexsort((mil, ids))

        if q == "id":
            q = self._diag_offload_array[arraylength*0:arraylength*1]
        elif q == "mileage":
            q = self._diag_offload_array[arraylength*1:arraylength*2]
        elif q == "r":
            q = self._diag_offload_array[arraylength*2:arraylength*3]
        elif q == "phi":
            q = self._diag_offload_array[arraylength*3:arraylength*4]
        elif q == "z":
            q = self._diag_offload_array[arraylength*4:arraylength*5]

        if ids is None:
            valid = ids[idx] > 0
        else:
            valid = ids[idx] == mrk
        return (np.array(q)[idx])[valid]
