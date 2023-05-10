"""
"""
import unyt
import ctypes
import numpy as np
from a5py.ascotpy import ascotpy2
from a5py.ascotpy.libascot import LibAscot, AscotpyInitException

from numpy.ctypeslib import ndpointer

class LibSimulate(LibAscot):


    def simulation_initoptions(self, opt):
        """
        Set simulation options from the given options dict.
        """

        # Simulation mode options
        self.sim.sim_mode    = opt["SIM_MODE"];
        self.sim.enable_ada  = opt["ENABLE_ADAPTIVE"]
        self.sim.record_mode = opt["RECORD_MODE"]

        # Time step
        self.sim.fix_usrdef_use    = opt["FIXEDSTEP_USE_USERDEFINED"]
        self.sim.fix_usrdef_val    = opt["FIXEDSTEP_USERDEFINED"]
        self.sim.fix_gyrodef_nstep = opt["FIXEDSTEP_GYRODEFINED"]
        self.sim.ada_tol_orbfol    = opt["ADAPTIVE_TOL_ORBIT"]
        self.sim.ada_tol_clmbcol   = opt["ADAPTIVE_TOL_CCOL"]
        self.sim.ada_max_drho      = opt["ADAPTIVE_MAX_DRHO"]
        self.sim.ada_max_dphi      = opt["ADAPTIVE_MAX_DPHI"]

        # Physics
        self.sim.enable_orbfol       = opt["ENABLE_ORBIT_FOLLOWING"]
        self.sim.enable_clmbcol      = opt["ENABLE_COULOMB_COLLISIONS"]
        self.sim.enable_mhd          = opt["ENABLE_MHD"]
        self.sim.disable_gctransform = opt["DISABLE_FIRSTORDER_GCTRANS"]
        self.sim.disable_energyccoll = opt["DISABLE_ENERGY_CCOLL"]
        self.sim.disable_pitchccoll  = opt["DISABLE_PITCH_CCOLL"]
        self.sim.disable_gcdiffccoll = opt["DISABLE_GCDIFF_CCOLL"]

        # Which end conditions are active
        self.sim.endcond_active = 0;
        if opt["ENDCOND_SIMTIMELIM"]:
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_tmax
        if opt["ENDCOND_CPUTIMELIM"]:
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_cpumax
        if opt["ENDCOND_RHOLIM"]:
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_rhomin
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_rhomax
        if opt["ENDCOND_ENERGYLIM"]:
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_emin
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_therm
        if opt["ENDCOND_WALLHIT"]:
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_wall
        if opt["ENDCOND_MAXORBS"]:
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_polmax
            self.sim.endcond_active = self.sim.endcond_active | ascotpy2.endcond_tormax
        self.sim.endcond_torandpol  = 0 + 1 * (opt["ENDCOND_MAXORBS"] == 2)

        # End condition parameters
        self.sim.endcond_max_simtime = opt["ENDCOND_MAX_SIMTIME"]
        self.sim.endcond_max_mileage = opt["ENDCOND_MAX_MILEAGE"]
        self.sim.endcond_max_cputime = opt["ENDCOND_MAX_CPUTIME"]
        self.sim.endcond_max_rho     = opt["ENDCOND_MAX_RHO"]
        self.sim.endcond_min_rho     = opt["ENDCOND_MIN_RHO"]
        self.sim.endcond_min_thermal = opt["ENDCOND_MIN_THERMAL"]
        self.sim.endcond_min_ekin    = opt["ENDCOND_MIN_ENERGY"]*unyt.e.base_value
        self.sim.endcond_max_polorb  = 2*np.pi * opt["ENDCOND_MAX_POLOIDALORBS"]
        self.sim.endcond_max_tororb  = 2*np.pi * opt["ENDCOND_MAX_TOROIDALORBS"]

        diag = self.sim.diag_offload_data
        diag.dist5D_collect    = opt["ENABLE_DIST_5D"] * 0 # Not implemented in output yet
        diag.dist6D_collect    = opt["ENABLE_DIST_6D"] * 0 # Not implemented in output yet
        diag.distrho5D_collect = opt["ENABLE_DIST_RHO5D"] * 0 # Not implemented in output yet
        diag.distrho6D_collect = opt["ENABLE_DIST_RHO6D"] * 0 # Not implemented in output yet
        diag.diagorb_collect   = opt["ENABLE_ORBITWRITE"]
        diag.diagtrcof_collect = opt["ENABLE_TRANSCOEF"] * 0 # Not implemented in output yet

        diagorb = diag.diagorb
        diagorb.mode          = opt["ORBITWRITE_MODE"]
        diagorb.Npnt          = opt["ORBITWRITE_NPOINT"]
        diagorb.writeInterval = opt["ORBITWRITE_INTERVAL"]

        diagorb.record_mode = self.sim.sim_mode
        if self.sim.record_mode and \
           (sim.sim_mode == ascotpy2.simulate_mode_fo or
            sim.sim_mode == ascotpy2.simulate_mode_hybrid):
            diagorb.record_mode = ascotpy2.simulate_mode_gc


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
        Initialize all required inputs and pack them.
        """
        self.init(bfield=True, efield=True, plasma=True, neutral=True,
                  wall=True, boozer=True, mhd=True)
        self.pack()


    def simulation_initmarkers(self, mrktype, mrk):
        """
        Create markers from the given marker dict.
        """
        if self._nmrk.value > 0:
            ascotpy2.libascot_deallocate(self.markers)

        nmrk = mrk["n"]
        pin = ascotpy2.libascot_allocate_input_particles(nmrk)
        prttypes = ascotpy2.input_particle_type__enumvalues

        # particle
        if mrktype == "particle":
            for i in range(nmrk):
                pin[i].type = ascotpy2.input_particle_type_p
                p = pin[i].c__SA_input_particle_0.p

                p.r       = mrk["r"][i]
                p.phi     = mrk["phi"][i]
                p.z       = mrk["z"][i]
                p.p_r     = mrk["p_r"][i]
                p.p_phi   = mrk["p_phi"][i]
                p.p_z     = mrk["p_z"][i]
                p.mass    = mrk["mass"][i]
                p.charge  = mrk["charge"][i]
                p.anum    = mrk["anum"][i]
                p.znum    = mrk["znum"][i]
                p.weight  = mrk["weight"][i]
                p.time    = mrk["time"][i]
                p.mileage = mrk["mileage"][i]
                p.id      = mrk["ids"][i]

        # particle gc
        if mrktype == "gc":
            for i in range(nmrk):
                pin[i].type = ascotpy2.input_particle_type_gc
                p = pin[i].c__SA_input_particle_0.p_gc

                p.r       = mrk["r"][i]
                p.phi     = mrk["phi"][i]
                p.z       = mrk["z"][i]
                p.energy  = mrk["energy"][i]
                p.pitch   = mrk["pitch"][i]
                p.zeta    = mrk["zeta"][i]
                p.mass    = mrk["mass"][i]
                p.charge  = mrk["charge"][i]
                p.anum    = mrk["anum"][i]
                p.znum    = mrk["znum"][i]
                p.weight  = mrk["weight"][i]
                p.time    = mrk["time"][i]
                p.id      = mrk["ids"][i]

        # particle fl
        if mrktype == "ml":
            for i in range(nmrk):
                pin[i].type = ascotpy2.input_particle_type_ml
                p = pin[i].c__SA_input_particle_0.p_ml

                p.r       = mrk["r"][i]
                p.phi     = mrk["phi"][i]
                p.z       = mrk["z"][i]
                p.pitch   = mrk["pitch"][i]
                p.weight  = mrk["weight"][i]
                p.time    = mrk["time"][i]
                p.id      = mrk["ids"][i]

        pout = ctypes.POINTER(ascotpy2.struct_c__SA_input_particle)()
        ascotpy2.prepare_markers(
            ctypes.byref(self.sim), self._mpi_size, self._mpi_rank, nmrk,
            ctypes.byref(pin), ctypes.byref(self.markers),
            ctypes.byref(self._nmrk), self.bfield_offload_array)


    def simulation_run(self, printsummary=True):
        """
        Run simulation.
        """
        if not self._offload_ready:
            raise AscotpyInitException(
                "Simulation can't be run before input is initialized")
        if self._nmrk.value == 0:
            raise AscotpyInitException(
                "Simulation can't be run before markers are initialized")
        if self._diag_occupied:
            raise AscotpyInitException(
                "Simulation can't be run while previous result is not free'd")

        pout = ctypes.POINTER(ascotpy2.struct_c__SA_particle_state)()

        ascotpy2.offload_and_simulate(
            ctypes.byref(self.sim), self._mpi_size, self._mpi_rank,
            self._mpi_root, self._nmrk, self._nmrk, self.markers,
            ctypes.byref(self.offload_data), self.offload_array,
            self.int_offload_array, ctypes.byref(pout),
            ctypes.byref(self.diag_offload_array))

        self.markers = pout
        self._diag_occupied = True
        if printsummary:
            ascotpy2.print_marker_summary(self.markers, self._nmrk)


    def simulation_freeoutput(self, markers=True, diagnostics=True):
        if markers and self._nmrk.value == 0:
            raise AscotpyInitException("No markers exist")
        if diagnostics and not self._diag_occupied:
            raise AscotpyInitException("No result exists")

        if markers:
            ascotpy2.libascot_deallocate(self.marker)
        if diagnostics:
            ascotpy2.libascot_deallocate(self.diag_offload_array)


    def simulation_getstate(self, q):
        if self._nmrk.value == 0:
            raise AscotpyInitException("No markers exist")

        arr = [ getattr(self.markers[i], q) for i in range(self._nmrk.value) ]
        return np.array(arr)


    def simulation_getorbit(self, q, ids=None):
        if not self._diag_occupied:
            raise AscotpyInitException("No result exists")

        if self.sim.diag_offload_data.diagorb_collect == 0:
            raise AscotpyInitException("No orbit diagnostics")

        arraylength = self.sim.diag_offload_data.diagorb.Nmrk \
            * self.sim.diag_offload_data.diagorb.Npnt

        mrk = ids
        ids = np.array(self.diag_offload_array[arraylength*0:arraylength*1])
        mil = np.array(self.diag_offload_array[arraylength*1:arraylength*2])
        idx = np.lexsort((mil, ids))

        if q == "id":
            q = self.diag_offload_array[arraylength*0:arraylength*1]
        elif q == "mileage":
            q = self.diag_offload_array[arraylength*1:arraylength*2]
        elif q == "r":
            q = self.diag_offload_array[arraylength*2:arraylength*3]
        elif q == "phi":
            q = self.diag_offload_array[arraylength*3:arraylength*4]
        elif q == "z":
            q = self.diag_offload_array[arraylength*4:arraylength*5]

        if ids is None:
            valid = ids[idx] > 0
        else:
            valid = ids == mrk
        return (np.array(q)[idx])[valid]
