"""Methods to init, run, and process interactive simulations.

Interactive simulations are simulations that are run via libascot.so.
They are equivalent to normal ascot5_main simulations, except that the
output data is not written to the HDF5 file but the C structs are accessed
directly from Python.
"""
import unyt
import ctypes
import numpy as np

from numpy.ctypeslib import ndpointer

from a5py.physlib import momentum_velocity
from a5py.routines.virtualrun import VirtualRun
from a5py.exceptions import *

from .libascot import _LIBASCOT
if _LIBASCOT:
    from . import ascot2py

class LibSimulate():
    """Mixin class that introduces methods for active simulations.
    """

    def simulation_initoptions(self, **opt):
        """Set simulation options for the interactive simulation.

        There is no need to free options before setting new ones.

        Parameters
        ----------
        **opt :
            Options in same format as accepted by :meth:`Opt.write_hdf5`.
        """
        self._virtualoptions = opt
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
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
        eV2J = unyt.e.base_value
        self._sim.endcond_max_simtime = opt["ENDCOND_MAX_SIMTIME"]
        self._sim.endcond_max_mileage = opt["ENDCOND_MAX_MILEAGE"]
        self._sim.endcond_max_cputime = opt["ENDCOND_MAX_CPUTIME"]
        self._sim.endcond_max_rho     = opt["ENDCOND_MAX_RHO"]
        self._sim.endcond_min_rho     = opt["ENDCOND_MIN_RHO"]
        self._sim.endcond_min_thermal = opt["ENDCOND_MIN_THERMAL"]
        self._sim.endcond_min_ekin    = opt["ENDCOND_MIN_ENERGY"]*eV2J
        self._sim.endcond_max_polorb  = 2*np.pi*opt["ENDCOND_MAX_POLOIDALORBS"]
        self._sim.endcond_max_tororb  = 2*np.pi*opt["ENDCOND_MAX_TOROIDALORBS"]

        diag = self._sim.diag_offload_data
        diag.dist5D_collect    = int(opt["ENABLE_DIST_5D"]) * 0    # Not impl.
        diag.dist6D_collect    = int(opt["ENABLE_DIST_6D"]) * 0    # Not impl.
        diag.distrho5D_collect = int(opt["ENABLE_DIST_RHO5D"]) * 0 # Not impl.
        diag.distrho6D_collect = int(opt["ENABLE_DIST_RHO6D"]) * 0 # Not impl.
        diag.diagtrcof_collect = int(opt["ENABLE_TRANSCOEF"]) * 0  # Not impl.
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
        for i in range(diagorb.ntoroidalplots):
            diagorb.toroidalangles[i] = torangs[i] * np.pi / 180

        diagorb.npoloidalplots = len(polangs)
        for i in range(diagorb.npoloidalplots):
            diagorb.poloidalangles[i] = polangs[i] * np.pi / 180

        diagorb.nradialplots = len(radials)
        for i in range(diagorb.nradialplots):
            diagorb.radialdistances[i] = radials[i]

    def simulation_initinputs(self):
        """Prepare input fields for the interactive simulation.

        Initializes simulation inputs. The inputs used in the simulation are
        those that are active.

        This method differs from :meth:`input_init` in that here the inputs are
        "packed" internally in a single (offload) array as this is what ASCOT5
        does. Inputs cannot be changed while the they are packed. The unpacking
        is done with :meth:`simulation_free`.

        This method must be called before running the simulation.
        """
        self.input_init(bfield=True, efield=True, plasma=True, neutral=True,
                        wall=True, boozer=True, mhd=True, switch=True)
        self._pack()


    def simulation_initmarkers(self, **mrk):
        """Create markers for the interactive simulations.

        Any existing markers are deallocated when new ones are created.

        Parameters
        ----------
        **mrk :
            Marker input (all marker types are supported) as accepted by
            :meth:`Prt.write_hdf5`, :meth:`GC.write_hdf5`, or
            :meth:`FL.write_hdf5`.
        """
        self._virtualmarkers = mrk
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        if self._nmrk.value > 0:
            ascot2py.libascot_deallocate(self._inistate)
            self._nmrk.value = 0

        nmrk = mrk["n"]
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

        #pout = ctypes.POINTER(ascot2py.struct_c__SA_input_particle)()
        ascot2py.prepare_markers(
            ctypes.byref(self._sim), self._mpi_size, self._mpi_rank, nmrk,
            ctypes.byref(pin), ctypes.byref(self._inistate),
            ctypes.byref(self._nmrk), self._bfield_offload_array)

    def simulation_run(self, printsummary=True):
        """Run the interactive simulation using inputs, options and markers that
        were set.

        Parameters
        ----------
        printsummary : bool, optional
            If True, summary of marker endstates is printed after the simulation
            completes.

        Returns
        -------
        run : :class:`VirtualRun`
            An object that acts almost exactly as :class:`RunGroup` except that
            the data is read from C arrays in the memory instead of HDF5 file.

            The run can be used until the output is freed with
            :meth:`simulation_free`. Previous data must be freed before
            rerunning the simulation.

        Raises
        ------
        AscotInitException
            If inputs are not packed, markers are not initialized or previous
            results have not been freed.
        """
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        if not self._offload_ready:
            raise AscotInitException(
                "Initialize inputs before running the simulation")
        if self._nmrk.value == 0:
            raise AscotInitException(
                "Initialize markers before running the simulation")
        if self._diag_occupied:
            raise AscotInitException(
                "Free previous results before running the simulation")

        # Copy inistate to "endstate" and give endstate to ASCOT5 as an input.
        # Otherwise the inistate would get overwritten as the code updates
        # the marker positions in the input array.
        self._endstate = ascot2py.libascot_allocate_particle_states(self._nmrk)
        for j in range(self._nmrk.value):
            for i in range(len(ascot2py.particle_state._fields_)):
                name = self._inistate[j]._fields_[i][0]
                val  = getattr(self._inistate[j], name)
                setattr(self._endstate[j], name, val)

        ascot2py.offload_and_simulate(
            ctypes.byref(self._sim), self._mpi_size, self._mpi_rank,
            self._mpi_root, self._nmrk, self._nmrk, self._endstate,
            ctypes.byref(self._offload_data), self._offload_array,
            self._int_offload_array, ctypes.byref(self._endstate),
            ctypes.byref(self._diag_offload_array))

        self._diag_occupied = True
        if printsummary:
            ascot2py.print_marker_summary(self._endstate, self._nmrk)

        class VirtualInput():
            """Wrapper for marker and options inputs.
            """

            def __init__(self, inp):
                self.inp = inp

            def read(self):
                return self.inp

        return VirtualRun(self, self._nmrk.value,
                          self._inistate, self._endstate,
                          self._sim.diag_offload_data.diagorb.Npnt,
                          self._sim.diag_offload_data.diagorb.record_mode,
                          self._sim.diag_offload_data.diagorb.mode,
                          VirtualInput(self._virtualoptions),
                          VirtualInput(self._virtualmarkers),
                          self._diag_offload_array)

    def simulation_free(self, inputs=False, markers=False, diagnostics=False):
        """Free resources used by the interactive simulation.

        If called without arguments, everything will be freed.

        Parameters
        ----------
        inputs : bool, optional
            If True, inputs are unpacked and freed.
        markers : bool, optional
            If True, markers are freed.
        diagnostics : bool, optional
            If True, diagnostics data is freed.
        """
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        if not inputs and not markers and not diagnostics:
            inputs      = True
            markers     = True
            diagnostics = True

        if inputs and self._offload_ready:
            self._unpack()
        if markers and self._nmrk.value > 0:
            self._nmrk.value = 0
            ascot2py.libascot_deallocate(self._inistate)
            self._virtualmarkers = None
        if diagnostics and self._diag_occupied:
            self._diag_occupied = False
            ascot2py.libascot_deallocate(self._diag_offload_array)
            ascot2py.libascot_deallocate(self._endstate)
