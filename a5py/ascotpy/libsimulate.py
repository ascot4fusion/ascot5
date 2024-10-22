"""Methods to init, run, and process interactive simulations.

Interactive simulations are simulations that are run via libascot.so.
They are equivalent to normal ascot5_main simulations, except that the
output data is not written to the HDF5 file but the C structs are accessed
directly from Python.
"""
import unyt
import ctypes
import numpy as np
import wurlitzer # For muting libascot.so

from numpy.ctypeslib import ndpointer

from a5py import physlib
from a5py.routines.virtualrun import VirtualRun, VirtualBBNBIRun
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

            If not given, the options are read from the HDF5 file.
        """
        self._virtualoptions = opt
        if not _LIBASCOT:
            # Raise exception if libascot.so is not found
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        if not opt:
            # Read options from HDF5 file
            opt = self.data.options.active.read()

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
        self._sim.enable_atomic       = int(opt["ENABLE_ATOMIC"])
        self._sim.disable_gctransform = int(opt["DISABLE_FIRSTORDER_GCTRANS"])
        self._sim.disable_energyccoll = int(opt["DISABLE_ENERGY_CCOLL"])
        self._sim.disable_pitchccoll  = int(opt["DISABLE_PITCH_CCOLL"])
        self._sim.disable_gcdiffccoll = int(opt["DISABLE_GCDIFF_CCOLL"])
        self._sim.reverse_time        = int(opt["REVERSE_TIME"])

        # Which end conditions are active
        self._sim.endcond_active = 0;
        if opt["ENDCOND_SIMTIMELIM"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_tlim
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
        if opt["ENDCOND_NEUTRALIZED"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_neutrz
        if opt["ENDCOND_IONIZED"]:
            self._sim.endcond_active = \
                self._sim.endcond_active | ascot2py.endcond_ioniz

        # End condition parameters
        eV2J = unyt.e.base_value
        self._sim.endcond_lim_simtime = opt["ENDCOND_LIM_SIMTIME"]
        self._sim.endcond_max_mileage = opt["ENDCOND_MAX_MILEAGE"]
        self._sim.endcond_max_cputime = opt["ENDCOND_MAX_CPUTIME"]
        self._sim.endcond_max_rho     = opt["ENDCOND_MAX_RHO"]
        self._sim.endcond_min_rho     = opt["ENDCOND_MIN_RHO"]
        self._sim.endcond_min_thermal = opt["ENDCOND_MIN_THERMAL"]
        self._sim.endcond_min_ekin    = opt["ENDCOND_MIN_ENERGY"]*eV2J
        self._sim.endcond_max_polorb  = 2*np.pi*opt["ENDCOND_MAX_POLOIDALORBS"]
        self._sim.endcond_max_tororb  = 2*np.pi*opt["ENDCOND_MAX_TOROIDALORBS"]

        # Setting options
        diag = self._sim.diag_data
        diag.dist5D_collect    = int(opt["ENABLE_DIST_5D"])
        diag.dist6D_collect    = int(opt["ENABLE_DIST_6D"]) * 0    # Not impl.
        diag.distrho5D_collect = int(opt["ENABLE_DIST_RHO5D"])*0
        diag.distrho6D_collect = int(opt["ENABLE_DIST_RHO6D"]) * 0 # Not impl.
        diag.distCOM_collect   = int(opt["ENABLE_DIST_COM"]) * 0   # Not impl.
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

        dist = diag.dist5D
        dist.min_r     = opt["DIST_MIN_R"]
        dist.max_r     = opt["DIST_MAX_R"]
        dist.n_r       = opt["DIST_NBIN_R"]
        dist.min_phi   = opt["DIST_MIN_PHI"] * np.pi / 180
        dist.max_phi   = opt["DIST_MAX_PHI"] * np.pi / 180
        dist.n_phi     = opt["DIST_NBIN_PHI"]
        dist.min_z     = opt["DIST_MIN_Z"]
        dist.max_z     = opt["DIST_MAX_Z"]
        dist.n_z       = opt["DIST_NBIN_Z"]
        dist.min_ppara = opt["DIST_MIN_PPA"]
        dist.max_ppara = opt["DIST_MAX_PPA"]
        dist.n_ppara   = opt["DIST_NBIN_PPA"]
        dist.min_pperp = opt["DIST_MIN_PPE"]
        dist.max_pperp = opt["DIST_MAX_PPE"]
        dist.n_pperp   = opt["DIST_NBIN_PPE"]
        dist.min_time  = opt["DIST_MIN_TIME"]
        dist.max_time  = opt["DIST_MAX_TIME"]
        dist.n_time    = opt["DIST_NBIN_TIME"]
        dist.min_q     = opt["DIST_MIN_CHARGE"]
        dist.max_q     = opt["DIST_MAX_CHARGE"]
        dist.n_q       = opt["DIST_NBIN_CHARGE"]

    def simulation_initinputs(self, bfield=True, efield=True, plasma=True,
                              neutral=True, wall=True, boozer=True, mhd=True,
                              asigma=True, switch=True):
        """Prepare input fields for the interactive simulation.

        Initializes simulation inputs. The inputs used in the simulation are
        those that are active.

        This method differs from :meth:`input_init` in that here the inputs are
        "packed" internally in a single (offload) array as this is what ASCOT5
        does. Inputs cannot be changed while the they are packed. The unpacking
        is done with :meth:`simulation_free`.

        This method must be called before running the simulation.
        """
        self.input_init(
            bfield=bfield, efield=efield, plasma=plasma, neutral=neutral,
            wall=wall, boozer=boozer, mhd=mhd, asigma=asigma, switch=switch)

    def simulation_initbbnbi(
            self,
            bfield=True,
            plasma=True,
            neutral=True,
            wall=True,
            asigma=True,
            nbi=True,
            switch=True):
        """Initialize inputs for BBNBI simulation.

        This method must be called before running BBNBI simulation.
        """
        self.input_init(
            bfield=bfield, plasma=plasma, neutral=neutral,
            wall=wall, asigma=asigma, switch=switch)
        self._init(self.data, nbi=getattr(self.data, "nbi").active.get_qid())

    def simulation_initmarkers(self, **mrk):
        """Create markers for the interactive simulations.

        Any existing markers are deallocated when new ones are created.

        Parameters
        ----------
        **mrk :
            Marker input (all marker types are supported) as accepted by
            :meth:`Prt.write_hdf5`, :meth:`GC.write_hdf5`, or
            :meth:`FL.write_hdf5`.

            If not given, markers are read from the HDF5 file.
        """
        if not _LIBASCOT:
            # If libascot.so is not found, raise exception
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        if self._nmrk.value > 0:
            # Deallocate previous markers
            ascot2py.libascot_deallocate(self._inistate)
            self._nmrk.value = 0
            self._virtualmarkers = None
        if not mrk:
            # Read markers from HDF5 file
            mrk = self.data.marker.active.read()

        self._virtualmarkers = mrk
        nmrk = mrk["n"]
        pin = ascot2py.libascot_allocate_input_particles(nmrk)
        prttypes = ascot2py.input_particle_type__enumvalues

        # particle
        if "vr" in mrk:
            @physlib.parseunits(r="m", phi="deg", z="m", time="s", mass="amu",
                                charge="e", vr="m/s", vphi="m/s", vz="m/s")
            def parse(r, phi, z, time, mass, charge, vr, vphi, vz, anum, znum,
                      weight, ids):
                return r, phi.to("rad"), z, time, mass.to("kg"), \
                    charge.to("C"), vr, vphi, vz, anum, znum, weight, ids

            r, phi, z, t, m, q, vr, vphi, vz, anum, znum, w, ids = parse(
                mrk["r"], mrk["phi"], mrk["z"], mrk["time"], mrk["mass"],
                mrk["charge"], mrk["vr"], mrk["vphi"], mrk["vz"],
                mrk["anum"], mrk["znum"], mrk["weight"], mrk["ids"])

            for i in range(nmrk):
                pin[i].type = ascot2py.input_particle_type_p
                p = pin[i].p

                vvec = np.array([vr[i], vphi[i], vz[i]])*unyt.m/unyt.s
                pvec = physlib.momentum_velocity(mrk["mass"][i], vvec)

                p.r       = r[i]
                p.phi     = phi[i]
                p.z       = z[i]
                p.p_r     = pvec[0]
                p.p_phi   = pvec[1]
                p.p_z     = pvec[2]
                p.mass    = m[i]
                p.charge  = q[i]
                p.anum    = anum[i]
                p.znum    = znum[i]
                p.weight  = w[i]
                p.time    = t[i]
                p.id      = ids[i]

        # particle gc
        elif "energy" in mrk:
            @physlib.parseunits(r="m", phi="deg", z="m", time="s", mass="amu",
                                charge="e", energy="eV", zeta="rad")
            def parse(r, phi, z, time, mass, charge, energy, pitch, zeta,
                      anum, znum, weight, ids):
                return r, phi.to("rad"), z, time, mass.to("kg"), \
                    charge.to("C"), energy.to("J"), pitch, zeta, anum, znum, \
                    weight, ids

            r, phi, z, t, m, q, energy, pitch, zeta, anum, znum, w, ids = parse(
                mrk["r"], mrk["phi"], mrk["z"], mrk["time"], mrk["mass"],
                mrk["charge"], mrk["energy"], mrk["pitch"], mrk["zeta"],
                mrk["anum"], mrk["znum"], mrk["weight"], mrk["ids"])

            for i in range(nmrk):
                pin[i].type = ascot2py.input_particle_type_gc
                p = pin[i].p_gc

                p.r       = r[i]
                p.phi     = phi[i]
                p.z       = z[i]
                p.energy  = energy[i]
                p.pitch   = pitch[i]
                p.zeta    = zeta[i]
                p.mass    = m[i]
                p.charge  = q[i]
                p.anum    = anum[i]
                p.znum    = znum[i]
                p.weight  = w[i]
                p.time    = t[i]
                p.id      = ids[i]

        # particle fl
        else:
            @physlib.parseunits(r="m", phi="deg", z="m", time="s",
                                energy="eV", zeta="rad")
            def parse(r, phi, z, time, pitch, weight, ids):
                return r, phi.to("rad"), z, time, pitch, weight, ids

            r, phi, z, t, pitch, w, ids = parse(
                mrk["r"], mrk["phi"], mrk["z"], mrk["time"], mrk["pitch"],
                mrk["weight"], mrk["ids"])

            for i in range(nmrk):
                pin[i].type = ascot2py.input_particle_type_ml
                p = pin[i].p_ml

                p.r       = r[i]
                p.phi     = phi[i]
                p.z       = z[i]
                p.pitch   = pitch[i]
                p.weight  = w[i]
                p.time    = t[i]
                p.id      = ids[i]

        def initmarkers():
            ps = ctypes.pointer(ascot2py.struct_c__SA_particle_state())
            self._nmrk.value = nmrk
            n_proc = ctypes.c_int32(0)
            ascot2py.prepare_markers(
                ctypes.byref(self._sim), self._nmrk, pin, ctypes.byref(ps),
                ctypes.byref(n_proc))

            ascot2py.mpi_gather_particlestate(
                ps, ctypes.byref(self._inistate), ctypes.byref(n_proc), self._nmrk,
                self._sim.mpi_rank, self._sim.mpi_size, self._sim.mpi_root)

            if self._sim.mpi_rank == self._sim.mpi_root:
                ascot2py.libascot_deallocate(ps)
            else:
                ascot2py.libascot_deallocate(self._inistate)
                self._inistate = ps

        if self._mute == "no":
            initmarkers()
        else:
            with wurlitzer.pipes() as (out, err):
                initmarkers()
            err = err.read()
            if self._mute == "err" and len(err) > 1: print(err)

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
        self._requireinit("bfield", "efield", "plasma", "wall", "boozer", "mhd",
                          "asigma", "neutral")
        if not _LIBASCOT:
            raise AscotInitException(
                "Python interface disabled as libascot.so is not found")
        if self._nmrk.value == 0:
            raise AscotInitException(
                "Initialize markers before running the simulation")
        if self._diag_occupied:
            raise AscotInitException(
                "Free previous results before running the simulation")

        # Make an inistate from markers belonging to this process
        n_proc = ctypes.c_int32(0)
        idx = ctypes.c_int32(0)
        ascot2py.mpi_my_particles(
            ctypes.byref(idx), ctypes.byref(n_proc), self._nmrk,
            self._sim.mpi_rank, self._sim.mpi_size)
        inistate = ascot2py.libascot_allocate_particle_states(n_proc.value)

        # Copy values from inistate to _inistate. Latter contains all markers
        for i in range(len(ascot2py.particle_state._fields_)):
            name = self._inistate[0]._fields_[i][0]
            for j in range(n_proc.value):
                val  = getattr(self._inistate[j], name)
                setattr(inistate[j], name, val)

        # Initialize diagnostics array and endstate
        self._endstate = ctypes.pointer(ascot2py.struct_c__SA_particle_state())
        ascot2py.diag_init(ctypes.byref(self._sim.diag_data), self._nmrk)
        self._diag_occupied = True

        # Simulate and print stdout/stderr if requested
        def runsim():
            ascot2py.offload_and_simulate(
                ctypes.byref(self._sim), self._nmrk, n_proc.value, inistate,
                ctypes.byref(self._nmrk), ctypes.byref(self._endstate))

        if self._mute == "no":
            runsim()
        else:
            with wurlitzer.pipes() as (out, err):
                runsim()
            err = err.read()
            if self._mute == "err" and len(err) > 1: print(err)

        # Print summary
        if self._sim.mpi_rank == self._sim.mpi_root and printsummary:
            ascot2py.print_marker_summary(self._endstate, self._nmrk)

        class VirtualInput():
            """Wrapper for marker and options inputs.
            """

            def __init__(self, inp):
                self.inp = inp

            def read(self):
                return self.inp

        diagorb = None
        if self._sim.diag_data.diagorb_collect:
            diagorb = self._sim.diag_data.diagorb
        dist5d = None
        if self._sim.diag_data.dist5D_collect:
            dist5d = self._sim.diag_data.dist5D
        return VirtualRun(self, self._nmrk.value,
                          self._inistate, self._endstate,
                          VirtualInput(self._virtualoptions),
                          VirtualInput(self._virtualmarkers),
                          diagorb=diagorb, dist5d=dist5d)

    def simulation_bbnbi(self, nprt, t1=0, t2=0, printsummary=True):
        """Run BBNBI simulation.

        Parameters
        ----------
        nprt : Number of markers to be launched in total.

            This number is the combined total for all injectors which are then
            distributed among the injectors in proportion to the injector power.
            For example, suppose there are 3 injectors with P_1 = 10 MW and
            P_2 = 5 MW and P_3 = 5 MW. If `nprt` = 10000, then 5000 markers are
            generated for injector 1 and 2500 markers for both injectors 2
            and 3.

            The number of ionized markers is either equal or smaller than `nprt`
            as some of the markers are lost as shinethrough.
        t1 : float, optional
            The time instant at which the injector is turned on.

            In the usual case where the background is time independent,
            the parameters `t1` and `t2` only affect the marker initial time
            which will be uniformly distributed in `[t1, t2]`.
        printsummary : bool, optional
            If True, summary of marker endstates is printed after the simulation
            completes.

        Returns
        -------
        run : :class:`VirtualBBNBI`
            An object that acts almost exactly as :class:`BBNBIGroup` except that
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
        if self._nmrk.value != 0:
            raise AscotInitException(
                "Free markers before rerunning the simulation")
        if self._diag_occupied:
            raise AscotInitException(
                "Free previous results before running the simulation")
        # Initialize diagnostics array and endstate
        self._endstate = ctypes.pointer(ascot2py.struct_c__SA_particle_state())
        ascot2py.diag_init(ctypes.byref(self._sim.diag_data), nprt)
        self._diag_occupied = True

        # Simulate and print stdout/stderr if requested
        def runsim():
            ascot2py.bbnbi_simulate(
                ctypes.byref(self._sim), nprt, t1, t2,
                ctypes.byref(self._endstate))

        if self._mute == "no":
            runsim()
        else:
            with wurlitzer.pipes() as (out, err):
                runsim()
            err = err.read()
            if self._mute == "err" and len(err) > 1: print(err)

        # Print summary
        self._nmrk.value = nprt
        if self._sim.mpi_rank == self._sim.mpi_root and printsummary:
            ascot2py.print_marker_summary(self._endstate, self._nmrk)

        class VirtualInput():
            """Wrapper for marker and options inputs.
            """

            def __init__(self, inp):
                self.inp = inp

            def read(self):
                return self.inp

        diagorb = None
        if self._sim.diag_data.diagorb_collect:
            diagorb = self._sim.diag_data.diagorb
        dist5d = None
        if self._sim.diag_data.dist5D_collect:
            dist5d = self._sim.diag_data.dist5D
        return VirtualBBNBIRun(
            self, nprt, self._endstate, VirtualInput(self._virtualoptions),
            dist5d=dist5d)

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

        if inputs:
            self.input_free()
        if markers and self._nmrk.value > 0:
            self._nmrk.value = 0
            ascot2py.libascot_deallocate(self._inistate)
            self._virtualmarkers = None
        if diagnostics and self._diag_occupied:
            self._diag_occupied = False
            ascot2py.diag_free(ctypes.byref(self._sim.diag_data))
            ascot2py.libascot_deallocate(self._endstate)
