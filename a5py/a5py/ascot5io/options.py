"""
Options IO.

File: options.py
"""
import h5py
import numpy as np
import textwrap

from . ascot5file import add_group, get_qid
from a5py.ascot5io.ascot5data import AscotData

def write_hdf5(fn, options, desc=None):
    """
    Write options.

    Unlike most other "write" functions, this one takes dictionary as an
    argument. The dictionary should have exactly the same format
    as given by the read_hdf5() or get_default() functions in this module.

    Args:
        fn : str <br>
            Full path to HDF5 file.
        options : dictionary <br>
            Options to be written in dictionary format.
        desc : str, optional <br>
            Description for this input.
    """

    parent = "options"
    group  = "opt"
    qid    = None

    info = []

    

    with h5py.File(fn, "a") as f:
        g    = add_group(f, parent, group, desc=desc)
        name = g.name

        for opt in options:
            if opt != "qid" and opt != "date" and opt != "description":
                data = options[opt]
                if type(data) is not np.array:
                    data = np.array(data)

                data =data.astype("f8")
                d = g.create_dataset(opt, (data.size,), data=data)
                d.attrs["info"] = np.string_(info[opt])

    return name


def read_hdf5(fn, qid):
    """
    Read options from HDF5 file.

    Args:
        fn : str <br>
            Full path to the HDF5 file.
        qid : str <br>
            qid of the options to be read.

    Returns:
        Dictionary containing options.
    """

    path = "options" + "/opt-" + qid

    with h5py.File(fn,"r") as f:
        out = {}

        # Metadata.
        out["qid"]  = qid
        out["date"] = f[path].attrs["date"]
        out["description"] = f[path].attrs["description"]

        # Actual data.
        for opt in f[path]:
            out[opt] = f[path][opt][:]

    return out


def get_default():
    """
    Get default option parameters.

    Returns:
        A dictionary containing default parameters.
    """

    #**************************************************************************#
    #*                      SIMULATION MODE AND TIME-STEP                     *#
    #*                                                                        *#
    #**************************************************************************#
    info.append(("SIM_MODE", textwrap.dedent("""\
        Simulation mode (1, 2, 3, 4)
        - 1 Gyro-orbit
        - 2 Guiding center
        - 3 Hybrid
        - 4 Magnetic field lines
    """)))
    info.append(("ENABLE_ADAPTIVE", textwrap.dedent("""\
        Use adaptive time-step (0, 1)
        This option is used only if SIM_MODE = 2 or 3. Gyro-orbit
        simulations are always done with fixed time-step and magnetic field line
        simulations with adaptive time-step.
        - 0 Use fixed time-step
        - 1 Use adaptive time-step
    """)))
    info.append(("RECORD_GO_AS_GC", textwrap.dedent("""\
        Record GOs as GCs in diagnostics (0, 1)
        - 0 Record GOs as GOs
        - 1 Record GOs as GCs
    """)))
    info.append(("FIXEDSTEP_USE_USERDEFINED", textwrap.dedent("""\
        Define fixed time-step value explicitly (0,1)
        Note: The adaptive scheme uses fixed time-step value as an initial step.
        - 0 Calculate time-step from FIXEDSTEP_NSTEPS_PER_GYROTIME
        - 1 Use opt.opt.FIXEDSTEP_USERDEFINED as a time-step
    """)))
    info.append(("FIXEDSTEP_USERDEFINED", textwrap.dedent("""\
        User-defined time-step [s]
    """)))
    info.append(("FIXEDSTEP_NSTEPS_PER_GYROTIME", textwrap.dedent("""\
        Time-step is 2pi / ( gyrofrequency * N ) where N is this parameter
    """)))
    info.append(("ADAPTIVE_TOL_ORBIT", textwrap.dedent("""\
        Relative error tolerance for orbit following in adaptive scheme
    """)))
    info.append(("ADAPTIVE_TOL_CCOL", textwrap.dedent("""\
        Relative error tolerance for Coulomb collisions in adaptive scheme
    """)))
    info.append(("ADAPTIVE_MAX_DRHO", textwrap.dedent("""\
        Maximum allowed change in rho during one time-step in adaptive scheme
    """)))
    info.append(("ADAPTIVE_MAX_DPHI", textwrap.dedent("""\
        Maximum allowed change in phi during one time-step in adaptive scheme
    """)))

    #**************************************************************************#
    #*                             END CONDITIONS                             *#
    #*                                                                        *#
    #**************************************************************************#
    info.append(("ENDCOND_SIMTIMELIM", textwrap.dedent("""\
        Terminate when marker's clock ("laboratory") time reaches a limit
        The limit is set by ENDCOND_MAX_SIM_TIME
    """)))
    info.append(("ENDCOND_CPUTIMELIM", textwrap.dedent("""\
        Terminate when marker's cpu time reaches a limit
        The limit is set by ENDCOND_MAX_CPU_TIME
    """)))
    info.append(("ENDCOND_RHOLIM", textwrap.dedent("""\
        Terminate if marker goes outside given rho boundaries
        rho boundaries are defined by ENDCOND_MAX_RHO and ENDCOND_MAX_RHO.
    """)))
    info.append(("ENDCOND_ENERGYLIM", textwrap.dedent("""\
        Terminate when marker energy is below a user-specified value
        The user specified values are ENDCOND_MIN_ENERGY and
        ENDCOND_MIN_ENERGY_TIMES_THERMAL. Marker is terminated when either
        of these limits is reached.
    """)))
    info.append(("ENDCOND_WALLHIT", textwrap.dedent("""\
        Terminate when marker impacts wall
    """)))
    info.append(("ENDCOND_MAXORBS", textwrap.dedent("""\
        Terminate when marker has completed user-specified number of orbits
        Limit ENDCOND_MAX_TOROIDALORBS is used for a number of toroidal
        and ENDCOND_MAX_POLOIDALORBS for poloidal orbits. Marker is
        terminated when either of these limits is reached.
    """)))
    info.append(("ENDCOND_MAX_SIM_TIME", textwrap.dedent("""\
        Maximum simulation time [s]
    """)))
    info.append(("ENDCOND_MAX_CPU_TIME", textwrap.dedent("""\
        Maximum cpu time [s]
    """)))
    info.append(("ENDCOND_MAX_RHO", textwrap.dedent("""\
        Maximum rho value
    """)))
    info.append(("ENDCOND_MIN_RHO", textwrap.dedent("""\
        Minimum rho value
    """)))
    info.append(("ENDCOND_MIN_ENERGY", textwrap.dedent("""\
        Minimum energy [eV]
    """)))
    info.append(("ENDCOND_MIN_ENERGY_TIMES_THERMAL", textwrap.dedent("""\
        Minimum energy limit local electron thermal energy times this value
    """)))
    info.append(("ENDCOND_MAX_TOROIDALORBS", textwrap.dedent("""\
        Maximum number of toroidal orbits
    """)))
    info.append(("ENDCOND_MAX_POLOIDALORBS", textwrap.dedent("""\
        Maximum number of poloidal orbits
    """)))

    #**************************************************************************#
    #*                               PHYSICS                                  *#
    #*                                                                        *#
    #**************************************************************************#

    info.append(("ENABLE_ORBIT_FOLLOWING", textwrap.dedent("""\
        Trace markers in an electromagnetic field
    """)))
    info.append(("ENABLE_COULOMB_COLLISIONS", textwrap.dedent("""\
        Markers experience Coulomb collisions with background plasma
    """)))
    info.append(("DISABLE_FIRSTORDER_GCTRANS", textwrap.dedent("""\
        Disable first order guiding center transformation in velocity space
    """)))
    info.append(("DISABLE_ENERGY_CCOLL", textwrap.dedent("""\
        Disable guiding center energy collisions
    """)))
    info.append(("DISABLE_PITCH_CCOLL", textwrap.dedent("""\
        Disable guiding center pitch collisions
    """)))
    info.append(("DISABLE_GCDIFF_CCOLL", textwrap.dedent("""\
        Disable guiding center spatial diffusion
    """)))

    #**************************************************************************#
    #*                            DISTRIBUTIONS                               *#
    #*                                                                        *#
    #**************************************************************************#

    info.append(
        ("ENABLE_DIST_5D",
         """\
         Collect distribution histogram in [R, phi, z, vpa, vpe, t, q]
         The coordinates are
           - R   major radius
           - phi toroidal angle
           - z   z-coordinate
           - vpa velocity component parallel to magnetic field
           - vpe velocity component perpendicular to magnetic field
           - t   time
           - q   charge
         """,
        )
    )
    info.append(
        ("ENABLE_DIST_6D",
         """\
         Collect distribution histogram in [R, phi, z, vR, vphi, vz, t, q]
         The coordinates are
             - R    major radius
             - phi  toroidal angle
             - z    z-coordinate
             - vR   velocity R-component
             - vphi velocity phi-component
             - vz   velocity z-component
             - t    time
             - q    charge
         """,
        )
    )
    info.append(
        ("ENABLE_DIST_rho5D",
         """\
         Collect distribution histogram in [rho, pol, phi, vpa, vpe, t, q]
         The coordinates are
             - rho  flux surface
             - pol  poloidal angle
             - phi  toroidal angle
             - z    z-coordinate
             - vpa  velocity component parallel to magnetic field
             - vpe  velocity component perpendicular to magnetic field
             - t    time
             - q    charge
         """,
        )
    )
    info.append(
        ("ENABLE_DIST_rho6D",
         """\
         Collect distribution histogram in [rho, pol, phi, vR, vphi, vz, t, q]
         The coordinates are
             - rho  flux surface
             - pol  poloidal angle
             - phi  toroidal angle
             - z    z-coordinate
             - vR   velocity R-component
             - vphi velocity phi-component
             - vz   velocity z-component
             - t    time
             - q    charge
         """,
        )
    )
    info.append(("DIST_MIN_R", textwrap.dedent("""\
        Minimum bin edge for major R coordinate [m]
    """)))
    info.append(("DIST_MAX_R", textwrap.dedent("""\
        Maximum bin edge for R coordinate [m]
    """)))
    info.append(("DIST_NBIN_R", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_R, DIST_MAX_R] is divided to
    """)))
    info.append(("DIST_MIN_phi", textwrap.dedent("""\
        Minimum bin edge for phi coordinate [deg]
    """)))
    info.append(("DIST_MAX_phi", textwrap.dedent("""\
        Maximum bin edge for phi coordinate [deg]
    """)))
    info.append(("DIST_NBIN_phi", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_phi, DIST_MAX_phi] is divided to
    """)))
    info.append(("DIST_MIN_z", textwrap.dedent("""\
        Minimum bin edge for z coordinate [m]
    """)))
    info.append(("DIST_MAX_z", textwrap.dedent("""\
        Maximum bin edge for z coordinate [m]
    """)))
    info.append(("DIST_NBIN_z", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_z, DIST_MAX_z] is divided to
    """)))
    info.append(("DIST_MIN_rho", textwrap.dedent("""\
        Minimum bin edge for rho coordinate
    """)))
    info.append(("DIST_MAX_rho", textwrap.dedent("""\
        Maximum bin edge for rho coordinate
    """)))
    info.append(("DIST_NBIN_rho", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_rho, DIST_MAX_rho] is divided to
    """)))
    info.append(("DIST_MIN_pol", textwrap.dedent("""\
        Minimum bin edge for pol coordinate [deg]
    """)))
    info.append(("DIST_MAX_pol", textwrap.dedent("""\
        Maximum bin edge for pol coordinate [deg]
    """)))
    info.append(("DIST_NBIN_pol", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_pol, DIST_MAX_pol] is divided to
    """)))
    info.append(("DIST_MIN_vpa", textwrap.dedent("""\
        Minimum bin edge for vpa coordinate [m/s]
    """)))
    info.append(("DIST_MAX_vpa", textwrap.dedent("""\
        Maximum bin edge for vpa coordinate [m/s]
    """)))
    info.append(("DIST_NBIN_vpa", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_vpa, DIST_MAX_vpa] is divided to
    """)))
    info.append(("DIST_MIN_vpe", textwrap.dedent("""\
        Minimum bin edge for vpe coordinate [m/s]
    """)))
    info.append(("DIST_MAX_vpe", textwrap.dedent("""\
        Maximum bin edge for vpe coordinate [m/s]
    """)))
    info.append(("DIST_NBIN_vpe", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_vpe, DIST_MAX_vpe] is divided to
    """)))
    info.append(("DIST_MIN_vR", textwrap.dedent("""\
        Minimum bin edge for vR coordinate [m/s]
    """)))
    info.append(("DIST_MAX_vR", textwrap.dedent("""\
        Maximum bin edge for vR coordinate [m/s]
    """)))
    info.append(("DIST_NBIN_vR", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_vR, DIST_MAX_vR] is divided to
    """)))
    info.append(("DIST_MIN_vphi", textwrap.dedent("""\
        Minimum bin edge for vphi coordinate [m/s]
    """)))
    info.append(("DIST_MAX_vphi", textwrap.dedent("""\
        Maximum bin edge for vphi coordinate [m/s]
    """)))
    info.append(("DIST_NBIN_vphi", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_vphi, DIST_MAX_vphi] is divided to
    """)))
    info.append(("DIST_MIN_vz", textwrap.dedent("""\
        Minimum bin edge for vz coordinate [m/s]
    """)))
    info.append(("DIST_MAX_vz", textwrap.dedent("""\
        Maximum bin edge for vz coordinate [m/s]
    """)))
    info.append(("DIST_NBIN_vz", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_vz, DIST_MAX_vz] is divided to
    """)))
    info.append(("DIST_MIN_t", textwrap.dedent("""\
        Minimum bin edge for t coordinate [s]
    """)))
    info.append(("DIST_MAX_t", textwrap.dedent("""\
        Maximum bin edge for t coordinate [s]
    """)))
    info.append(("DIST_NBIN_t", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_t, DIST_MAX_t] is divided to
    """)))
    info.append(("DIST_MIN_q", textwrap.dedent("""\
        Minimum bin edge for q coordinate [e]
    """)))
    info.append(("DIST_MAX_q", textwrap.dedent("""\
        Maximum bin edge for q coordinate [e]
    """)))
    info.append(("DIST_NBIN_q", textwrap.dedent("""\
        Number of bins the interval [DIST_MIN_q, DIST_MAX_q] is divided to
    """)))

    #**************************************************************************#
    #*                             ORBIT WRITE                                *#
    #*                                                                        *#
    #**************************************************************************#

    info.append(
        ("ENABLE_ORBITWRITE",
         """\
         Enable diagnostics that store marker orbit
             - 0 Marker orbit diagnostics are not collected
             - 1 Marker orbit diagnostics are collected
         """,
         1)
    )
    info.append(
        ("ORBITWRITE_MODE",
         """\
         What kind of marker orbit diagnostics are collected
         These are only used if ENABLE_ORBITWRITE is active.
             - 0 When marker crosses a plane (Poincare-plot)
             - 1 Between given time intervals
         """,
         1)
    )
    info.append(
        ("ORBITWRITE_MAXPOINTS",
         """\
         Maximum number of points (per marker) to be written
         If this number is exceeded when marker is being simulated, the oldest
         points will be replaced as long as the simulation continues. Thus, this
         parameter is effectively the number of marker's last positions that are
         stored.
         """,
         10)
    )
    info.append(
        ("ORBITWRITE_TOROIDALANGLES",
         """\
         Poloidal angles of toroidal planes where toroidal plots are collected
         Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 0.
         """,
         [0, 180])
    )
    info.append(
        ("ORBITWRITE_POLOIDALANGLES",
         """\
         Toroidal angles of poloidal planes where poloidal plots are collected
         Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 0.
         """,
         [0, 180])
    )
    info.append(
        ("ORBITWRITE_INTERVAL",
         """\
         Time interval for writing marker state [s]
         Used when ENABLE_ORBITWRITE = 1 and ORBITWRITE_MODE = 1.
         """,
         1e-8)
    )

    opt = {}
    #**************************************************************************#
    #*                      SIMULATION MODE AND TIME-STEP                     *#
    #*                                                                        *#
    #**************************************************************************#

    opt["SIM_MODE"]        = 3
    opt["ENABLE_ADAPTIVE"] = 1
    opt["RECORD_GO_AS_GC"] = 0
    opt["FIXEDSTEP_USE_USERDEFINED"]     = 0
    opt["FIXEDSTEP_USERDEFINED"]         = 1e-8
    opt["FIXEDSTEP_NSTEPS_PER_GYROTIME"] = 20
    opt["ADAPTIVE_TOL_ORBIT"] = 1e-8
    opt["ADAPTIVE_TOL_CCOL"]  = 1e-1
    opt["ADAPTIVE_MAX_DRHO"]  = 0.1
    opt["ADAPTIVE_MAX_DPHI"]  = 5

    #**************************************************************************#
    #*                             END CONDITIONS                             *#
    #*                                                                        *#
    #**************************************************************************#

    opt["ENDCOND_SIMTIMELIM"]       = 1
    opt["ENDCOND_CPUTIMELIM"]       = 1
    opt["ENDCOND_RHOLIM"]           = 0
    opt["ENDCOND_ENERGYLIM"]        = 1
    opt["ENDCOND_WALLHIT"]          = 1
    opt["ENDCOND_MAXORBS"]          = 0
    opt["ENDCOND_MAX_SIM_TIME"]     = 1.0
    opt["ENDCOND_MAX_CPU_TIME"]     = 3600
    opt["ENDCOND_MAX_RHO"]          = 0
    opt["ENDCOND_MIN_RHO"]          = 1
    opt["ENDCOND_MIN_ENERGY"]       = 1e3
    opt["ENDCOND_MIN_ENERGY_TIMES_THERMAL"] = 1
    opt["ENDCOND_MAX_TOROIDALORBS"] = 100
    opt["ENDCOND_MAX_POLOIDALORBS"] = 100

    #**************************************************************************#
    #*                               PHYSICS                                  *#
    #*                                                                        *#
    #**************************************************************************#

    opt["ENABLE_ORBIT_FOLLOWING"]     = 1
    opt["ENABLE_COULOMB_COLLISIONS"]  = 1
    opt["DISABLE_FIRSTORDER_GCTRANS"] = 0
    opt["DISABLE_ENERGY_CCOLL"]       = 0
    opt["DISABLE_PITCH_CCOLL"]        = 0
    opt["DISABLE_GCDIFF_CCOLL"]       = 0

    #**************************************************************************#
    #*                            DISTRIBUTIONS                               *#
    #*                                                                        *#
    #**************************************************************************#

    opt["ENABLE_DIST_5D"]    = 1
    opt["ENABLE_DIST_6D"]    = 1
    opt["ENABLE_DIST_rho5D"] = 1
    opt["ENABLE_DIST_rho6D"] = 1
    opt["DIST_MIN_R"]     = 0.1
    opt["DIST_MAX_R"]     = 10
    opt["DIST_NBIN_R"]    = 10
    opt["DIST_MIN_phi"]   = 0
    opt["DIST_MAX_phi"]   = 360
    opt["DIST_NBIN_phi"]  = 10
    opt["DIST_MIN_z"]     = -5
    opt["DIST_MAX_z"]     = 5
    opt["DIST_NBIN_z"]    = 10
    opt["DIST_MIN_rho"]   = 0
    opt["DIST_MAX_rho"]   = 2
    opt["DIST_NBIN_rho"]  = 10
    opt["DIST_MIN_pol"]   = 0
    opt["DIST_MAX_pol"]   = 360
    opt["DIST_NBIN_pol"]  = 10
    opt["DIST_MIN_vpa"]   = -3e8
    opt["DIST_MAX_vpa"]   = 3e8
    opt["DIST_NBIN_vpa"]  = 10
    opt["DIST_MIN_vpe"]   = 0
    opt["DIST_MAX_vpe"]   = 3e8
    opt["DIST_NBIN_vpe"]  = 10
    opt["DIST_MIN_vR"]    = 0
    opt["DIST_MAX_vR"]    = 3e8
    opt["DIST_NBIN_vR"]   = 10
    opt["DIST_MIN_vphi"]  = 0
    opt["DIST_MAX_vphi"]  = 3e8
    opt["DIST_NBIN_vphi"] = 10
    opt["DIST_MIN_vz"]    = 0
    opt["DIST_MAX_vz"]    = 3e8
    opt["DIST_NBIN_vz"]   = 10
    opt["DIST_MIN_t"]     = 0
    opt["DIST_MAX_t"]     = 1
    opt["DIST_NBIN_t"]    = 10
    opt["DIST_MIN_q"]     = -100
    opt["DIST_MAX_q"]     = 100
    opt["DIST_NBIN_q"]    = 10


class Opt(AscotData):

    def read(self):
        return read_hdf5(self._file, self.get_qid())
