/**
 * @file hdf5_options.c
 * @brief Read options from HDF5 file
 */
#include <stdlib.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../diag.h"
#include "../diag/diag_orb.h"
#include "../diag/dist_5D.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_rho5D.h"
#include "../diag/dist_rho6D.h"
#include "../endcond.h"
#include "../math.h"
#include "../simulate.h"
#include "hdf5_helpers.h"
#include "hdf5_options.h"

#define OPTPATH /**< Macro that is used to store paths to data groups */

int hdf5_options_read_dist5D(hid_t file, dist_5D_offload_data* dist,
                             char* qid);
int hdf5_options_read_dist6D(hid_t file, dist_6D_offload_data* dist,
                             char* qid);
int hdf5_options_read_distrho5D(hid_t file, dist_rho5D_offload_data* dist,
                                char* qid);
int hdf5_options_read_distrho6D(hid_t file, dist_rho6D_offload_data* dist,
                                char* qid);
int hdf5_options_read_diagorb(hid_t file, diag_orb_offload_data* diagorb,
                              char* qid);

/**
 * @brief Read options and diagnostics settings from HDF5 file
 *
 * This function reads options with given qid from HDF5 file. The file is opened
 * and closed outside of this function.
 *
 * The options are read directly to simulation offload data.
 *
 * @param file the file where data is read
 * @param sim pointer to simulation offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read(hid_t file, sim_offload_data* sim, char* qid){

    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    real tempfloat; // For reading float data that is converted to int.

    if( hdf5_read_double(OPTPATH "SIM_MODE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->sim_mode = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_ADAPTIVE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->enable_ada = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "RECORD_GO_AS_GC", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->record_GOasGC = (int)tempfloat;


    if( hdf5_read_double(OPTPATH "FIXEDSTEP_USE_USERDEFINED", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->fix_usrdef_use = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "FIXEDSTEP_USERDEFINED", &sim->fix_usrdef_val,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "FIXEDSTEP_NSTEPS_PER_GYROTIME", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->fix_stepsPerGO = (int)tempfloat;


    if( hdf5_read_double(OPTPATH "ADAPTIVE_TOL_ORBIT", &sim->ada_tol_orbfol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ADAPTIVE_TOL_CCOL", &sim->ada_tol_clmbcol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ADAPTIVE_MAX_DRHO", &sim->ada_max_drho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ADAPTIVE_MAX_DPHI", &sim->ada_max_dphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}


    if( hdf5_read_double(OPTPATH "ENABLE_ORBIT_FOLLOWING", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->enable_orbfol = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_COULOMB_COLLISIONS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->enable_clmbcol = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "DISABLE_FIRSTORDER_GCTRANS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->disable_gctransform = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "DISABLE_ENERGY_CCOLL", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->disable_energyccoll = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "DISABLE_PITCH_CCOLL", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->disable_pitchccoll = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "DISABLE_GCDIFF_CCOLL", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->disable_gcdiffccoll = (int)tempfloat;

    int ec;
    sim->endcond_active = 0;

    if( hdf5_read_double(OPTPATH "ENDCOND_SIMTIMELIM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_tmax * (ec > 0);
    if( hdf5_read_double(OPTPATH "ENDCOND_CPUTIMELIM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_cpumax * (ec > 0);
    if( hdf5_read_double(OPTPATH "ENDCOND_RHOLIM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_rhomin * (ec > 0);
    sim->endcond_active = sim->endcond_active | endcond_rhomax * (ec > 0);
    if( hdf5_read_double(OPTPATH "ENDCOND_ENERGYLIM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_emin * (ec > 0);
    sim->endcond_active = sim->endcond_active | endcond_therm * (ec > 0);
    if( hdf5_read_double(OPTPATH "ENDCOND_WALLHIT", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_wall * (ec > 0);
    if( hdf5_read_double(OPTPATH "ENDCOND_MAXORBS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_polmax * (ec > 0);
    sim->endcond_active = sim->endcond_active | endcond_tormax * (ec > 0);


    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_SIM_TIME",
                         &sim->endcond_maxSimTime,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_CPU_TIME",
                         &sim->endcond_maxCpuTime,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_RHO",
                         &sim->endcond_maxRho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MIN_RHO",
                         &sim->endcond_minRho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MIN_ENERGY",
                         &sim->endcond_minEkin,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MIN_ENERGY_TIMES_THERMAL",
                         &sim->endcond_minEkinPerTi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_minEkin *= CONST_E; // eV -> J

    int temp;
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_POLOIDALORBS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    temp = (int)tempfloat;
    sim->endcond_maxPolOrb = temp * 2 *CONST_PI;

    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_TOROIDALORBS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    temp = (int)tempfloat;
    sim->endcond_maxTorOrb = temp * 2 *CONST_PI;


    /* See which diagnostics are active */
    diag_offload_data* diag = &sim->diag_offload_data;

    if( hdf5_read_double(OPTPATH "ENABLE_DIST_5D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->dist5D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_6D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->dist6D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_rho5D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->distrho5D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_rho6D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->distrho6D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_ORBITWRITE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->diagorb_collect = (int)tempfloat;

    /* Read individual diagnostics data */
    if(diag->dist5D_collect) {
        if( hdf5_options_read_dist5D(file, &diag->dist5D, qid) ) {
            return 1;
        }
    }
    if(diag->dist6D_collect) {
        if( hdf5_options_read_dist6D(file, &diag->dist6D, qid) ) {
            return 1;
        }
    }
    if(diag->distrho5D_collect) {
        if( hdf5_options_read_distrho5D(file, &diag->distrho5D, qid) ) {
            return 1;
        }
    }
    if(diag->distrho6D_collect) {
        if( hdf5_options_read_distrho6D(file, &diag->distrho6D, qid) ) {
            return 1;
        }
    }
    if(diag->diagorb_collect) {
        diag->diagorb.record_mode = sim->sim_mode;
        if(sim->record_GOasGC && (sim->sim_mode == simulate_mode_fo ||
                                  sim->sim_mode == simulate_mode_hybrid) ) {
            diag->diagorb.record_mode = simulate_mode_gc;
        }

        if( hdf5_options_read_diagorb(file, &diag->diagorb, qid) ) {
            return 1;
        }
    }

    return 0;
}

/**
 * @brief Helper function to read dist5D settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param dist pointer to dist5D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_dist5D(hid_t file, dist_5D_offload_data* dist,
                             char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_R", &dist->min_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_R", &dist->max_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_R", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_r = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_phi", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_z", &dist->min_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_z", &dist->max_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_z", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_z = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpa", &dist->min_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpa", &dist->max_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vpa", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vpara = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpe", &dist->min_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpe", &dist->max_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vpe", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vperp = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_t", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_q", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_q", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_q", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_q = (int)tempfloat;

    return 0;
}

/**
 * @brief Helper function to read dist6D settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param dist pointer to dist6D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_dist6D(hid_t file, dist_6D_offload_data* dist,
                             char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_R", &dist->min_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_R", &dist->max_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_R", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_r = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_phi", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_z", &dist->min_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_z", &dist->max_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_z", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_z = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vR", &dist->min_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vR", &dist->max_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vR", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vr = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vphi", &dist->min_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vphi", &dist->max_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vphi", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vphi = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vz", &dist->min_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vz", &dist->max_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vz", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vz = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_t", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_q", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_q", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_q", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_q = (int)tempfloat;

    return 0;
}

/**
 * @brief Helper function to read distrho5D settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param dist pointer to distrho5D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_distrho5D(hid_t file, dist_rho5D_offload_data* dist,
                                char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_rho", &dist->min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_rho", &dist->max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_rho", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_rho = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_phi", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_theta", &dist->min_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_theta", &dist->max_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_theta", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_theta = (int)tempfloat;
    dist->min_theta = math_deg2rad(dist->min_theta);
    dist->max_theta = math_deg2rad(dist->max_theta);

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpa", &dist->min_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpa", &dist->max_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vpa", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vpara = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpe", &dist->min_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpe", &dist->max_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vpe", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vperp = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_t", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_q", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_q", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_q", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_q = (int)tempfloat;

    return 0;
}

/**
 * @brief Helper function to read dist5D settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param dist pointer to distrho6D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_distrho6D(hid_t file, dist_rho6D_offload_data* dist,
                                char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_rho", &dist->min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_rho", &dist->max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_rho", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_rho = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_phi", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_theta", &dist->min_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_theta", &dist->max_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_theta", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_theta = (int)tempfloat;
    dist->min_theta = math_deg2rad(dist->min_theta);
    dist->max_theta = math_deg2rad(dist->max_theta);

    if( hdf5_read_double(OPTPATH "DIST_MIN_vR", &dist->min_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vR", &dist->max_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vR", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vr = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vphi", &dist->min_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vphi", &dist->max_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vphi", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vphi = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_vz", &dist->min_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vz", &dist->max_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_vz", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_vz = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_t", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_q", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_q", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_q", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_q = (int)tempfloat;

    return 0;
}

/**
 * @brief Helper function to read orbit diagnostics settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param diagorb pointer to orbit diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_diagorb(hid_t file, diag_orb_offload_data* diagorb,
                              char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "ORBITWRITE_MODE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diagorb->mode = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "ORBITWRITE_MAXPOINTS", &tempfloat,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    diagorb->Npnt = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ORBITWRITE_INTERVAL",
                         &(diagorb->writeInterval),
                         file, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        diagorb->toroidalangles[i] = 361;
        diagorb->poloidalangles[i] = 361;
    }

    if( hdf5_read_double(OPTPATH "ORBITWRITE_TOROIDALANGLES",
                         (diagorb->toroidalangles),
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ORBITWRITE_POLOIDALANGLES",
                         (diagorb->poloidalangles),
                         file, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        if(diagorb->toroidalangles[i] == 361){
            diagorb->ntoroidalplots = i;
            break;
        }
    }
    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        if(diagorb->poloidalangles[i] == 361){
            diagorb->npoloidalplots = i;
            break;
        }
    }

    for(int i=0; i < diagorb->ntoroidalplots; i++) {
        diagorb->toroidalangles[i] = diagorb->toroidalangles[i]*CONST_PI/180;
    }
    for(int i=0; i < diagorb->npoloidalplots; i++) {
        diagorb->poloidalangles[i] = diagorb->poloidalangles[i]*CONST_PI/180;
    }

    return 0;
}
