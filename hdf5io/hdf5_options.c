/**
 * @file hdf5_options.c
 * @brief Read options from HDF5 file
 */
#include <stdlib.h>
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
#include <hdf5.h>
#include "hdf5_helpers.h"
#include <hdf5_hl.h>
#include "hdf5_options.h"

int hdf5_options_read_dist5D(hid_t file, dist_5D_offload_data* dist,
                             char* qid);
int hdf5_options_read_dist6D(hid_t file, dist_6D_offload_data* dist,
                             char* qid);
int hdf5_options_read_distrho5D(hid_t file, dist_rho5D_offload_data* dist,
                                char* qid);
int hdf5_options_read_distrho6D(hid_t file, dist_rho6D_offload_data* dist,
                                char* qid);
int hdf5_options_read_orbits(hid_t file, diag_orb_offload_data* orbits,
                             char* qid);

/**
 * @brief Read options and diagnostics settings from HDF5 file
 *
 * This function reads options with given qid from HDF5 file. The file is opened
 * and closed outside of this function.
 *
 * The options are read directly to simulation offload data.
 *
 * @param f the file where data is read
 * @param sim pointer to simulation offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read(hid_t file, sim_offload_data* sim, char* qid){

    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    if( hdf5_read_int(OPTPATH "SIM_MODE", &sim->sim_mode,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ENABLE_ADAPTIVE", &sim->enable_ada,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "RECORD_GO_AS_GC", &sim->record_GOasGC,
                      file, qid, __FILE__, __LINE__) ) {return 1;}


    if( hdf5_read_int(OPTPATH "FIXEDSTEP_USE_USERDEFINED", &sim->fix_usrdef_use,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "FIXEDSTEP_USERDEFINED", &sim->fix_usrdef_val,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "FIXEDSTEP_NSTEPS_PER_GYROTIME",
                      &sim->fix_stepsPerGO,
                      file, qid, __FILE__, __LINE__) ) {return 1;}


    if( hdf5_read_double(OPTPATH "ADAPTIVE_TOL_ORBIT", &sim->ada_tol_orbfol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ADAPTIVE_TOL_CCOL", &sim->ada_tol_clmbcol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ADAPTIVE_MAX_DRHO", &sim->ada_max_drho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ADAPTIVE_MAX_DPHI", &sim->ada_max_dphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}


    if( hdf5_read_int(OPTPATH "ENABLE_ORBIT_FOLLOWING", &sim->enable_orbfol,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ENABLE_COULOMB_COLLISIONS",
                      &sim->enable_clmbcol,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DISABLE_FIRSTORDER_GCTRANS",
                      &sim->disable_gctransform,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    int ec;
    sim->endcond_active = 0;

    if( hdf5_read_int(OPTPATH "ENDCOND_SIMTIMELIM", &ec,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_active = sim->endcond_active | endcond_tmax * (ec > 0);
    if( hdf5_read_int(OPTPATH "ENDCOND_CPUTIMELIM", &ec,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_active = sim->endcond_active | endcond_cpumax * (ec > 0);
    if( hdf5_read_int(OPTPATH "ENDCOND_RHOLIM", &ec,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_active = sim->endcond_active | endcond_rhomin * (ec > 0);
    sim->endcond_active = sim->endcond_active | endcond_rhomax * (ec > 0);
    if( hdf5_read_int(OPTPATH "ENDCOND_ENERGYLIM", &ec,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_active = sim->endcond_active | endcond_emin * (ec > 0);
    sim->endcond_active = sim->endcond_active | endcond_therm * (ec > 0);
    if( hdf5_read_int(OPTPATH "ENDCOND_WALLHIT", &ec,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_active = sim->endcond_active | endcond_wall * (ec > 0);
    if( hdf5_read_int(OPTPATH "ENDCOND_MAXORBS", &ec,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
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
                         &sim->endcond_minEkinPerTe,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_minEkin *= CONST_E; // eV -> J

    int temp;
    if( hdf5_read_int(OPTPATH "ENDCOND_MAX_POLOIDALORBS", &temp,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_maxPolOrb = temp * 2 *CONST_PI;

    if( hdf5_read_int(OPTPATH "ENDCOND_MAX_TOROIDALORBS", &temp,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_maxTorOrb = temp * 2 *CONST_PI;


    /* See which diagnostics are active */
    diag_offload_data* diag = &sim->diag_offload_data;

    if( hdf5_read_int(OPTPATH "ENABLE_R_phi_z_vpa_vpe_t_q_DIST",
                      &diag->dist5D_collect,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST",
                      &diag->dist6D_collect,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ENABLE_rho_pol_phi_vpa_vpe_t_q_DIST",
                      &diag->distrho5D_collect,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ENABLE_rho_pol_phi_vR_vphi_vz_t_q_DIST",
                      &diag->distrho6D_collect,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ENABLE_ORBITWRITE", &diag->orb_collect,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

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
    if(diag->orb_collect) {
        if( hdf5_options_read_orbits(file, &diag->orbits, qid) ) {
            return 1;
        }
    }

    return 0;
}

/**
 * @brief Helper function to read dist5D settings from HDF5 file
 *
 * @param f the file where settings are read
 * @param dist pointer to dist5D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_dist5D(hid_t file, dist_5D_offload_data* dist,
                             char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    if( hdf5_read_double(OPTPATH "DIST_MIN_R", &dist->min_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_R", &dist->max_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_R", &dist->n_r,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_phi", &dist->n_phi,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_z", &dist->min_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_z", &dist->max_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_z", &dist->n_z,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpa", &dist->min_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpa", &dist->max_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vpa", &dist->n_vpara,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpe", &dist->min_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpe", &dist->max_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vpe", &dist->n_vperp,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_t", &dist->n_time,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    int charge;
    if( hdf5_read_int(OPTPATH "DIST_MIN_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_MAX_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->max_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_NBIN_q", &dist->n_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Helper function to read dist6D settings from HDF5 file
 *
 * @param f the file where settings are read
 * @param dist pointer to dist6D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_dist6D(hid_t file, dist_6D_offload_data* dist,
                             char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    if( hdf5_read_double(OPTPATH "DIST_MIN_R", &dist->min_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_R", &dist->max_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_R", &dist->n_r,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_phi", &dist->n_phi,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_z", &dist->min_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_z", &dist->max_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_z", &dist->n_z,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vR", &dist->min_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vR", &dist->max_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vR", &dist->n_vr,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vphi", &dist->min_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vphi", &dist->max_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vphi", &dist->n_vphi,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vz", &dist->min_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vz", &dist->max_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vz", &dist->n_vz,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_t", &dist->n_time,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    int charge;
    if( hdf5_read_int(OPTPATH "DIST_MIN_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_MAX_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->max_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_NBIN_q", &dist->n_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Helper function to read distrho5D settings from HDF5 file
 *
 * @param f the file where settings are read
 * @param dist pointer to distrho5D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_distrho5D(hid_t file, dist_rho5D_offload_data* dist,
                                char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    if( hdf5_read_double(OPTPATH "DIST_MIN_rho", &dist->min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_rho", &dist->max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_rho", &dist->n_rho,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_phi", &dist->n_phi,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_pol", &dist->min_pol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_pol", &dist->max_pol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_pol", &dist->n_pol,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_pol = math_deg2rad(dist->min_pol);
    dist->max_pol = math_deg2rad(dist->max_pol);

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpa", &dist->min_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpa", &dist->max_vpara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vpa", &dist->n_vpara,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vpe", &dist->min_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vpe", &dist->max_vperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vpe", &dist->n_vperp,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_t", &dist->n_time,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    int charge;
    if( hdf5_read_int(OPTPATH "DIST_MIN_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_MAX_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->max_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_NBIN_q", &dist->n_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Helper function to read dist5D settings from HDF5 file
 *
 * @param f the file where settings are read
 * @param dist pointer to distrho6D diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_distrho6D(hid_t file, dist_rho6D_offload_data* dist,
                                char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    if( hdf5_read_double(OPTPATH "DIST_MIN_rho", &dist->min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_rho", &dist->max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_rho", &dist->n_rho,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_phi", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_phi", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_phi", &dist->n_phi,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_pol", &dist->min_pol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_pol", &dist->max_pol,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_pol", &dist->n_pol,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_pol = math_deg2rad(dist->min_pol);
    dist->max_pol = math_deg2rad(dist->max_pol);

    if( hdf5_read_double(OPTPATH "DIST_MIN_vR", &dist->min_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vR", &dist->max_vr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vR", &dist->n_vr,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vphi", &dist->min_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vphi", &dist->max_vphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vphi", &dist->n_vphi,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_vz", &dist->min_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_vz", &dist->max_vz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_vz", &dist->n_vz,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    if( hdf5_read_double(OPTPATH "DIST_MIN_t", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_t", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "DIST_NBIN_t", &dist->n_time,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    int charge;
    if( hdf5_read_int(OPTPATH "DIST_MIN_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->min_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_MAX_q", &charge,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->max_q = (real)charge;
    if( hdf5_read_int(OPTPATH "DIST_NBIN_q", &dist->n_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    return 0;
}

/**
 * @brief Helper function to read orbit diagnostics settings from HDF5 file
 *
 * @param f the file where settings are read
 * @param dist pointer to orbit diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_orbits(hid_t file, diag_orb_offload_data* orbits,
                             char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt-XXXXXXXXXX/"

    if( hdf5_read_int(OPTPATH "ORBITWRITE_MODE",
                      &(orbits->mode),
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ORBITWRITE_NTOROIDALPLOTS",
                      &(orbits->ntoroidalplots),
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ORBITWRITE_TOROIDALANGLES",
                         (orbits->toroidalangles),
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ORBITWRITE_NPOLOIDALPLOTS",
                      &(orbits->npoloidalplots),
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ORBITWRITE_POLOIDALANGLES",
                         (orbits->poloidalangles),
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ORBITWRITE_INTERVAL",
                         &(orbits->writeInterval),
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(OPTPATH "ORBITWRITE_LASTNPOINTS",
                      &(orbits->writeNlast),
                      file, qid, __FILE__, __LINE__) ) {return 1;}

    for(int i=0; i<orbits->ntoroidalplots; i++) {
        orbits->toroidalangles[i] = orbits->toroidalangles[i]*CONST_PI/180;
    }
    for(int i=0; i<orbits->npoloidalplots; i++) {
        orbits->poloidalangles[i] = orbits->poloidalangles[i]*CONST_PI/180;
    }

    return 0;
}
