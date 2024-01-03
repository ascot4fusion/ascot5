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
#include "../diag/dist_com.h"
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
int hdf5_options_read_distCOM(hid_t file, dist_COM_offload_data* dist,
                              char* qid);
int hdf5_options_read_diagorb(hid_t file, diag_orb_offload_data* diagorb,
                              char* qid);
int hdf5_options_read_diagtrcof(hid_t file,
                                diag_transcoef_offload_data* diagtrcof,
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
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat; // For reading float data that is converted to int.

    if( hdf5_read_double(OPTPATH "SIM_MODE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->sim_mode = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_ADAPTIVE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->enable_ada = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "RECORD_MODE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->record_mode = (int)tempfloat;


    if( hdf5_read_double(OPTPATH "FIXEDSTEP_USE_USERDEFINED", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->fix_usrdef_use = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "FIXEDSTEP_USERDEFINED", &sim->fix_usrdef_val,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "FIXEDSTEP_GYRODEFINED", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->fix_gyrodef_nstep = (int)tempfloat;


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
    if( hdf5_read_double(OPTPATH "ENABLE_MHD", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->enable_mhd = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_ATOMIC", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->enable_atomic = (int)tempfloat;
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
    if( hdf5_read_double(OPTPATH "REVERSE_TIME", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->reverse_time = (int)tempfloat;

    int ec;
    sim->endcond_active = 0;

    if( hdf5_read_double(OPTPATH "ENDCOND_SIMTIMELIM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_tlim * (ec > 0);
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
    if( hdf5_read_double(OPTPATH "ENDCOND_NEUTRALIZED", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_neutr * (ec > 0);
    if( hdf5_read_double(OPTPATH "ENDCOND_IONIZED", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    ec = (int)tempfloat;
    sim->endcond_active = sim->endcond_active | endcond_ioniz * (ec > 0);

    sim->endcond_torandpol = 0;
    if( ec == 2) {
        sim->endcond_torandpol = 1;
    }

    if( hdf5_read_double(OPTPATH "ENDCOND_LIM_SIMTIME",
                         &sim->endcond_lim_simtime,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_MILEAGE",
                         &sim->endcond_max_mileage,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_CPUTIME",
                         &sim->endcond_max_cputime,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_RHO",
                         &sim->endcond_max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MIN_RHO",
                         &sim->endcond_min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MIN_ENERGY",
                         &sim->endcond_min_ekin,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ENDCOND_MIN_THERMAL",
                         &sim->endcond_min_thermal,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->endcond_min_ekin *= CONST_E; // eV -> J

    int temp;
    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_POLOIDALORBS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    temp = (int)tempfloat;
    sim->endcond_max_polorb = temp * 2 *CONST_PI;

    if( hdf5_read_double(OPTPATH "ENDCOND_MAX_TOROIDALORBS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    temp = (int)tempfloat;
    sim->endcond_max_tororb = temp * 2 *CONST_PI;

    /* BMC options */
    if( hdf5_read_double(OPTPATH "BMC_TIMEDEPENDENT", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->bmc_timedependent = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "BMC_ORBIT_SUBCYCLES", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->bmc_orbit_subcycles = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "BMC_TIMESTEP", &sim->bmc_timestep,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "BMC_TSTART", &sim->bmc_tstart,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "BMC_TSTOP", &sim->bmc_tstop,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "BMC_MASS", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->bmc_mass = tempfloat * CONST_U;
    if( hdf5_read_double(OPTPATH "BMC_CHARGE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->bmc_charge = tempfloat * CONST_E;
    if( hdf5_read_double(OPTPATH "BMC_ANUM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->bmc_znum = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "BMC_ZNUM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    sim->bmc_znum = (int)tempfloat;

    /* See which diagnostics are active */
    diag_offload_data* diag = &sim->diag_offload_data;

    if( hdf5_read_double(OPTPATH "ENABLE_DIST_5D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->dist5D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_6D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->dist6D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_RHO5D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->distrho5D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_RHO6D", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->distrho6D_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_DIST_COM", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->distCOM_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_ORBITWRITE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->diagorb_collect = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ENABLE_TRANSCOEF", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diag->diagtrcof_collect = (int)tempfloat;

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
    if(diag->distCOM_collect) {
        if( hdf5_options_read_distCOM(file, &diag->distCOM, qid) ) {
            return 1;
        }
    }
    if(diag->diagorb_collect) {
        diag->diagorb.record_mode = sim->sim_mode;
        if(sim->record_mode && (sim->sim_mode == simulate_mode_fo ||
                                sim->sim_mode == simulate_mode_hybrid) ) {
            diag->diagorb.record_mode = simulate_mode_gc;
        }

        if( hdf5_options_read_diagorb(file, &diag->diagorb, qid) ) {
            return 1;
        }
    }

    if(diag->diagtrcof_collect) {
        if( hdf5_options_read_diagtrcof(file, &diag->diagtrcof, qid) ) {
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
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_R", &dist->min_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_R", &dist->max_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_R", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_r = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PHI", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PHI", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PHI", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_Z", &dist->min_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_Z", &dist->max_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_Z", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_z = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PPA", &dist->min_ppara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PPA", &dist->max_ppara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PPA", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_ppara = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PPE", &dist->min_pperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PPE", &dist->max_pperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PPE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pperp = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_TIME", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_TIME", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_TIME", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_CHARGE", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_CHARGE", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_CHARGE", &tempfloat,
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
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_R", &dist->min_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_R", &dist->max_r,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_R", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_r = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PHI", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PHI", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PHI", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_Z", &dist->min_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_Z", &dist->max_z,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_Z", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_z = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PR", &dist->min_pr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PR", &dist->max_pr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PR", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pr = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PPHI", &dist->min_pphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PPHI", &dist->max_pphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PPHI", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pphi = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PZ", &dist->min_pz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PZ", &dist->max_pz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PZ", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pz = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_TIME", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_TIME", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_TIME", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_CHARGE", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_CHARGE", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_CHARGE", &tempfloat,
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
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_RHO", &dist->min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_RHO", &dist->max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_RHO", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_rho = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PHI", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PHI", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PHI", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_THETA", &dist->min_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_THETA", &dist->max_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_THETA", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_theta = (int)tempfloat;
    dist->min_theta = math_deg2rad(dist->min_theta);
    dist->max_theta = math_deg2rad(dist->max_theta);

    if( hdf5_read_double(OPTPATH "DIST_MIN_PPA", &dist->min_ppara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PPA", &dist->max_ppara,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PPA", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_ppara = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PPE", &dist->min_pperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PPE", &dist->max_pperp,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PPE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pperp = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_TIME", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_TIME", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_TIME", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_CHARGE", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_CHARGE", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_CHARGE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_q = (int)tempfloat;

    return 0;
}

/**
 * @brief Helper function to read distrho6D settings from HDF5 file
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
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_RHO", &dist->min_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_RHO", &dist->max_rho,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_RHO", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_rho = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PHI", &dist->min_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PHI", &dist->max_phi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PHI", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_phi = (int)tempfloat;
    dist->min_phi = math_deg2rad(dist->min_phi);
    dist->max_phi = math_deg2rad(dist->max_phi);

    if( hdf5_read_double(OPTPATH "DIST_MIN_THETA", &dist->min_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_THETA", &dist->max_theta,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_THETA", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_theta = (int)tempfloat;
    dist->min_theta = math_deg2rad(dist->min_theta);
    dist->max_theta = math_deg2rad(dist->max_theta);

    if( hdf5_read_double(OPTPATH "DIST_MIN_PR", &dist->min_pr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PR", &dist->max_pr,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PR", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pr = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PPHI", &dist->min_pphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PPHI", &dist->max_pphi,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PPHI", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pphi = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PZ", &dist->min_pz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PZ", &dist->max_pz,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PZ", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_pz = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_TIME", &dist->min_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_TIME", &dist->max_time,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_TIME", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_time = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_CHARGE", &dist->min_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_CHARGE", &dist->max_q,
                      file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_CHARGE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_q = (int)tempfloat;

    return 0;
}

/**
 * @brief Helper function to read COM dist settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param dist pointer to distCOM diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_distCOM(hid_t file, dist_COM_offload_data* dist,
                              char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat;
    if( hdf5_read_double(OPTPATH "DIST_MIN_MU", &dist->min_mu,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_MU", &dist->max_mu,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_MU", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_mu = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_EKIN", &dist->min_Ekin,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_EKIN", &dist->max_Ekin,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_EKIN", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_Ekin = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "DIST_MIN_PTOR", &dist->min_Ptor,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_MAX_PTOR", &dist->max_Ptor,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "DIST_NBIN_PTOR", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    dist->n_Ptor = (int)tempfloat;

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
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    real tempfloat;

    if( hdf5_read_double(OPTPATH "ORBITWRITE_MODE", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diagorb->mode = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "ORBITWRITE_NPOINT", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diagorb->Npnt = (int)tempfloat;
    if( hdf5_read_double(OPTPATH "ORBITWRITE_INTERVAL",
                         &(diagorb->writeInterval),
                         file, qid, __FILE__, __LINE__) ) {return 1;}



    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        diagorb->toroidalangles[i] = TOROIDAL_ANGLE_FILLER_VALUE;
        diagorb->poloidalangles[i] = POLOIDAL_ANGLE_FILLER_VALUE;
        diagorb->radialdistances[i] = RADIAL_FILLER_VALUE;
    }

    if( hdf5_read_double(OPTPATH "ORBITWRITE_TOROIDALANGLES",
                         (diagorb->toroidalangles),
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ORBITWRITE_POLOIDALANGLES",
                         (diagorb->poloidalangles),
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(OPTPATH "ORBITWRITE_RADIALDISTANCES",
                         (diagorb->radialdistances),
                         file, qid, __FILE__, __LINE__) ) {return 1;}


    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        if(diagorb->toroidalangles[0] < 0) {
            /* Negative angle means plane is disabled */
            diagorb->ntoroidalplots = 0;
            break;
        }
        if(diagorb->toroidalangles[i] == TOROIDAL_ANGLE_FILLER_VALUE) {
            diagorb->ntoroidalplots = i;
            break;
        }
    }

    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        if(diagorb->poloidalangles[0] < 0) {
            /* Negative angle means plane is disabled */
            diagorb->npoloidalplots = 0;
            break;
        }
        if(diagorb->poloidalangles[i] == POLOIDAL_ANGLE_FILLER_VALUE) {
            diagorb->npoloidalplots = i;
            break;
        }
    }

    for(int i=0; i < DIAG_ORB_MAXPOINCARES; i++) {
        if(diagorb->radialdistances[0] < 0) {
            /* Negative angle means plane is disabled */
            diagorb->nradialplots = 0;
            break;
        }
        if(diagorb->radialdistances[i] == RADIAL_FILLER_VALUE) {
            diagorb->nradialplots = i;
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

/**
 * @brief Helper function to read transport coefficient settings from HDF5 file
 *
 * @param file the file where settings are read
 * @param diagtrcof pointer to orbit diagnostics offload data
 * @param qid QID of the options to be read
 *
 * @return zero if reading succeeded.
 */
int hdf5_options_read_diagtrcof(
    hid_t file, diag_transcoef_offload_data* diagtrcof, char* qid) {
    #undef OPTPATH
    #define OPTPATH "/options/opt_XXXXXXXXXX/"

    if( hdf5_read_double(OPTPATH "TRANSCOEF_INTERVAL", &diagtrcof->interval,
                         file, qid, __FILE__, __LINE__) ) {return 1;}

    real tempfloat;
    if( hdf5_read_double(OPTPATH "TRANSCOEF_NAVG", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diagtrcof->Navg = (int)tempfloat;

    if( hdf5_read_double(OPTPATH "TRANSCOEF_RECORDRHO", &tempfloat,
                         file, qid, __FILE__, __LINE__) ) {return 1;}
    diagtrcof->recordrho = (int)tempfloat;

    return 0;
}
