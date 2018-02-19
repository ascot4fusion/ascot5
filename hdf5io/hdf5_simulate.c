#include <stdlib.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../diag.h"
#include "../diag_orb.h"
#include "../dist_5D.h"
#include "../dist_6D.h"
#include "../endcond.h"
#include "../simulate.h"
#include "hdf5.h"
#include "hdf5_helpers.h"
#include "hdf5_hl.h"


int hdf5_simulate(hid_t f, sim_offload_data* sim){
    herr_t err;

    #if VERBOSE > 0
        printf("Reading options from the HDF5 file...\n");
    #endif

    err = hdf5_find_group(f, "/options/");
    if(err < 0) {
        return -1;
    }

    char active[11];
    err = H5LTget_attribute_string(f, "/options/", "active", active);
    if(err < 0) {
        return -1;
    }
    active[10] = '\0';

    #if VERBOSE > 0
        printf("Active qid is %s\n", active);
    #endif

    char path[256];

    /* End conditions */
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/SIM_MODE", active, path), &sim->sim_mode);
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_ADAPTIVE", active, path), &sim->enable_ada);
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/RECORD_GO_AS_GC", active, path), &sim->record_GOasGC);

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/FIXEDSTEP_USE_USERDEFINED", active, path), &sim->fix_usrdef_use);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/FIXEDSTEP_USERDEFINED", active, path), &sim->fix_usrdef_val);
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/FIXEDSTEP_NSTEPS_PER_GYROTIME", active, path), &sim->fix_stepsPerGO);

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_TOL_ORBIT", active, path), &sim->ada_tol_orbfol);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_TOL_CCOL", active, path), &sim->ada_tol_clmbcol);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_MAX_DRHO", active, path), &sim->ada_max_drho);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_MAX_DPHI", active, path), &sim->ada_max_dphi);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_MAX_ACC", active, path), &sim->ada_max_acc);

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_ORBIT_FOLLOWING", active, path), &sim->enable_orbfol);
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_COULOMB_COLLISIONS", active, path), &sim->enable_clmbcol);

    /* End conditions */
    int ec;
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_SIMTIMELIM", active, path), &ec);
    sim->endcond_active = 0;
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_tmax;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_SIM_TIME", active, path), &sim->endcond_maxSimTime);

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_CPUTIMELIM", active, path), &ec);
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_cpumax;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_CPU_TIME", active, path), &sim->endcond_maxCpuTime);

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_RHOLIM", active, path), &ec);
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_rhomax | endcond_rhomin;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_RHO", active, path), &sim->endcond_maxRho);
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MIN_RHO", active, path), &sim->endcond_minRho);

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_ENERGYLIM", active, path), &ec);
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_emin | endcond_therm;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MIN_ENERGY", active, path), &sim->endcond_minEkin);
    sim->endcond_minEkin = sim->endcond_minEkin*CONST_E;
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MIN_ENERGY_TIMES_THERMAL", active, path), &sim->endcond_minEkinPerTe);

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_WALLHIT", active, path), &ec);
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_wall;
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAXORBS", active, path), &ec);
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_polmax | endcond_tormax;
    }
    int temp[1];
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_POLOIDALORBS", active, path), temp);
    sim->endcond_maxPolOrb = temp[0] * 2 *CONST_PI;
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_TOROIDALORBS", active, path), temp);
    sim->endcond_maxTorOrb = temp[0] * 2 *CONST_PI;

    /* Diagnostics */
    diag_offload_data* diag = &sim->diag_offload_data;

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_R_phi_z_vpa_vpe_t_q_DIST", active, path), &diag->dist5D_collect);
    if(diag->dist5D_collect) {
        dist_5D_offload_data* dist = &diag->dist5D;

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_R", active, path), &dist->min_r);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_R", active, path), &dist->max_r);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_R", active, path), &dist->n_r);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_phi", active, path), &dist->min_phi);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_phi", active, path), &dist->max_phi);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_phi", active, path), &dist->n_phi);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_z", active, path), &dist->min_z);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_z", active, path), &dist->max_z);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_z", active, path), &dist->n_z);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_vpa", active, path), &dist->min_vpara);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_vpa", active, path), &dist->max_vpara);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_vpa", active, path), &dist->n_vpara);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_vpe", active, path), &dist->min_vperp);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_vpe", active, path), &dist->max_vperp);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_vpe", active, path), &dist->n_vperp);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_t", active, path), &dist->min_time);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_t", active, path), &dist->max_time);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_t", active, path), &dist->n_time);

        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MIN_q", active, path), &dist->min_q);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_MAX_q", active, path), &dist->max_q);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vpa_vpe_t_q_NBIN_q", active, path), &dist->n_q);
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST", active, path), &diag->dist6D_collect);
    if(diag->dist6D_collect) {
        dist_6D_offload_data* dist = &diag->dist6D;

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_R", active, path), &dist->min_r);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_R", active, path), &dist->max_r);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_R", active, path), &dist->n_r);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_phi", active, path), &dist->min_phi);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_phi", active, path), &dist->max_phi);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_phi", active, path), &dist->n_phi);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_z", active, path), &dist->min_z);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_z", active, path), &dist->max_z);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_z", active, path), &dist->n_z);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_vR", active, path), &dist->min_vr);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_vR", active, path), &dist->max_vr);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_vR", active, path), &dist->n_vr);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_vphi", active, path), &dist->min_vphi);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_vphi", active, path), &dist->max_vphi);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_vphi", active, path), &dist->n_vphi);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_vz", active, path), &dist->min_vz);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_vz", active, path), &dist->max_vz);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_vz", active, path), &dist->n_vz);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_t", active, path), &dist->min_time);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_t", active, path), &dist->max_time);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_t", active, path), &dist->n_time);

        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MIN_q", active, path), &dist->min_q);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_MAX_q", active, path), &dist->max_q);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_R_phi_z_vR_vphi_vz_t_q_NBIN_q", active, path), &dist->n_q);
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_ORBITWRITE", active, path), &diag->orb_collect);
    if(diag->orb_collect) {
        diag_orb_offload_data* orbits = &diag->orbits;

        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_MODE", active, path), &orbits->mode);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_NTOROIDALPLOTS", active, path), &orbits->ntoroidalplots);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_TOROIDALANGLES", active, path), orbits->toroidalangles);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_NPOLOIDALPLOTS", active, path), &orbits->npoloidalplots);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_POLOIDALANGLES", active, path), orbits->poloidalangles);
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_INTERVAL", active, path), &orbits->writeInterval);
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_LASTNPOINTS", active, path), &orbits->writeNlast);

        for(int i=0; i<orbits->ntoroidalplots; i++) {
            orbits->toroidalangles[i] = orbits->toroidalangles[i]*CONST_PI/180;
        }
        for(int i=0; i<orbits->npoloidalplots; i++) {
            orbits->poloidalangles[i] = orbits->poloidalangles[i]*CONST_PI/180;
        }
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_DEBUGDIST", active, path), &diag->debug_collect);

    return 1;
}
