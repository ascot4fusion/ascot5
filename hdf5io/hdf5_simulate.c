#include <stdlib.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../diag.h"
#include "../diag_orb.h"
#include "../dist_5D.h"
#include "../dist_6D.h"
#include "../endcond.h"
#include "../math.h"
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
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_ADAPTIVE", active, path), &sim->enable_ada);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/RECORD_GO_AS_GC", active, path), &sim->record_GOasGC);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/FIXEDSTEP_USE_USERDEFINED", active, path), &sim->fix_usrdef_use);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/FIXEDSTEP_USERDEFINED", active, path), &sim->fix_usrdef_val);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/FIXEDSTEP_NSTEPS_PER_GYROTIME", active, path), &sim->fix_stepsPerGO);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_TOL_ORBIT", active, path), &sim->ada_tol_orbfol);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_TOL_CCOL", active, path), &sim->ada_tol_clmbcol);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_MAX_DRHO", active, path), &sim->ada_max_drho);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_MAX_DPHI", active, path), &sim->ada_max_dphi);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ADAPTIVE_MAX_ACC", active, path), &sim->ada_max_acc);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_ORBIT_FOLLOWING", active, path), &sim->enable_orbfol);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_COULOMB_COLLISIONS", active, path), &sim->enable_clmbcol);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    /* End conditions */
    int ec;
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_SIMTIMELIM", active, path), &ec);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    sim->endcond_active = 0;
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_tmax;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_SIM_TIME", active, path), &sim->endcond_maxSimTime);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_CPUTIMELIM", active, path), &ec);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_cpumax;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_CPU_TIME", active, path), &sim->endcond_maxCpuTime);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_RHOLIM", active, path), &ec);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_rhomax | endcond_rhomin;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_RHO", active, path), &sim->endcond_maxRho);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MIN_RHO", active, path), &sim->endcond_minRho);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_ENERGYLIM", active, path), &ec);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_emin | endcond_therm;
    }
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MIN_ENERGY", active, path), &sim->endcond_minEkin);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    sim->endcond_minEkin = sim->endcond_minEkin*CONST_E;
    err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MIN_ENERGY_TIMES_THERMAL", active, path), &sim->endcond_minEkinPerTe);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_WALLHIT", active, path), &ec);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_wall;
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAXORBS", active, path), &ec);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(ec){
        sim->endcond_active = sim->endcond_active | endcond_polmax | endcond_tormax;
    }
    int temp[1];
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_POLOIDALORBS", active, path), temp);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    sim->endcond_maxPolOrb = temp[0] * 2 *CONST_PI;
    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENDCOND_MAX_TOROIDALORBS", active, path), temp);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    sim->endcond_maxTorOrb = temp[0] * 2 *CONST_PI;

    /* Diagnostics */
    diag_offload_data* diag = &sim->diag_offload_data;

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_R_phi_z_vpa_vpe_t_q_DIST", active, path), &diag->dist5D_collect);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(diag->dist5D_collect) {
        dist_5D_offload_data* dist = &diag->dist5D;

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_R", active, path), &dist->min_r);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_R", active, path), &dist->max_r);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_R", active, path), &dist->n_r);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_phi", active, path), &dist->min_phi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_phi", active, path), &dist->max_phi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_phi", active, path), &dist->n_phi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        dist->min_phi = math_deg2rad(dist->min_phi);
        dist->max_phi = math_deg2rad(dist->max_phi);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_z", active, path), &dist->min_z);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_z", active, path), &dist->max_z);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_z", active, path), &dist->n_z);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_vpa", active, path), &dist->min_vpara);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_vpa", active, path), &dist->max_vpara);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_vpa", active, path), &dist->n_vpara);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_vpe", active, path), &dist->min_vperp);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_vpe", active, path), &dist->max_vperp);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_vpe", active, path), &dist->n_vperp);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_t", active, path), &dist->min_time);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_t", active, path), &dist->max_time);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_t", active, path), &dist->n_time);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

	int charge;
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_q", active, path), &charge);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
	dist->min_q = (real)charge;
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_q", active, path), &charge);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
	dist->max_q = (real)charge;
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_q", active, path), &dist->n_q);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_R_phi_z_vR_vphi_vz_t_q_DIST", active, path), &diag->dist6D_collect);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(diag->dist6D_collect) {
        dist_6D_offload_data* dist = &diag->dist6D;

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_R", active, path), &dist->min_r);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_R", active, path), &dist->max_r);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_R", active, path), &dist->n_r);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_phi", active, path), &dist->min_phi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_phi", active, path), &dist->max_phi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_phi", active, path), &dist->n_phi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        dist->min_phi = math_deg2rad(dist->min_phi);
        dist->max_phi = math_deg2rad(dist->max_phi);

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_z", active, path), &dist->min_z);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_z", active, path), &dist->max_z);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_z", active, path), &dist->n_z);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_vR", active, path), &dist->min_vr);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_vR", active, path), &dist->max_vr);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_vR", active, path), &dist->n_vr);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_vphi", active, path), &dist->min_vphi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_vphi", active, path), &dist->max_vphi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_vphi", active, path), &dist->n_vphi);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_vz", active, path), &dist->min_vz);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_vz", active, path), &dist->max_vz);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_vz", active, path), &dist->n_vz);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_t", active, path), &dist->min_time);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_t", active, path), &dist->max_time);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_t", active, path), &dist->n_time);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

	int charge;
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MIN_q", active, path), &charge);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
	dist->min_q = (real)charge;
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_MAX_q", active, path), &charge);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
	dist->max_q = (real)charge;
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/DIST_NBIN_q", active, path), &dist->n_q);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_ORBITWRITE", active, path), &diag->orb_collect);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
    if(diag->orb_collect) {
        diag_orb_offload_data* orbits = &diag->orbits;

        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_MODE", active, path), &orbits->mode);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_NTOROIDALPLOTS", active, path), &orbits->ntoroidalplots);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_TOROIDALANGLES", active, path), orbits->toroidalangles);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_NPOLOIDALPLOTS", active, path), &orbits->npoloidalplots);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_POLOIDALANGLES", active, path), orbits->poloidalangles);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_double(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_INTERVAL", active, path), &orbits->writeInterval);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}
        err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ORBITWRITE_LASTNPOINTS", active, path), &orbits->writeNlast);
	if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

        for(int i=0; i<orbits->ntoroidalplots; i++) {
            orbits->toroidalangles[i] = orbits->toroidalangles[i]*CONST_PI/180;
        }
        for(int i=0; i<orbits->npoloidalplots; i++) {
            orbits->poloidalangles[i] = orbits->poloidalangles[i]*CONST_PI/180;
        }
    }

    err = H5LTread_dataset_int(f, hdf5_generate_qid_path("/options/opt-XXXXXXXXXX/ENABLE_DEBUGDIST", active, path), &diag->debug_collect);
    if(err) {printf("Error while reading HDF5 data at %s line %d", __FILE__, __LINE__); return -1;}

    return 1;
}
