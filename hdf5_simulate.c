#include <stdlib.h>
#include "ascot5.h"
#include "simulate.h"
#include "endcond.h"
#include "distributions.h"
#include "diag.h"
#include "diag_orb.h"
#include "hdf5.h"
#include "hdf5_helpers.h"
#include "hdf5_hl.h"


void hdf5_simulate(hid_t f, sim_offload_data* sim){
    herr_t err;	
    
    /* End conditions */
    err = H5LTread_dataset_int(f, "/options/SIM_MODE", &sim->sim_mode);
    err = H5LTread_dataset_int(f, "/options/ENABLE_ADAPTIVE", &sim->enable_ada);
    err = H5LTread_dataset_int(f, "/options/RECORD_GO_AS_GC", &sim->record_GOasGC);
	
    err = H5LTread_dataset_int(f, "/options/FIXEDSTEP_USE_USERDEFINED", &sim->fix_usrdef_use);
    err = H5LTread_dataset_double(f, "/options/FIXEDSTEP_USERDEFINED", &sim->fix_usrdef_val);
    err = H5LTread_dataset_int(f, "/options/FIXEDSTEP_NSTEPS_PER_GYROTIME", &sim->fix_stepsPerGO);
	
    err = H5LTread_dataset_double(f, "/options/ADAPTIVE_TOL_ORBIT", &sim->ada_tol_orbfol);
    err = H5LTread_dataset_double(f, "/options/ADAPTIVE_TOL_CCOL", &sim->ada_tol_clmbcol);
    err = H5LTread_dataset_double(f, "/options/ADAPTIVE_MAX_DRHO", &sim->ada_max_drho);
    err = H5LTread_dataset_double(f, "/options/ADAPTIVE_MAX_DPHI", &sim->ada_max_dphi);
    err = H5LTread_dataset_double(f, "/options/ADAPTIVE_MAX_ACC", &sim->ada_max_acc);
	
    err = H5LTread_dataset_int(f, "/options/ENABLE_ORBIT_FOLLOWING", &sim->enable_orbfol);
    err = H5LTread_dataset_int(f, "/options/ENABLE_COULOMB_COLLISIONS", &sim->enable_clmbcol);

    

    /* End conditions */
    int ec;
    err = H5LTread_dataset_int(f, "/options/ENDCOND_SIMTIMELIM", &ec);
    if(ec){
	sim->endcond_active = sim->endcond_active | endcond_tmax;
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MAX_SIM_TIME", &sim->endcond_maxSimTime);
    }
	
    err = H5LTread_dataset_int(f, "/options/ENDCOND_CPUTIMELIM", &ec);
    if(ec){
	sim->endcond_active = sim->endcond_active | endcond_cpumax;
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MAX_CPU_TIME", &sim->endcond_maxCpuTime);
    }
	
    err = H5LTread_dataset_int(f, "/options/ENDCOND_RHOLIM", &ec);
    if(ec){
	sim->endcond_active = sim->endcond_active | endcond_rhomax | endcond_rhomin;
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MAX_RHO", &sim->endcond_maxRho);
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MIN_RHO", &sim->endcond_minRho);
    }
	
    err = H5LTread_dataset_int(f, "/options/ENDCOND_ENERGYLIM", &ec);
    if(ec){
	sim->endcond_active = sim->endcond_active | endcond_emin | endcond_therm;
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MIN_ENERGY", &sim->endcond_minEkin);
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MIN_ENERGY_TIMES_THERMAL", &sim->endcond_minEkinPerTe);
    }
	
    err = H5LTread_dataset_int(f, "/options/ENDCOND_WALLHIT", &ec);
    if(ec){
	sim->endcond_active = sim->endcond_active | endcond_wall;
    }
	
    err = H5LTread_dataset_int(f, "/options/ENDCOND_MAXORBS", &ec);
    if(ec){
	sim->endcond_active = sim->endcond_active | endcond_polmax | endcond_tormax;
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MAX_POLOIDALORBS", &sim->endcond_maxTorOrb);
	err = H5LTread_dataset_double(f, "/options/ENDCOND_MAX_TOROIDALORBS", &sim->endcond_maxPolOrb);
    }

    /* Diagnostics */
    diag_offload_data* diag = &sim->diag_offload_data;
    int enable;
    
    err = H5LTread_dataset_int(f, "/options/ENABLE_RZVparaVperp_DIST", &diag->dist4D_collect);
    if(diag->dist4D_collect) {
	dist_rzvv_offload_data* dist_rzvpavpe = &diag->dist4D;

	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MIN_R", &dist_rzvpavpe->min_r);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MAX_R", &dist_rzvpavpe->max_r);
	err = H5LTread_dataset_int(f, "/options/DIST_RZVparaVperp_BIN_R", &dist_rzvpavpe->n_r);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MIN_Z", &dist_rzvpavpe->min_z);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MAX_Z", &dist_rzvpavpe->max_z);
	err = H5LTread_dataset_int(f, "/options/DIST_RZVparaVperp_BIN_Z", &dist_rzvpavpe->n_z);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MIN_VPARA", &dist_rzvpavpe->min_vpara);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MAX_VPARA", &dist_rzvpavpe->max_vpara);
	err = H5LTread_dataset_int(f, "/options/DIST_RZVparaVperp_BIN_VPARA", &dist_rzvpavpe->n_vpara);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MIN_VPERP", &dist_rzvpavpe->min_vperp);
	err = H5LTread_dataset_double(f, "/options/DIST_RZVparaVperp_MAX_VPERP", &dist_rzvpavpe->max_vperp);
	err = H5LTread_dataset_int(f, "/options/DIST_RZVparaVperp_BIN_VPERP", &dist_rzvpavpe->n_vperp);
    }
    
    err = H5LTread_dataset_int(f, "/options/ENABLE_ORBITWRITE", &diag->orb_collect);
    if(diag->orb_collect) {
	diag_orb_data* orbits = &diag->orbits;
    }

/*
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_MODE", &sim->orbitwrite_mode);
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_NTOROIDALPLOTS", &sim->orbitwrite_ntoroidalplots);
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_TOROIDALANGLES", &sim->orbitwrite_toroidalangles);
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_NPOLOIDALPLOTS", &sim->orbitwrite_npoloidalplots);
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_POLOIDALANGLES", &sim->orbitwrite_poloidalangles);
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_INTERVAL", &sim->orbitwrite_interval);
	err = H5LTread_dataset_double(f, "/options/ORBITWRITE_LASTNPOINTS", &sim->orbitwrite_lastnpoints);
	err = H5LTread_dataset_double(f, "/options/ENABLE_DEBUGDIST", &sim->enable_debugdist);
    */
}
