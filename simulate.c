#include <stdlib.h>
#include "particle.h"
#include "simulate.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"

void simulate(int id, int n_particles, input_particle* p,
              sim_offload_data* offload_data,
              real* B_offload_array,
              real* E_offload_array,
              real* plasma_offload_array,
              real* wall_offload_array,
              real* diag_offload_array) {
    sim_data sim;

    
}

void sim_init(sim_data* sim, sim_offload_data* offload_data) {
    sim->sim_mode             = offload_data->sim_mode;
    sim->enable_ada           = offload_data->enable_ada;
    sim->record_GOasGC        = offload_data->record_GOasGC;

    sim->fix_usrdef_use       = offload_data->fix_usrdef_use;
    sim->fix_usrdef_val       = offload_data->fix_usrdef_val;
    sim->fix_stepsPerGO       = offload_data->fix_stepsPerGO;

    sim->ada_tol_orbfol       = offload_data->ada_tol_orbfol;
    sim->ada_tol_clmbcol      = offload_data->ada_tol_clmbcol;
    sim->ada_max_drho         = offload_data->ada_max_drho;
    sim->ada_max_dphi         = offload_data->ada_max_dphi;
    sim->ada_max_acc          = offload_data->ada_max_acc;

    sim->enable_orbfol        = offload_data->enable_orbfol;
    sim->enable_clmbcol       = offload_data->enable_clmbcol;

    sim->endcond_active       = offload_data->endcond_active;
    sim->endcond_maxSimTime   = offload_data->endcond_maxSimTime;
    sim->endcond_maxCpuTime   = offload_data->endcond_maxCpuTime;
    sim->endcond_minRho       = offload_data->endcond_minRho;
    sim->endcond_maxRho       = offload_data->endcond_maxRho;
    sim->endcond_minEkin      = offload_data->endcond_minEkin;
    sim->endcond_minEkinPerTe = offload_data->endcond_minEkinPerTe;
    sim->endcond_maxTorOrb    = offload_data->endcond_maxTorOrb;
    sim->endcond_maxPolOrb    = offload_data->endcond_maxPolOrb;
}

void simulate_begin(int id, int n_particles, input_particle* p,
		    sim_offload_data* offload_data, sim_data* sim,
		    real* B_offload_array,
		    real* E_offload_array,
		    real* plasma_offload_array,
		    real* wall_offload_array,
		    real* diag_offload_array) {

    
    sim_init(sim, offload_data);
    wall_init(&sim->wall_data, &offload_data->wall_offload_data,
              wall_offload_array);
    B_field_init(&sim->B_data, &offload_data->B_offload_data, B_offload_array);
    E_field_init(&sim->E_data, &offload_data->E_offload_data, E_offload_array);
    plasma_1d_init(&sim->plasma_data, &offload_data->plasma_offload_data,
                   plasma_offload_array);
    diag_init(&sim->diag_data, &offload_data->diag_offload_data,
              diag_offload_array);
    
    if(sim->sim_mode == 1) {
	/* FO simulation */
	for(int i = 0; i < n_particles; i++) {
	    particle_marker_to_state(p, i, &sim->B_data, 1);
	}
    }
    else if(sim->sim_mode == 2 || sim->sim_mode == 3) {
	/* GC or hybrid simulation */
	for(int i = 0; i < n_particles; i++) {
	    particle_marker_to_state(p, i, &sim->B_data, 2);
	}
    }
    else if(sim->sim_mode == 4) {
	/* ML simulation */
	for(int i = 0; i < n_particles; i++) {
	    particle_marker_to_state(p, i, &sim->B_data, 4);
	}
    }

}


void simulate_continue(int id, int n_particles, input_particle* p,
		       sim_data* sim) {

    particle_queue p_fo;
    particle_queue p_gc;
    particle_queue p_ml;

    p_fo.n = 0;
    p_gc.n = 0;
    p_ml.n = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].type == input_particle_type_ps) {
            p_fo.n++;
        } else if(p[i].type == input_particle_type_gcs) {
            p_gc.n++;
        } else if(p[i].type == input_particle_type_mls) {
            p_ml.n++;
        }
    }

    p_fo.p = (particle_state**) malloc(p_fo.n * sizeof(particle_state*));
    p_gc.p = (particle_state**) malloc(p_gc.n * sizeof(particle_state*));
    p_ml.p = (particle_state**) malloc(p_ml.n * sizeof(particle_state*));

    p_fo.next = 0;
    p_gc.next = 0;
    p_ml.next = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].type == input_particle_type_ps) {
            p_fo.p[p_fo.next++] = &p[i].p_s;
        } else if(p[i].type == input_particle_type_gcs) {
            p_gc.p[p_gc.next++] = &p[i].p_s;
        } else if(p[i].type == input_particle_type_mls) {
            p_ml.p[p_ml.next++] = &p[i].p_s;
        }
    }

    p_fo.next = 0;
    p_gc.next = 0;
    p_ml.next = 0;
    

    if(p_gc.n > 0) {
	if(sim->enable_ada) {
	    #pragma omp parallel
	    {
	        simulate_gc_adaptive(&p_gc, sim);
	    }
	}
	else {
	    #pragma omp parallel
	    {
	        simulate_gc_fixed(&p_gc, sim);
	    }
	}
    }
    else if(p_fo.n > 0) {
	#pragma omp parallel
	{
	    simulate_fo_fixed(&p_fo, sim);
	}
    }
    else if(p_ml.n > 0) {
	#pragma omp parallel
	{
	    simulate_ml_adaptive(&p_ml, sim);
	}
    }

    free(p_fo.p);
    free(p_gc.p);
    free(p_ml.p);

}

void simulate_hybrid(int id, int n_particles, input_particle* p,
		     sim_data* sim) {

}

void simulate_end(sim_data* sim) {
    diag_clean(&sim->diag_data);
}

