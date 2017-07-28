#include <stdlib.h>
#include "particle.h"
#include "simulate.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"
#include "hdf5io/hdf5_orbits.h"
#include "offload.h"

void simulate(int id, int n_particles, particle_state* p,
              sim_offload_data* sim_offload,
              offload_package* offload_data,
              real* offload_array,
              real* diag_offload_array) {
    sim_data sim;
    sim_init(&sim, sim_offload);

    real* ptr;
    ptr = offload_unpack(offload_data, offload_array,
                         sim_offload->B_offload_data.offload_array_length);
    B_field_init(&sim.B_data, &sim_offload->B_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
                         sim_offload->E_offload_data.offload_array_length);
    E_field_init(&sim.E_data, &sim_offload->E_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
                         sim_offload->plasma_offload_data.offload_array_length);
    plasma_1d_init(&sim.plasma_data, &sim_offload->plasma_offload_data, ptr);

    ptr = offload_unpack(offload_data, offload_array,
                         sim_offload->wall_offload_data.offload_array_length);
    wall_init(&sim.wall_data, &sim_offload->wall_offload_data, ptr);

    diag_init(&sim.diag_data, &sim_offload->diag_offload_data,
              diag_offload_array);
    
    particle_queue pq;
    particle_queue pqhyb;

    pq.n = 0;
    pqhyb.n = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].endcond == 0) {
            pq.n++;
        } 
	/* TODO check here for the hybrid endcond
	else if() {
            pqhyb.n++;
        } 
	*/
    }

    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pqhyb.p = (particle_state**) malloc(pqhyb.n * sizeof(particle_state*));

    pq.next = 0;
    pqhyb.next = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].endcond == 0) {
            pq.p[pq.next++] = &p[i];
        } 
	/* hybrid
	else if() {
	    p[i].p_s.endcond = 0;
            pq.p[pqhyb.next++] = &p[i];
        } 
	*/
    }

    pq.next = 0;
    pqhyb.next = 0;

    if(pq.n > 0 && (sim.sim_mode == 2 || sim.sim_mode == 3) ) {
	sim.diag_data.orbits.type = diag_orb_type_gc;

	if(sim.enable_ada) {
	    #pragma omp parallel
	    {
	        simulate_gc_adaptive(&pq, &sim);
	    }
	}
	else {
	    #pragma omp parallel
	    {
	        simulate_gc_fixed(&pq, &sim);
	    }
	}
    }
    else if(pq.n > 0 && sim.sim_mode == 1) {
	if(sim.record_GOasGC) {
	    sim.diag_data.orbits.type = diag_orb_type_gc;
	}
	else {
	    sim.diag_data.orbits.type = diag_orb_type_fo;
	}
	#pragma omp parallel
	{
	    simulate_fo_fixed(&pq, &sim);
	}
    }
    else if(pq.n > 0 && sim.sim_mode == 4) {
	sim.diag_data.orbits.type = diag_orb_type_ml;
	#pragma omp parallel
	{
	    simulate_ml_adaptive(&pq, &sim);
	}
    }
    else if(pqhyb.n > 0 && sim.sim_mode == 3) {
        /* fo simulation for the hybrid mode */
	#pragma omp parallel
	{
	    simulate_fo_fixed(&pqhyb, &sim);
	}
    }

    free(pq.p);
    free(pqhyb.p);

    // Temporary solution
    #ifdef NOTARGET
        hdf5_orbits_write(&sim, sim_offload->hdf5_out);
    #endif

    diag_clean(&sim.diag_data);
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

void simulate_continue(int id, int n_particles, input_particle* p,
		       sim_data* sim) {

}

void simulate_hybrid(int id, int n_particles, input_particle* p,
		     sim_data* sim) {

}

