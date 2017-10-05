#include <stdlib.h>
#include <string.h>
#include "endcond.h"
#include "hdf5io/hdf5_orbits.h"
#include "offload.h"
#include "particle.h"
#include "simulate.h"
#include "simulate/simulate_ml_adaptive.h"
#include "simulate/simulate_gc_adaptive.h"
#include "simulate/simulate_gc_fixed.h"
#include "simulate/simulate_fo_fixed.h"

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
    particle_queue pq_hybrid;

    pq.n = 0;
    pq_hybrid.n = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].endcond == 0) {
            pq.n++;
        }
        else if(p[i].endcond == endcond_hybrid) {
            pq_hybrid.n++;
        }
    }

    pq.p = (particle_state**) malloc(pq.n * sizeof(particle_state*));
    pq_hybrid.p = (particle_state**)malloc(pq_hybrid.n*sizeof(particle_state*));

    pq.next = 0;
    pq_hybrid.next = 0;
    for(int i = 0; i < n_particles; i++) {
        if(p[i].endcond == 0) {
            pq.p[pq.next++] = &p[i];
        }
        else if(p[i].endcond == endcond_hybrid) {
            pq_hybrid.p[pq_hybrid.next++] = &p[i];
        }

    }
    pq.next = 0;

    #if VERBOSE >= 1
    printf("All fields initialized. Simulation begins.\n");
    #endif

    if(pq.n > 0 && (sim.sim_mode == simulate_mode_gc
            || sim.sim_mode == simulate_mode_hybrid)) {
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
    else if(pq.n > 0 && sim.sim_mode == simulate_mode_fo) {
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
    else if(pq.n > 0 && sim.sim_mode == simulate_mode_ml) {
        sim.diag_data.orbits.type = diag_orb_type_ml;
        #pragma omp parallel
        {
            simulate_ml_adaptive(&pq, &sim);
        }
    }

    /* Finish simulating hybrid particles with fo */
    if(sim.sim_mode == simulate_mode_hybrid) {

        /* Determine the number markers that should be run 
	 * in fo after previous gc simulation */
        int n_new = 0;
        for(int i = 0; i < pq.n; i++) {
            if(pq.p[i]->endcond == endcond_hybrid) {
                n_new++;
            }
        }

        if(n_new > 0) {
	    /* Reallocate and add "old" hybrid particles to the hybrid queue */
            particle_state** tmp = pq_hybrid.p;
            pq_hybrid.p = (particle_state**) malloc((pq_hybrid.n + n_new)
                    * sizeof(particle_state*));
            memcpy(pq_hybrid.p, tmp, pq_hybrid.n * sizeof(particle_state*));
            free(tmp);

	    /* Add "new" hybrid particles and reset their end condition */
	    pq_hybrid.n += n_new;
            for(int i = 0; i < pq.n; i++) {
                if(pq.p[i]->endcond == endcond_hybrid) {
		    pq.p[i]->endcond = 0;
		    pq_hybrid.p[pq_hybrid.next++] = pq.p[i];
                }
            }
        }
        pq_hybrid.next = 0;

	sim.record_GOasGC = 1; // Make sure we don't collect fos in gc diagnostics
	sim.diag_data.orbits.type = diag_orb_type_gc;
        #pragma omp parallel
        {
            simulate_fo_fixed(&pq_hybrid, &sim);
        }
    }

    free(pq.p);
    free(pq_hybrid.p);

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
