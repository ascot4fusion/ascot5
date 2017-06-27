/**
 * @file simulate_gc_fixed.c
 * @brief Simulate particles with guiding center using RK4 integrator
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>
#include "ascot5.h"
#include "step_gc_rk4.h"
#include "wall.h"
#include "distributions.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma_1d.h"
#include "simulate.h"
#include "math.h"
#include "diag.h"
#include "particle.h"
#include "endcond.h"
#include "simulate_gc_fixed.h"
#include "orbit_write.h"
#include "mccc/mccc.h"

void simulate_gc_fixed(int id, int n_particles, particle* particles,
			  sim_offload_data sim_offload,
			  real* B_offload_array,
			  real* E_offload_array,
			  real* plasma_offload_array,
			  real* wall_offload_array,
			  real* dist_offload_array) {
    sim_data sim;


/* BACKGROUND INITIALIZATION */

    /* Simulation data */
    sim_init(&sim, &sim_offload);

    /* Wall data */
    wall_init(&sim.wall_data, &sim_offload.wall_offload_data,
              wall_offload_array);

    /* Magnetic field data */
    B_field_init(&sim.B_data, &sim_offload.B_offload_data, B_offload_array);

    /* Electric field data */
    E_field_init(&sim.E_data, &sim_offload.E_offload_data, E_offload_array);

    /* Plasma data */
    plasma_1d_init(&sim.plasma_data, &sim_offload.plasma_offload_data,
                   plasma_offload_array);

    /* Diagnostics data */
    diag_init(&sim.diag_data,&sim_offload.diag_offload_data);
    dist_rzvv_init(&sim.dist_data, &sim_offload.dist_offload_data,
                   dist_offload_array);
	
   
    int i_next_prt = 0;

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
    #pragma omp parallel
    {
	int err[NSIMD];
	real hin[NSIMD];

        particle_simd_gc p;  // This array holds current states
	particle_simd_gc p0; // This array stores previous states
        int i, i_prt;

/** MARKER INITIALIZATION */
        for(i = 0; i < NSIMD; i++) {
            #pragma omp critical
            i_prt = i_next_prt++;
            if(i_prt < n_particles) {
		/* Guiding center transformation */
                particle_to_gc(&particles[i_prt], i_prt, &p, i, &sim.B_data);

		// Determine initial time step
		// TODO get this one from options
		hin[i] = 1.e-8;
            }
            else {
		/* Dummy marker to fill NSIMD when we ran out of actual particles */
                particle_to_gc_dummy(&p, i);
            }

	    /* Init dummy particles here, the (required) fields are initialized 
	     * separately at each time step */
	    particle_to_gc_dummy(&p0, i); 
        }

	
	

/* MAIN SIMULATION LOOP 
 * - Store current state
 * - Integrate motion due to bacgkround EM-field (orbit-following)
 * - Integrate scattering due to Coulomb collisions
 * - Advance time
 * - Check for end condition(s)
 * - Update diagnostics
 * - 
 * */
        int n_running = 0;
        do {
	    
            #pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		/* Store marker states */
	        p0.r[i]        = p.r[i];
		p0.phi[i]      = p.phi[i];
		p0.z[i]        = p.z[i];
		p0.vpar[i]     = p.vpar[i];
		p0.mu[i]       = p.mu[i];
		p0.theta[i]    = p.theta[i];
		p0.time[i]     = p.time[i];
		p0.running[i]  = p.running[i];
		p0.endcond[i]  = p.endcond[i];
		p0.walltile[i] = p.walltile[i];
	    }
	    
	    if(sim.enable_orbfol) {
	        step_gc_rk4(&p, hin, &sim.B_data, &sim.E_data);
	    }

	    if(sim.enable_clmbcol) {
	        mccc_step_gc_fixed(&p, &sim.B_data, &sim.plasma_data, hin, err);
            }
 
		

	    #pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		if(p.running[i]){
		    p.time[i] = p.time[i] + hin[i];
		}
	    }
	    
            endcond_check_gc(&p, &p0, &sim);

	    diag_update_gc(&sim.diag_data, &p, &p0);

            /* update number of running particles */
            n_running = 0;
            int k;
            for(k = 0; k < NSIMD; k++) {
                if( !p.running[k] && p.id[k] >= 0) {

                    gc_to_particle(&p, k, &particles[p.index[k]]);

                    #pragma omp critical
                    i_prt = i_next_prt++;
                    if(i_prt < n_particles) {
                        particle_to_gc(&particles[i_prt], i_prt, &p, k,
				       &sim.B_data);
		
			// Determine initial time step
			// TODO get this one from options
			hin[k] = 1.e-8;
						
                    }
                    else {
                        p.id[k] = -1;
                    }
	        }
            }

	    #pragma omp simd reduction(+:n_running)
	    for(k = 0; k < NSIMD; k++) {
		n_running += p.running[k];
	    }
	   

        } while(n_running > 0);
	
        

    }

    // TODO Move these to main program
    diag_write(&sim.diag_data);
    diag_clean(&sim.diag_data);

}

