/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file simulate_ml_adaptive.c
 * @brief Simulate magnetic field-lines using adaptive time-step
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>
#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"
#include "../wall.h"
#include "../distributions.h"
#include "../diag.h"
#include "../B_field.h"
#include "../E_field.h"
#include "../plasma_1d.h"
#include "simulate_ml_adaptive.h"
#include "step/step_ml_cashkarp.h"
#include "../endcond.h"
#include "../math.h"

/**
 * @brief Simulates magnetic field-lines using adaptive time-step
 *
 * The simulation includes:
 * - orbit-following with Cash-Karp method
 * 
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error 
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
 *
 * Note simulation time is defined by assuming field-lines
 * "travel" at the speed of light.
 *
 * @param id 
 * @param n_particles number of markers to be simulated
 * @param particles pointer to marker struct
 * @param sim_offload simulation offload data struct
 * @param B_offload_array offload array of magnetic field data
 * @param E_offload_array offload array of electric field data
 * @param plasma_offload_array offload array of plasma data
 * @param wall_offload_array offload array of wall data
 * @param diag_offload_array offload array of diagnostics data
 *
 * @todo See simulate_gc_adaptive.c
 * @todo Define simulation time by assuming markers travel at the speed of light
 */
void simulate_ml_adaptive(int id, int n_particles, particle* particles,
			  sim_offload_data sim_offload,
			  real* B_offload_array,
			  real* E_offload_array,
			  real* plasma_offload_array,
			  real* wall_offload_array,
			  real* diag_offload_array) {
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
    diag_init(&sim.diag_data,&sim_offload.diag_offload_data,
	diag_offload_array);
	
   
    int i_next_prt = 0;

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
    #pragma omp parallel
    {
	/* Arrays needed for the adaptive time step */
	real hin[NSIMD];
	real hout[NSIMD];
	real hnext[NSIMD];
	real tol = 1.0e-1;// TODO move to sim


        particle_simd_ml p;  // This array holds current states
	particle_simd_ml p0; // This array stores previous states
        int i, i_prt;

/** MARKER INITIALIZATION */
        for(i = 0; i < NSIMD; i++) {
            #pragma omp critical
            i_prt = i_next_prt++;
            if(i_prt < n_particles) {
		/* Guiding center transformation */
                particle_to_ml(&particles[i_prt], i_prt, &p, i, &sim.B_data);

		// Determine initial time step
		// TODO get this one from physics
		hin[i] = 1.e-4;

		
            }
            else {
		/* Dummy marker to fill NSIMD when we ran out of actual particles */
                particle_to_ml_dummy(&p, i);
            }

	    /* Init dummy particles here, the (required) fields are initialized 
	     * separately at each time step */
	    particle_to_ml_dummy(&p0, i); 
        }


    

/* MAIN SIMULATION LOOP 
 * - Store current state
 * - Integrate motion due to bacgkround EM-field (orbit-following)
 * - Integrate scattering due to Coulomb collisions
 * - Check whether time step was accepted
 *   - NO:  revert to initial state and ignore the end of the loop 
 *          (except CPU_TIME_MAX end condition if this is implemented)
 *   - YES: update particle time, clean redundant Wiener processes, and proceed
 * - Check wall collisions
 * - Check for end condition(s)
 * - Update diagnostics
 * - 
 * */
        int n_running = 0;
        do {
            #pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		/* Store marker states in case time step will be rejected */
	        p0.r[i]        = p.r[i];
		p0.phi[i]      = p.phi[i];
		p0.z[i]        = p.z[i];
		p0.distance[i] = p.distance[i];
		p0.running[i]  = p.running[i];
		p0.endcond[i]  = p.endcond[i];
		p0.walltile[i] = p.walltile[i];

		hout[i] = 1.0;
		hnext[i] = 1.0; 
	    }

	    
	    
            if(sim.enable_orbfol) {
	        step_ml_cashkarp(&p, hin, hout, tol,&sim.B_data);
		/* Check whether time step was rejected */
                #pragma omp simd
	        for(i = 0; i < NSIMD; i++) {
		    if(p.running[i] && hout[i] < 0){
			p.running[i] = 0;
			hnext[i] = hout[i];
			
		    }
	        }
	    }
 
		

            #pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		/* Retrieve marker states in case time step was rejected */
		if(hnext[i] < 0){
		    p.r[i]        = p0.r[i];
		    p.phi[i]      = p0.phi[i];
		    p.z[i]        = p0.z[i];
		    p.distance[i] = p0.distance[i];
		    p.running[i]  = p0.running[i];
		    p.endcond[i]  = p0.endcond[i];
		    p.walltile[i] = p0.walltile[i];
		    hin[i] = -hnext[i];
		    
		}
		else{
		    if(p.running[i]){
			
			p.distance[i] = p.distance[i] + hin[i];
			
			/* Determine next time step */
			if(hnext[i] > hout[i]) {
			    hnext[i] = hout[i];
			}
			if(hnext[i] == 1.0) {
			    hnext[i] = hin[i];
			}
			hin[i] = hnext[i];
			
		    }
		}
	    }

            endcond_check_ml(&p, &p0, &sim);
	    
	    diag_update_ml(&sim.diag_data, &p, &p0);

            /* update number of running particles */
            n_running = 0;
            int k;
            for(k = 0; k < NSIMD; k++) {
                if( !p.running[k] && p.id[k] >= 0) {

                    ml_to_particle(&p, k, &particles[p.index[k]]);

                    #pragma omp critical
                    i_prt = i_next_prt++;
                    if(i_prt < n_particles) {
                        particle_to_ml(&particles[i_prt], i_prt, &p, k,
				       &sim.B_data);
		
			// Determine initial time step
			// TODO get this one from physics
			hin[k] = 1.e-4;
						
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

    diag_clean(&sim.diag_data);

}

