/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file simulate_gc_fixed.c
 * @brief Simulate guiding centers using fixed time-step
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
#include "simulate_gc_fixed.h"
#include "step/step_gc_rk4.h"
#include "mccc/mccc.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"

real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);

/**
 * @brief Simulates guiding centers using fixed time-step
 *
 * The simulation includes:
 * - orbit-following with RK4 method
 * - Coulomb collisions with Euler-Maruyama method
 * 
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined either directly or as a fraction of
 * gyrotime.
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
 */
void simulate_gc_fixed(int id, int n_particles, particle* particles,
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
		/* Guiding center to SIMD transformation */
                particle_to_gc(&particles[i_prt], i_prt, &p, i, &sim.B_data);

		/* Determine simulation time-step */
		hin[i] = simulate_gc_fixed_inidt(&sim, &p, i);
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
		
			/* Determine simulation time-step */
			hin[k] = simulate_gc_fixed_inidt(&sim, &p, k);
						
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

/**
 * @brief Calculates time step value
 */
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i) {
    /* Value defined directly by user */
    if(sim->fix_usrdef_use) {
	return sim->fix_usrdef_val;
    }

    /* Value calculated from gyrotime */
    real B = sqrt(p->B_r[i]*p->B_r[i] + p->B_phi[i]*p->B_phi[i] + p->B_z[i]*p->B_z[i]);
    real gamma = 1; // TODO relativistic
    real gyrotime = fabs( CONST_2PI * p->mass[i] * gamma / ( p->charge[i] * B ) );

    return gyrotime/sim->fix_stepsPerGO;
}
