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
#include "../consts.h"

real simulate_ml_adaptive_inidt(sim_data* sim, particle_simd_ml* p, int i);

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
 * @param particles pointer to marker struct
 *
 * @todo See simulate_gc_adaptive.c
 * @todo Define simulation time by assuming markers travel at the speed of light
 */
void simulate_ml_adaptive(particle_queue* pq, sim_data* sim) {



	

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
    #pragma omp parallel
    {
	/* Arrays needed for the adaptive time step */
	real hin[NSIMD];
	real hout[NSIMD];
	real hnext[NSIMD];
	int cycle[NSIMD];
	real tol = sim->ada_tol_orbfol;
	int i;

        particle_simd_ml p;  // This array holds current states
	particle_simd_ml p0; // This array stores previous states

	for(i=0; i< NSIMD; i++) {
	    p.id[i] = -1;
	    p.running[i] = 0;
	}

	/* Initialize running particles */
	int n_running = particle_cycle_ml(pq, &p, &sim->B_data, cycle);
	
	#pragma omp simd
	for(i = 0; i < NSIMD; i++) {
	    if(cycle[i] > 0) {
		/* Determine initial time-step */
		hin[i] = simulate_ml_adaptive_inidt(sim, &p, i);
	    }
	}


    

/* MAIN SIMULATION LOOP 
 * - Store current state
 * - Integrate motion due to bacgkround EM-field (orbit-following)
 * - Check whether time step was accepted
 *   - NO:  revert to initial state and ignore the end of the loop 
 *          (except CPU_TIME_MAX end condition if this is implemented)
 *   - YES: update particle time, clean redundant Wiener processes, and proceed
 * - Check for end condition(s)
 * - Update diagnostics
 * - 
 * */
        do {
            #pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		/* Store marker states in case time step will be rejected */
	        p0.r[i]        = p.r[i];
		p0.phi[i]      = p.phi[i];
		p0.z[i]        = p.z[i];
		p0.time[i]     = p.time[i];
		p0.running[i]  = p.running[i];
		p0.endcond[i]  = p.endcond[i];
		p0.walltile[i] = p.walltile[i];

		hout[i] = 1.0;
		hnext[i] = 1.0; 
	    }

	    
	    
            if(sim->enable_orbfol) {
	        step_ml_cashkarp(&p, hin, hout, tol, &sim->B_data);
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
		    p.time[i]     = p0.time[i];
		    p.running[i]  = p0.running[i];
		    p.endcond[i]  = p0.endcond[i];
		    p.walltile[i] = p0.walltile[i];
		    hin[i] = -hnext[i];
		    
		}
		else{
		    if(p.running[i]){
			
			p.time[i] = p.time[i] + hin[i]/CONST_C;
			
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

            endcond_check_ml(&p, &p0, sim);
	    
	    diag_update_ml(&sim->diag_data, &p, &p0);

            /* update number of running particles */
	    n_running = particle_cycle_ml(pq, &p, &sim->B_data, cycle);
	
	    #pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
	        if(cycle[i] > 0) {
		    /* Determine initial time-step */
	            hin[i] = simulate_ml_adaptive_inidt(sim, &p, i);
	        }
	    }
	   

        } while(n_running > 0);
        
    }

}

/**
 * @brief Calculates time step value
 */
real simulate_ml_adaptive_inidt(sim_data* sim, particle_simd_ml* p, int i) {

    /* Value calculated from speed of light 
     * assuming initial time step is of the order of 1 cm / c*/
    /* Define this with a compiler parameter */

    return 1.0e-2;

}
