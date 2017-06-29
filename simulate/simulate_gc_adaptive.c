/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file simulate_gc_adaptive.c
 * @brief Simulate guiding centers using adaptive time-step
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
#include "simulate_gc_adaptive.h"
#include "step/step_gc_cashkarp.h"
#include "mccc/mccc.h"
#include "mccc/mccc_wiener.h"
#include "../endcond.h"
#include "../math.h"
#include "../consts.h"

real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i);

/**
 * @brief Simulates guiding centers using adaptive time-step
 *
 * The simulation includes:
 * - orbit-following with Cash-Karp method
 * - Coulomb collisions with Milstein method
 * 
 * The simulation is carried until all marker have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The adaptive time-step is determined by integrator error 
 * tolerances as well as user-defined limits for how much
 * marker state can change during a single time-step.
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
 * @todo what is that id field?
 * @todo time step limits for how much a marker travels in rho or phi
 * @todo integrators are not updating the marker B-field fields as they should 
 */
void simulate_gc_adaptive(int id, int n_particles, particle* particles,
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
	mccc_wienarr *wienarr[NSIMD];
	real hin[NSIMD];
	real hout_orb[NSIMD];
	real hout_col[NSIMD];
	real hnext[NSIMD];
	int err[NSIMD];
	int windex[NSIMD];
	real tol_col = sim.ada_tol_clmbcol;
	real tol_orb = sim.ada_tol_orbfol;


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

		if(sim.enable_clmbcol) {
		    /* Allocate array storing the Wiener processes */
		    wienarr[i] = mccc_wiener_allocate(5,WIENERSLOTS,p.time[i]);
	        }

		/* Determine initial time step */
		hin[i] = simulate_gc_adaptive_inidt(&sim, &p, i);

		
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
 * - Check whether time step was accepted
 *   - NO:  revert to initial state and ignore the end of the loop 
 *          (except CPU_TIME_MAX end condition if this is implemented)
 *   - YES: update particle time, clean redundant Wiener processes, and proceed
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
		p0.vpar[i]     = p.vpar[i];
		p0.mu[i]       = p.mu[i];
		p0.theta[i]    = p.theta[i];
		p0.time[i]     = p.time[i];
		p0.running[i]  = p.running[i];
		p0.endcond[i]  = p.endcond[i];
		p0.walltile[i] = p.walltile[i];
		// Just use some large value here
		hout_orb[i] = 1.0;
		hout_col[i] = 1.0;
		hnext[i] = 1.0; 
	    }

	    
	    if(sim.enable_orbfol) {
	        step_gc_cashkarp(&p, hin, hout_orb, tol_orb,
		                 &sim.B_data, &sim.E_data);
		
		/* Check whether time step was rejected */
                #pragma omp simd
	        for(i = 0; i < NSIMD; i++) {
		    if(p.running[i] && hout_orb[i] < 0){
			p.running[i] = 0;
			hnext[i] = hout_orb[i];
		    }
	        }
	    }

            if(sim.enable_clmbcol) {
	        mccc_step_gc_adaptive(&p, &sim.B_data, &sim.plasma_data,
				      hin, hout_col, wienarr, tol_col, err);
		
		/* Check whether time step was rejected */
                #pragma omp simd
	        for(i = 0; i < NSIMD; i++) {
		    if(p.running[i] && hout_col[i] < 0){
			p.running[i] = 0;
			hnext[i] = hout_col[i];
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
		    p.vpar[i]     = p0.vpar[i];
		    p.mu[i]       = p0.mu[i];
		    p.theta[i]    = p0.theta[i];
		    p.time[i]     = p0.time[i];
		    p.running[i]  = p0.running[i];
		    p.endcond[i]  = p0.endcond[i];
		    p.walltile[i] = p0.walltile[i];
		    hin[i] = -hnext[i];
		    
		}
		else{
		    if(p.running[i]){
			
			p.time[i] = p.time[i] + hin[i];
			
			/* Determine next time step */
			if(hnext[i] > hout_orb[i]) {
			    hnext[i] = hout_orb[i];
			}
			if(hnext[i] > hout_col[i]) {
			    hnext[i] = hout_col[i];
			}
			if(hnext[i] == 1.0) {
			    hnext[i] = hin[i];
			}
			hin[i] = hnext[i];
			if(sim.enable_clmbcol) {
			    /* Clear wiener processes */
			    mccc_wiener_clean(wienarr[i], p.time[i], &err[i]);
		        }
			
		    }
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

		    if(sim.enable_clmbcol) {
		        /* Free the associated Wiener array */
		        mccc_wiener_deallocate(wienarr[k]);
		    }

                    #pragma omp critical
                    i_prt = i_next_prt++;
                    if(i_prt < n_particles) {
                        particle_to_gc(&particles[i_prt], i_prt, &p, k,
				       &sim.B_data);

		        if(sim.enable_clmbcol) {
			    wienarr[k] = mccc_wiener_allocate(5,WIENERSLOTS,p.time[k]);
		        }
		
			/* Determine initial time step */
			hin[k] = simulate_gc_adaptive_inidt(&sim, &p, k);
						
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
real simulate_gc_adaptive_inidt(sim_data* sim, particle_simd_gc* p, int i) {

    /* Just use some large value if no physics are defined */
    real h = 1.0;

    /* Value calculated from gyrotime */
    if(sim->enable_orbfol) {
	real B = sqrt(p->B_r[i]*p->B_r[i] + p->B_phi[i]*p->B_phi[i] + p->B_z[i]*p->B_z[i]);
	real gamma = 1; // TODO relativistic
	real gyrotime = fabs( CONST_2PI * p->mass[i] * gamma / ( p->charge[i] * B ) );
	if(h > gyrotime) {h=gyrotime;}
    }

    /* Value calculated from collision frequency */
    if(sim->enable_clmbcol) {
	real clogab[NSIMD*MAX_SPECIES];
	real Dpara[NSIMD*MAX_SPECIES];
	real DX[NSIMD*MAX_SPECIES];
	real K[NSIMD*MAX_SPECIES];
	real nub[NSIMD*MAX_SPECIES];
	real dQ[NSIMD*MAX_SPECIES];
	real dDparab[NSIMD*MAX_SPECIES];
	mccc_update_gc(p,&sim->B_data,&sim->plasma_data,clogab,Dpara,DX,K,nub,dQ,dDparab);
	real nu = 0;
	for(int k = 0; k < MAX_SPECIES; k++) {
	    nu = nu + nub[MAX_SPECIES*i + k];
	}
	/* Only small angle collisions so divide this by 100 */
	real colltime = 1/(100*nu);
	if(h > colltime) {h=colltime;}
    }

    return h;
}

