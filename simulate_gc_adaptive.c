/**
 * @file simulate_gc_rk4.c
 * @brief Simulate particles with guiding center using RK4 integrator
 */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>
#include "ascot5.h"
#include "step_gc_cashkarp.h"
#include "wall.h"
#include "distributions.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma_1d.h"
#include "simulate.h"
#include "math.h"
#include "particle.h"
#include "endcond.h"
#include "simulate_gc_adaptive.h"
#include "orbit_write.h"
#include "mccc/mccc.h"
#include "mccc/mccc_wiener.h"

void simulate_gc_adaptive(int id, int n_particles, particle* particles,
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
    dist_rzvv_init(&sim.dist_data, &sim_offload.dist_offload_data,
                   dist_offload_array);
    writeorbit_guidingcenter* writedata = malloc( sizeof(writeorbit_guidingcenter));
    int writeMarker[NSIMD];
    int writeN;
    int slot = 0;

    /* Arrays needed for the adaptive time step */
    mccc_wienarr *wienarr[NSIMD];
    real hin[NSIMD];
    real hout[NSIMD];
    real hnext[NSIMD];
    int err[NSIMD];
    int windex[NSIMD];
    real tol_col = 1.0;
    real tol_orb = 1.0e-8;
	
   
    int i_next_prt = 0;

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
#pragma omp parallel
    {
        particle_simd_gc p;
	particle_simd_gc p0; // This array stores previous states
        int i, i_prt;
/** MARKER INITIALIZATION */
        for(i = 0; i < NSIMD; i++) {
#pragma omp critical
            i_prt = i_next_prt++;
            if(i_prt < n_particles) {
		/* Guiding center transformation */
                particle_to_gc(&particles[i_prt], i_prt, &p, i, &sim.B_data);

		/* Allocate wienarr. NDIM = 5 because we are simulating 
		    guiding centers */
		wienarr[i] = mccc_wiener_allocate(5,WIENERSLOTS,p.time[i]);
		

		// Determine initial time step
		// TODO get this one from physics
		hin[i] = 1.e-6;
            }
            else {
		/* Dummy marker to fill NSIMD when we ran out of actual particles */
                particle_to_gc_dummy(&p, i);
            }

	    /* Init dummy particles here, the (required) fields are initialized 
	     * separately at each time step */
	    particle_to_gc_dummy(&p0, i); 
        }
	
        int n_running = 0;
/* MAIN SIMULATION LOOP */
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
		hnext[i] = 0.0;
	    }

	    
	    
//#if ORBITFOLLOWING == 1
	    step_gc_cashkarp(&p, p.time, hin, hout, tol_orb,
	    		     &sim.B_data, &sim.E_data);
#pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		if(hnext[i] == 0 || hnext[i] > hout[i]){
		    hnext[i] = hout[i];
		}
	    }
//#endif

//#if COULOMBCOLL == 1

	    mccc_step_gc_adaptive(&p, &sim.B_data, &sim.plasma_data,
				  hin, hout, wienarr, tol_col, err);
#pragma omp simd
	    for(i = 0; i < NSIMD; i++) {
		if(hnext[i] == 0 || hnext[i] > hout[i]){
		    hnext[i] = hout[i];
		}
	    }
	    
                
//#endif
 
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
			// Clear wiener processes
			p.time[i] = p.time[i] + hin[i];
			
			mccc_wiener_clean(wienarr[i], p.time[i], &err[i]);
			
			hin[i] = hnext[i];
		    }
		}
	    }
            endcond_check(&p, &sim);
	    
            //Not yet compatible with time step rejections
	    //dist_rzvv_update_gc(&sim.dist_data, &p, sim.tstep);

            /* update number of running particles */
            n_running = 0;
            int k;
            for(k = 0; k < NSIMD; k++) {
		writeMarker[k] = -1;
                if( p.running[k] == 0 && p.id[k] >= 0) {
		    mccc_wiener_deallocate(wienarr[k]);
                    gc_to_particle(&p, k, &particles[p.index[k]]);
#pragma omp critical
                    i_prt = i_next_prt++;
                    if(i_prt < n_particles) {
                        particle_to_gc(&particles[i_prt], i_prt, &p, k,
				       &sim.B_data);
			/* Allocate wienarr. NDIM = 5 because we are simulating 
			   guiding centers */
			wienarr[k] = mccc_wiener_allocate(5,WIENERSLOTS,p.time[k]);
		
			// Determine initial time step
			// TODO get this one from physics
			hin[i] = 1.e-8;
						
                    }
                    else {
                        p.id[k] = -1;
                    }
                }
                else {
#pragma omp critical
		    writeMarker[k] = n_running;
#pragma omp critical
                    n_running += p.running[k];		
                }
            }
	    //printf("asda");
	    writeorbit_store_gc2guidingcenter(p, writedata, writeMarker, n_running, slot);
	    slot = slot + n_running;
	   

        } while(n_running > 0);
	
	FILE* fn = fopen("gcorbits.test","w");
	while(slot > 0){
	    fprintf(fn,"%d, %le, %le, %le, %le, %le, %le\n",
	    writedata->id[slot], writedata->time[slot], writedata->r[slot], writedata->phi[slot], writedata->z[slot],
	    writedata->vpar[slot], writedata->mu[slot]);

	    slot = slot -1;
		
	}
	fclose(fn);
	free(writedata);
    }

}

