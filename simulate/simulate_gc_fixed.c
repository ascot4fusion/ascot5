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

#pragma omp declare target
real simulate_gc_fixed_inidt(sim_data* sim, particle_simd_gc* p, int i);
#pragma omp end declare target

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
 * @param pq particles to be simulated
 * @param sim simulation data
 *
 * @todo See simulate_gc_adaptive.c
 */
void simulate_gc_fixed(particle_queue* pq, sim_data* sim) {
   
    int cycle[NSIMD];
    int err[NSIMD];
    real hin[NSIMD];
    real cputime_last[NSIMD];
    real cputime;

    particle_simd_gc p;  // This array holds current states
    particle_simd_gc p0; // This array stores previous states


    // This is diagnostic specific data which is declared 
    // here to make it thread safe
    diag_storage* diag_strg;
    diag_storage_aquire(&sim->diag_data, &diag_strg);

    for(int i=0; i< NSIMD; i++) {
	p.id[i] = -1;
	p.running[i] = 0;
    }

    /* Initialize running particles */
    int n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

    /* Determine simulation time-step */
    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
	if(cycle[i] > 0) {
	    hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
	    cputime_last[i] = A5_WTIME;
	}
    }

/* MAIN SIMULATION LOOP 
 * - Store current state
 * - Integrate motion due to bacgkround EM-field (orbit-following)
 * - Integrate scattering due to Coulomb collisions
 * - Advance time
 * - Check for end condition(s)
 * - Update diagnostics
 * */
    while(n_running > 0) {
        /* Store marker states */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            p0.r[i]        = p.r[i];
            p0.phi[i]      = p.phi[i];
            p0.z[i]        = p.z[i];
            p0.vpar[i]     = p.vpar[i];
            p0.mu[i]       = p.mu[i];
            p0.theta[i]    = p.theta[i];

            p0.time[i]       = p.time[i];
	    p0.cputime[i]    = p.cputime[i];
	    p0.rho[i]        = p.rho[i];
	    p0.weight[i]     = p.weight[i];
	    p0.cputime[i]    = p.cputime[i]; 
	    p0.rho[i]        = p.rho[i];      
	    p0.pol[i]        = p.pol[i]; 

	    p0.mass[i]       = p.mass[i];
	    p0.charge[i]     = p.charge[i];

            p0.running[i]    = p.running[i];
            p0.endcond[i]    = p.endcond[i];
            p0.walltile[i]   = p.walltile[i];

	    p0.B_r[i]        = p.B_r[i];
	    p0.B_phi[i]      = p.B_phi[i];
	    p0.B_z[i]        = p.B_z[i];

	    p0.B_r_dr[i]     = p.B_r_dr[i];
	    p0.B_r_dphi[i]   = p.B_r_dphi[i];
	    p0.B_r_dz[i]     = p.B_r_dz[i];

	    p0.B_phi_dr[i]   = p.B_phi_dr[i];
	    p0.B_phi_dphi[i] = p.B_phi_dphi[i];
	    p0.B_phi_dz[i]   = p.B_z_dz[i];

	    p0.B_z_dr[i]     = p.B_z_dr[i];
	    p0.B_z_dphi[i]   = p.B_z_dphi[i];
	    p0.B_z_dz[i]     = p.B_z_dz[i];

        }
        
        if(sim->enable_orbfol) {
            step_gc_rk4(&p, hin, &sim->B_data, &sim->E_data);
        }

        if(sim->enable_clmbcol) {
            mccc_step_gc_fixed(&p, &sim->B_data, &sim->plasma_data, hin, err);
        }

	cputime = A5_WTIME;
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
            if(p.running[i]) {
                p.time[i] = p.time[i] + hin[i];
		p.cputime[i] += cputime - cputime_last[i];
		cputime_last[i] = cputime;
            }
        }
        
        endcond_check_gc(&p, &p0, sim);

        diag_update_gc(&sim->diag_data, diag_strg, &p, &p0);

        /* Update running particles */
        n_running = particle_cycle_gc(pq, &p, &sim->B_data, cycle);

        /* Determine simulation time-step */
        #pragma omp simd
        for(int i = 0; i < NSIMD; i++) {
	    if(cycle[i] > 0) {
		hin[i] = simulate_gc_fixed_inidt(sim, &p, i);
		cputime_last[i] = A5_WTIME;
	    }
        }

    }

    diag_storage_discard(diag_strg);
    
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
