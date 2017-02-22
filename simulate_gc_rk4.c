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
#include "step_gc_rk4.h"
#include "wall.h"
#include "distributions.h"
#include "B_field.h"
#include "E_field.h"
#include "plasma_1d.h"
#include "interact.h"
#include "simulate.h"
#include "math.h"
#include "particle.h"
#include "endcond.h"
#include "simulate_gc_rk4.h"

#pragma omp declare target
void print_orbit_gc(particle_simd_gc* p, double t);
#pragma omp end declare target

void simulate_gc_rk4(int id, int n_particles, particle* particles,
		     sim_offload_data sim_offload,
		     real* B_offload_array,
		     real* E_offload_array,
		     real* plasma_offload_array,
		     real* wall_offload_array,
		     real* dist_offload_array) {
    sim_data sim;
    sim_init(&sim, &sim_offload);

    double w_start=0.0;
#ifdef _OMP
    w_start= omp_get_wtime();
#endif    
    wall_init(&sim.wall_data, &sim_offload.wall_offload_data,
              wall_offload_array);
    double w_end = 0.0;
#ifdef _OMP
    w_end = omp_get_wtime();
#endif 

    #if VERBOSE >= 1
    if(sim.wall_data.type == 3) {
        printf("%d: Initialized wall tree in %.2f s, %.1f MB.\n", id,
               w_end - w_start,
               sim.wall_data.w3d.tree_array_size*sizeof(int) / (1024.0*1024.0));
    }
    #endif

    B_field_init(&sim.B_data, &sim_offload.B_offload_data, B_offload_array);

    E_field_init(&sim.E_data, &sim_offload.E_offload_data, E_offload_array);

    plasma_1d_init(&sim.plasma_data, &sim_offload.plasma_offload_data,
                   plasma_offload_array);

    dist_rzvv_init(&sim.dist_data, &sim_offload.dist_offload_data,
                   dist_offload_array);

    #ifdef REPORTORBIT
    /** @todo fixme */
    /*print_orbit_gc(&p[0], 0.0);*/
    #endif 

    double tt_start = 0.0;
#ifdef _OMP
    tt_start = omp_get_wtime();
#endif
    int i_next_prt = 0;

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
    #pragma omp parallel
    {
        particle_simd_gc p;
        int i, i_prt;
        for(i = 0; i < NSIMD; i++) {
            #pragma omp critical
            i_prt = i_next_prt++;
            if(i_prt < n_particles) {
                particle_to_gc(&particles[i_prt], i_prt, &p, i, &sim.B_data);
            }
            else {
                particle_to_gc_dummy(&p, i);
            }
        }

        int collstepdivisor = round(sim.tcollstep / sim.tstep);
        int orbsteps = collstepdivisor;

        int n_running = 0;

        /* Main simulation loop */
        do {
            step_gc_rk4(&p, 0, sim.tstep, &sim.B_data, &sim.E_data);

            #if COULOMBCOLL == 1
            orbsteps++;
            if(orbsteps >= collstepdivisor) {
                interact_step_gc_euler(&p[i], t, orbsteps*sim.tstep,
                                       &sim.B_data, &sim.plasma_data);
                orbsteps = 0;
            }
            #endif

            endcond_check(&p, &sim);

            dist_rzvv_update_gc(&sim.dist_data, &p, sim.tstep);

            /* update number of running particles */
            n_running = 0;
            int k;
            #pragma omp simd reduction(+:n_running)
            for(k = 0; k < NSIMD; k++) {
                if(!p.running[k] && p.id[k] >= 0) {
                    gc_to_particle(&p, k, &particles[p.index[k]]);
                    #pragma omp critical
                    i_prt = i_next_prt++;
                    if(i_prt < n_particles) {
                        particle_to_gc(&particles[i_prt], i_prt, &p, k,
                            &sim.B_data);
                    }
                    else {
                        p.id[k] = -1;
                    }
                }
                else {
                    n_running += p.running[k];
                }
            }
        } while(n_running > 0);
    }
    double tt_end = 0.0;
#ifdef _OMP
    tt_end = omp_get_wtime();
#endif

    #if VERBOSE >= 1
    printf("%d: %d particles done in %lf s.\n", id, n_particles,
           tt_end-tt_start);
    #endif
}

void print_orbit_gc(particle_simd_gc* p, double t) {
    #ifdef REPORTORBIT
    int k;
    for(k = 0; k < NSIMD; k++) {
        #if REPORTORBIT != -1
        if(p->id[k] == REPORTORBIT)
        #else
        if(p->id[k] != -1)
        #endif
        {
            real absB = math_normc(p->B_r[k], p->B_phi[k], p->B_z[k]);
            real v2 = p->vpar[k]*p->vpar[k] + 2*absB*p->mu[k]/p->mass[k];
            printf("%d, %le, %le, %le, %le, %le, %le, %le, %le, %le, %le\n",
            (int) p->id[k], t,
            p->r[k], p->phi[k], p->z[k],
            p->vpar[k], p->mu[k],
            p->mass[k], p->charge[k],
            0.5*p->mass[k]*v2, p->vpar[k]/sqrt(v2));
        }
    }
    #endif
}
