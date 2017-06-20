/**
 * @file simulate_fo_lf.c
 * @brief Simulate particles with full orbits using leap frog integrator
 */
#include <stdio.h>
#include <omp.h>
#include <immintrin.h>
#include <math.h>
#include "ascot5.h"
#include "step_fo_lf.h"
#include "step_gc_rk4.h"
#include "wall.h"
#include "distributions.h"
#include "B_field.h"
#include "plasma_1d.h"
#include "interact.h"
#include "simulate.h"
#include "simulate_fo_lf.h"
#include "math.h"
#include "particle.h"

#pragma omp declare target
void print_orbit_fo(particle_simd_fo* p, double t);
#pragma omp end declare target

void simulate_fo_lf(int id, int n_particles, particle* particles,
                    sim_offload_data sim_offload,
                    real* B_offload_array,
                    real* plasma_offload_array,
                    real* wall_offload_array,
                    real* dist_offload_array) {
    sim_data sim;
    sim_init(&sim, &sim_offload);
#ifdef _OPENMP
    double w_start = omp_get_wtime();
#endif
    wall_init(&sim.wall_data, &sim_offload.wall_offload_data,
              wall_offload_array);
#ifdef _OPENMP
    double w_end = omp_get_wtime();
#endif
    #if VERBOSE >= 1
    if(sim.wall_data.type == 3) {
#ifdef _OPENMP
        printf("%d: Initialized wall tree in %.2f s, %.1f MB.\n", id,
               w_end - w_start,
               sim.wall_data.w3d.tree_array_size*sizeof(int) / (1024.0*1024.0));
#endif
    }
    #endif

    B_field_init(&sim.B_data, &sim_offload.B_offload_data, B_offload_array);

    plasma_1d_init(&sim.plasma_data, &sim_offload.plasma_offload_data,
                   plasma_offload_array);

    dist_rzvv_init(&sim.dist_data, &sim_offload.dist_offload_data,
                   dist_offload_array);

    #ifdef REPORTORBIT
    /** @todo fixme */
    /*for(i = 0; i < n; i++) {
        print_orbit_fo(&p[i], 0.0);
    }*/
    #endif 

#ifdef _OPENMP
    double tt_start = omp_get_wtime();
#endif
    int i_next_prt = 0;

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
    #pragma omp parallel
    {
        particle_simd_fo p;
        int i_prt;
        for(int i = 0; i < NSIMD; i++) {
            #pragma omp critical
            i_prt = i_next_prt++;
            if(i_prt < n_particles) {
                particle_to_fo(&particles[i_prt], i_prt, &p, i, &sim.B_data,
                    &sim.E_data);
            }
            else {
                particle_to_fo_dummy(&p, i);
            }
        }

        int collstepdivisor = round(sim.tcollstep / sim.tstep);
        int orbsteps = collstepdivisor;

        int n_running = 0;

        /* Main simulation loop */
        do {
            step_fo_lf(&p, 0, sim.tstep, &sim.B_data);

            #if COULOMBCOLL == 1
            orbsteps++;
            if(orbsteps == collstepdivisor) {
                interact_step_fo_euler(&p, 0, orbsteps*sim.tstep,
                                       &sim.B_data, &sim.plasma_data);
                orbsteps = 0;
            }
            #endif

//            endcond_check_fo(&p, &sim);

            dist_rzvv_update_fo(&sim.dist_data, &p, sim.tstep);

            /* update number of running particles */
            n_running = 0;
            int k;
            #pragma omp simd reduction(+:n_running)
            for(k = 0; k < NSIMD; k++) {
                if(!p.running[k] && p.id[k] >= 0) {
                    fo_to_particle(&p, k, &particles[p.index[k]]);
                    i_prt = i_next_prt++;
                    if(i_prt < n_particles) {
                        particle_to_fo(&particles[i_prt], i_prt, &p, k,
                            &sim.B_data, &sim.E_data);
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
#ifdef _OPENMP
    double tt_end = omp_get_wtime();
#endif
    #if VERBOSE >= 1
#ifdef _OPENMP
    printf("%d: %d particles done in %lf s.\n", id, n_particles,
           tt_end-tt_start);
#endif
    #endif

    /* copy particles back to individual particles */
//    for(i = 0; i < n_particles; i++) {
//    }
}

void print_orbit_fo(particle_simd_fo* p, double t) {
    #ifdef REPORTORBIT
    int k;
    for(k = 0; k < NSIMD; k++) {
        #if REPORTORBIT != -1
        if(p->id[k] == REPORTORBIT)
        #else
        if(p->id[k] != -1)
        #endif
            printf("%d, %le, %le, %le, %le, %le, %le, %le, %le, %le\n",
                (int) p->id[k], t,
                p->r[k], p->phi[k], p->z[k],
                p->rdot[k], p->phidot[k], p->zdot[k],
                p->mass[k], p->charge[k]);
    }
    #endif
}
