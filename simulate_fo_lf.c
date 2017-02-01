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

    double w_start = omp_get_wtime();
    wall_init(&sim.wall_data, &sim_offload.wall_offload_data,
              wall_offload_array);
    double w_end = omp_get_wtime();
    #if VERBOSE >= 1
    if(sim.wall_data.type == 3) {
        printf("%d: Initialized wall tree in %.2f s, %.1f MB.\n", id,
               w_end - w_start,
               sim.wall_data.w3d.tree_array_size*sizeof(int) / (1024.0*1024.0));
    }
    #endif

    B_field_init(&sim.B_data, &sim_offload.B_offload_data, B_offload_array);

    plasma_1d_init(&sim.plasma_data, &sim_offload.plasma_offload_data,
                   plasma_offload_array);

    dist_rzvv_init(&sim.dist_data, &sim_offload.dist_offload_data,
                   dist_offload_array);

    int n = n_particles / NSIMD;
    if(n_particles % NSIMD != 0) {
        n++;
    }
    particle_simd_fo* p = (particle_simd_fo*) _mm_malloc(
                                                n*sizeof(particle_simd_fo), 64);
    int i;
    for(i = 0; i < n; i++) {
        int j;
        for(j = 0; j < NSIMD; j++) {
            if(i*NSIMD + j < n_particles) {
                p[i].r[j] = particles[i*NSIMD+j].r;
                p[i].phi[j] = particles[i*NSIMD+j].phi;
                p[i].z[j] = particles[i*NSIMD+j].z;
                p[i].rdot[j] = particles[i*NSIMD+j].rdot;
                p[i].phidot[j] = particles[i*NSIMD+j].phidot;
                p[i].zdot[j] = particles[i*NSIMD+j].zdot;
                p[i].mass[j] = particles[i*NSIMD+j].mass;
                p[i].charge[j] = particles[i*NSIMD+j].charge;
                p[i].weight[j] = particles[i*NSIMD+j].weight;
                p[i].id[j] = particles[i*NSIMD+j].id;
                p[i].running[j] = 1;
            }
            else {
                /* if there particles don't fill the last struct, remaining
                   particles are dummies that won't be simulated */
                p[i].r[j] = 1;
                p[i].phi[j] = 1;
                p[i].z[j] = 1;
                p[i].rdot[j] = 1;
                p[i].phidot[j] = 1;
                p[i].zdot[j] = 1;
                p[i].mass[j] = 1;
                p[i].charge[j] = 1;
                p[i].weight[j] = 0;
                p[i].id[j] = -1;
                p[i].running[j] = 0;
            }
        }
    }

    #ifdef REPORTORBIT
    /** @todo fixme */
    /*for(i = 0; i < n; i++) {
        print_orbit_fo(&p[i], 0.0);
    }*/
    #endif 

    double tt_start = omp_get_wtime();

    /* SIMD particle structs will be computed in parallel with the maximum
     * number of threads available on the platform */
    #pragma omp parallel for 
    for(i = 0; i < n; i++) {

        double t;
        double trprev = sim.t0;

        int orbsteps = 0;
        int collstepdivisor = round(sim.tcollstep / sim.tstep);

        /* Main simulation loop */
        for(t = sim.t0; t < sim.tmax; t += sim.tstep) {
            int k;
            #pragma omp simd
            for(k = 0; k < NSIMD; k++) {
                p[i].prev_r[k] = p[i].r[k];
                p[i].prev_phi[k] = p[i].phi[k];
                p[i].prev_z[k] = p[i].z[k];
            }

            /* Perform a single rk4 step. Vectorization magic happens inside
             * the function. */
            step_fo_lf(&p[i], t, sim.tstep, &sim.B_data);

            /* Check for wall collisions. This is performed as vector op */
            #pragma omp simd
            for(k = 0; k < NSIMD; k++) {
                if(wall_hit_wall(p[i].prev_r[k], p[i].prev_phi[k],
                                 p[i].prev_z[k], p[i].r[k], p[i].phi[k],
                                 p[i].z[k], &sim.wall_data)) 
                {
                    p[i].running[k] = 0;
                }
            }

            /* Perform collision step */
            #if COULOMBCOLL == 1
            orbsteps++;
            if(orbsteps == collstepdivisor) {
                interact_step_fo_euler(&p[i], t, orbsteps*sim.tstep,
                                       &sim.B_data, &sim.plasma_data);
                orbsteps = 0;
            }
            #endif

            /* Print particle orbit */
            #ifdef REPORTORBIT
            if((t-trprev) > sim.trstep) {
                print_orbit_fo(&p[i], t+sim.tstep);
                trprev = t;
            }
            #endif 

            /* Update histogram, vectorization inside */
            dist_rzvv_update_fo(&sim.dist_data, &p[i], sim.tstep);
        }
    }
    double tt_end = omp_get_wtime();

    #if VERBOSE >= 1
    printf("%d: %d particles done in %lf s.\n", id, n_particles,
           tt_end-tt_start);
    #endif

    /* copy particles back to individual particles */
    for(i = 0; i < n_particles; i++) {
        particles[i].r = p[i/NSIMD].r[i % NSIMD];
        particles[i].phi = p[i/NSIMD].phi[i % NSIMD];
        particles[i].z = p[i/NSIMD].z[i % NSIMD];
        particles[i].rdot = p[i/NSIMD].rdot[i % NSIMD];
        particles[i].phidot = p[i/NSIMD].phidot[i % NSIMD];
        particles[i].zdot = p[i/NSIMD].zdot[i % NSIMD];
        particles[i].mass = p[i/NSIMD].mass[i % NSIMD];
        particles[i].charge = p[i/NSIMD].charge[i % NSIMD];
        particles[i].weight = p[i/NSIMD].weight[i % NSIMD];
        particles[i].id = p[i/NSIMD].id[i % NSIMD];
        particles[i].running = p[i/NSIMD].running[i % NSIMD];
    }
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
