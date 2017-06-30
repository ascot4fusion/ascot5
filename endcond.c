/**
 * @file endcond.c
 * @brief Check for end conditions
 */
#include <math.h>
#include "endcond.h"
#include "particle.h"
#include "simulate.h"

void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i, sim_data* sim) {
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        /* Max time */
        if(sim->endcond_active & endcond_tmax
            && p_f->time[i] > sim->endcond_maxSimTime) {
            p_f->endcond[i] |= endcond_tmax;
            p_f->running[i] = 0;
        }

        /* Wall hit */
        if(sim->endcond_active & endcond_wall) {
            int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                p_f->r[i], p_f->phi[i], p_f->z[i], &sim->wall_data);
            if(tile > 0) {
                p_f->walltile[i] = tile;
                p_f->endcond[i] |= endcond_wall;
                p_f->running[i] = 0;
            }
        }

        real energy = p_f->mu[i]*sqrt(p_f->B_r[i]*p_f->B_r[i]
            + p_f->B_phi[i]*p_f->B_phi[i] + p_f->B_z[i]*p_f->B_z[i])
            + 0.5*p_f->mass[i]*p_f->vpar[i]*p_f->vpar[i];

        /* Min energy */
        if(sim->endcond_active & endcond_emin
            && energy < sim->endcond_minEkin) {
	    printf("%d %d %d %d\n", sim->endcond_active, endcond_tmax, endcond_wall, endcond_emin);
            p_f->endcond[i] |= endcond_emin;
            p_f->running[i] = 0;
        }
    }
}

void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i, sim_data* sim) {
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        /* Max time */
        if(sim->endcond_active & endcond_tmax
            && p_f->time[i] > sim->endcond_maxSimTime) {
            p_f->endcond[i] |= endcond_tmax;
            p_f->running[i] = 0;
        }

        /* Wall hit */
        if(sim->endcond_active & endcond_wall) {
            int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                p_f->r[i], p_f->phi[i], p_f->z[i], &sim->wall_data);
            if(tile > 0) {
                p_f->walltile[i] = tile;
                p_f->endcond[i] |= endcond_wall;
                p_f->running[i] = 0;
            }
        }

        real energy = 10.0e7; // TODO fix me

        /* Min energy */
        if(sim->endcond_active & endcond_emin
            && energy < sim->endcond_minEkin) {
            p_f->endcond[i] |= endcond_emin;
            p_f->running[i] = 0;
        }
    }
}

void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i, sim_data* sim) {
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        /* Max time */
        if(sim->endcond_active & endcond_tmax
            && p_f->distance[i] > sim->endcond_maxSimTime) {
            p_f->endcond[i] |= endcond_tmax;
            p_f->running[i] = 0;
        }

        /* Wall hit */
        if(sim->endcond_active & endcond_wall) {
            int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                p_f->r[i], p_f->phi[i], p_f->z[i], &sim->wall_data);
            if(tile > 0) {
                p_f->walltile[i] = tile;
                p_f->endcond[i] |= endcond_wall;
                p_f->running[i] = 0;
            }
        }
    }
}
