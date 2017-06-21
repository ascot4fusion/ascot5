/**
 * @file endcond.c
 * @brief Check for end conditions
 */
#include <math.h>
#include "endcond.h"
#include "particle.h"
#include "simulate.h"

void endcond_check_gc(particle_simd_gc* p, sim_data* sim) {
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        /* Max time */
        if(sim->active_endcond & endcond_tmax
            && p->time[i] > sim->tmax) {
            p->endcond[i] |= endcond_tmax;
            p->running[i] = 0;
        }

        /* Wall hit */
        if(sim->active_endcond & endcond_wall) {
            int tile = wall_hit_wall(p->prev_r[i], p->prev_phi[i], p->prev_z[i],
                p->r[i], p->phi[i], p->z[i], &sim->wall_data);
            if(tile > 0) {
                p->walltile[i] = tile;
                p->endcond[i] |= endcond_wall;
                p->running[i] = 0;
            }
        }

        real energy = p->mu[i]*sqrt(p->B_r[i]*p->B_r[i]
            + p->B_phi[i]*p->B_phi[i] + p->B_z[i]*p->B_z[i])
            + 0.5*p->mass[i]*p->vpar[i]*p->vpar[i];

        /* Min energy */
        if(sim->active_endcond & endcond_emin
            && energy < sim->emin) {
            p->endcond[i] |= endcond_emin;
            p->running[i] = 0;
        }
    }
}

void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i, sim_data* sim) {
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        /* Max time */
        if(sim->active_endcond & endcond_tmax
            && p_f->time[i] > sim->tmax) {
            p_f->endcond[i] |= endcond_tmax;
            p_f->running[i] = 0;
        }

        /* Wall hit */
        if(sim->active_endcond & endcond_wall) {
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
        if(sim->active_endcond & endcond_emin
            && energy < sim->emin) {
            p_f->endcond[i] |= endcond_emin;
            p_f->running[i] = 0;
        }
    }
}
