/**
 * @file distributions.c
 * @brief Distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "dist_5D.h"
#include "../particle.h"

/**
 * @brief Internal function calculating the index in the histogram array
 */
#pragma omp declare target
unsigned long dist_5D_index(int i_r, int i_phi, int i_z, int i_vpara,
                            int i_vperp, int i_time, int i_q, int n_phi,
                            int n_z, int n_vpara, int n_vperp, int n_time,
                            int n_q) {
    return i_r    * (n_phi * n_z * n_vpara * n_vperp * n_time * n_q)
        + i_phi   * (n_z * n_vpara * n_vperp * n_time * n_q)
        + i_z     * (n_vpara * n_vperp * n_time * n_q)
        + i_vpara * (n_vperp * n_time * n_q)
        + i_vperp * (n_time * n_q)
        + i_time  * (n_q)
        + i_q;
}
#pragma omp end declare target

/**
 * @brief Frees the offload data
 */
void dist_5D_free_offload(dist_5D_offload_data* offload_data) {
    offload_data->n_r = 0;
    offload_data->min_r = 0;
    offload_data->max_r = 0;
    offload_data->n_phi = 0;
    offload_data->min_phi = 0;
    offload_data->max_phi = 0;
    offload_data->n_z = 0;
    offload_data->min_z = 0;
    offload_data->max_z = 0;
    offload_data->n_vpara = 0;
    offload_data->min_vpara = 0;
    offload_data->max_vpara = 0;
    offload_data->n_vperp = 0;
    offload_data->min_vperp = 0;
    offload_data->max_vperp = 0;
}

/**
 * @brief Initializes distribution from offload data
 */
void dist_5D_init(dist_5D_data* dist_data, dist_5D_offload_data* offload_data,
                  real* offload_array) {
    dist_data->n_r   = offload_data->n_r;
    dist_data->min_r = offload_data->min_r;
    dist_data->max_r = offload_data->max_r;

    dist_data->n_phi   = offload_data->n_phi;
    dist_data->min_phi = offload_data->min_phi;
    dist_data->max_phi = offload_data->max_phi;

    dist_data->n_z   = offload_data->n_z;
    dist_data->min_z = offload_data->min_z;
    dist_data->max_z = offload_data->max_z;

    dist_data->n_vpara   = offload_data->n_vpara;
    dist_data->min_vpara = offload_data->min_vpara;
    dist_data->max_vpara = offload_data->max_vpara;

    dist_data->n_vperp   = offload_data->n_vperp;
    dist_data->min_vperp = offload_data->min_vperp;
    dist_data->max_vperp = offload_data->max_vperp;

    dist_data->n_time   = offload_data->n_time;
    dist_data->min_time = offload_data->min_time;
    dist_data->max_time = offload_data->max_time;

    dist_data->n_q   = offload_data->n_q;
    dist_data->min_q = offload_data->min_q;
    dist_data->max_q = offload_data->max_q;

    dist_data->histogram = &offload_array[0];
}

/**
 * @brief Update the histogram from full-orbit particles
 *
 * This function updates the histogram from the particle data. Bins are
 * calculated as vector op and histogram is updates as an atomic operation to
 * avoid race conditions.
 *
 * @param dist pointer to distribution parameter struct
 * @param p pointer to SIMD particle struct
 */
void dist_5D_update_fo(dist_5D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i) {
    real phi[NSIMD];
    real vpara[NSIMD];
    real vperp[NSIMD];

    int i_r[NSIMD];
    int i_phi[NSIMD];
    int i_z[NSIMD];
    int i_vpara[NSIMD];
    int i_vperp[NSIMD];
    int i_time[NSIMD];
    int i_q[NSIMD];

    int ok[NSIMD];
    real weight[NSIMD];

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            i_r[i] = floor((p_f->r[i] - dist->min_r)
                     / ((dist->max_r - dist->min_r)/dist->n_r));

            phi[i] = fmod(p_f->phi[i], 2*CONST_PI);
            if(phi[i] < 0) {
                phi[i] = phi[i] + 2*CONST_PI;
            }
            i_phi[i] = floor((phi[i] - dist->min_phi)
                       / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            i_z[i] = floor((p_f->z[i] - dist->min_z)
                     / ((dist->max_z - dist->min_z) / dist->n_z));

            vpara[i] = (p_f->rdot[i] * p_f->B_r[i] +
                        (p_f->phidot[i] * p_f->r[i])
                        * p_f->B_phi[i] + p_f->zdot[i] * p_f->B_z[i])
                       / sqrt(p_f->B_r[i]*p_f->B_r[i]
                              +p_f->B_phi[i]*p_f->B_phi[i]
                              + p_f->B_z[i]*p_f->B_z[i]);
            i_vpara[i] = floor((vpara[i] - dist->min_vpara)
                       / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));

            vperp[i] = sqrt(p_f->rdot[i]*p_f->rdot[i] + (p_f->phidot[i]
                                            *p_f->phidot[i]*p_f->r[i]*p_f->r[i])
                            + p_f->zdot[i]*p_f->zdot[i] - vpara[i]*vpara[i]);
            i_vperp[i] = floor((vperp[i] - dist->min_vperp)
                       / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i] - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_r[i]     >= 0  &&  i_r[i]     <= dist->n_r - 1      &&
               i_phi[i]   >= 0  &&  i_phi[i]   <= dist->n_phi - 1    &&
               i_z[i]     >= 0  &&  i_z[i]     <= dist->n_z - 1      &&
               i_vpara[i] >= 0  &&  i_vpara[i] <= dist->n_vpara - 1  &&
               i_vperp[i] >= 0  &&  i_vperp[i] <= dist->n_vperp - 1  &&
               i_time[i]  >= 0  &&  i_time[i]  <= dist->n_time - 1   &&
               i_q[i]     >= 0  &&  i_q[i]     <= dist->n_q - 1        ) {
                ok[i] = 1;
                weight[i] = p_f->weight[i] * (p_f->time[i] - p_i->time[i]);
            }
            else {
                ok[i] = 0;
            }
        }
    }

    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i] && ok[i]) {
            unsigned long index = dist_5D_index(i_r[i], i_phi[i], i_z[i],
                                                i_vpara[i], i_vperp[i],
                                                i_time[i], i_q[i],
                                                dist->n_phi, dist->n_z,
                                                dist->n_vpara, dist->n_vperp,
                                                dist->n_time, dist->n_q);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}

/**
 * @brief Update the histogram from full-orbit particles
 *
 * This function updates the histogram from the particle data. Bins are
 * calculated as vector op and histogram is updates as an atomic operation to
 * avoid race conditions.
 *
 * @param dist pointer to distribution parameter struct
 * @param p pointer to SIMD particle struct
 */
void dist_5D_update_gc(dist_5D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i) {
    real phi[NSIMD];
    real vperp[NSIMD];

    int i_r[NSIMD];
    int i_phi[NSIMD];
    int i_z[NSIMD];
    int i_vpara[NSIMD];
    int i_vperp[NSIMD];
    int i_time[NSIMD];
    int i_q[NSIMD];

    int ok[NSIMD];
    real weight[NSIMD];

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            i_r[i] = floor((p_f->r[i] - dist->min_r)
                     / ((dist->max_r - dist->min_r)/dist->n_r));

            phi[i] = fmod(p_f->phi[i], 2*CONST_PI);
            if(phi[i] < 0) {
                phi[i] = phi[i] + 2*CONST_PI;
            }
            i_phi[i] = floor((phi[i] - dist->min_phi)
                       / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            i_z[i] = floor((p_f->z[i] - dist->min_z)
                    / ((dist->max_z - dist->min_z) / dist->n_z));

            i_vpara[i] = floor((p_f->vpar[i] - dist->min_vpara)
                       / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));

            vperp[i] = sqrt(2 * sqrt(p_f->B_r[i]*p_f->B_r[i]
                                     +p_f->B_phi[i]*p_f->B_phi[i]
                                     +p_f->B_z[i]*p_f->B_z[i])
                            * p_f->mu[i] / p_f->mass[i]);
            i_vperp[i] = floor((vperp[i] - dist->min_vperp)
                       / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i] - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_r[i]     >= 0  &&  i_r[i]     <= dist->n_r - 1      &&
               i_phi[i]   >= 0  &&  i_phi[i]   <= dist->n_phi - 1    &&
               i_z[i]     >= 0  &&  i_z[i]     <= dist->n_z - 1      &&
               i_vpara[i] >= 0  &&  i_vpara[i] <= dist->n_vpara - 1  &&
               i_vperp[i] >= 0  &&  i_vperp[i] <= dist->n_vperp - 1  &&
               i_time[i]  >= 0  &&  i_time[i]  <= dist->n_time - 1   &&
               i_q[i]     >= 0  &&  i_q[i]     <= dist->n_q - 1        ) {
                ok[i] = 1;
                weight[i] = p_f->weight[i] * (p_f->time[i] - p_i->time[i]);
            }
            else {
                ok[i] = 0;
            }
        }
    }

    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i] && ok[i]) {
            unsigned long index = dist_5D_index(i_r[i], i_phi[i], i_z[i],
                                                i_vpara[i], i_vperp[i],
                                                i_time[i], i_q[i],
                                                dist->n_phi,  dist->n_z,
                                                dist->n_vpara, dist->n_vperp,
                                                dist->n_time, dist->n_q);

            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}


void dist_5D_sum(int start, int stop, real* array1, real* array2) {
    for(int i = start; i < stop; i++) {
        array1[i] += array2[i];
    }
}
