/**
 * @file dist_5D.c
 * @brief 5D distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../physlib.h"
#include "dist_5D.h"
#include "../particle.h"

/**
 * @brief Function for calculating the index in the histogram array
 */
unsigned long dist_5D_index(int i_r, int i_phi, int i_z, int i_ppara,
                            int i_pperp, int i_time, int i_q, int n_phi,
                            int n_z, int n_ppara, int n_pperp, int n_time,
                            int n_q) {
    return i_r    * (n_phi * n_z * n_ppara * n_pperp * n_time * n_q)
        + i_phi   * (n_z * n_ppara * n_pperp * n_time * n_q)
        + i_z     * (n_ppara * n_pperp * n_time * n_q)
        + i_ppara * (n_pperp * n_time * n_q)
        + i_pperp * (n_time * n_q)
        + i_time  * (n_q)
        + i_q;
}

/**
 * @brief Frees the offload data
 *
 * @param offload_data pointer to offload data struct
 */
void dist_5D_free_offload(dist_5D_offload_data* offload_data) {
    offload_data->n_r       = 0;
    offload_data->min_r     = 0;
    offload_data->max_r     = 0;
    offload_data->n_phi     = 0;
    offload_data->min_phi   = 0;
    offload_data->max_phi   = 0;
    offload_data->n_z       = 0;
    offload_data->min_z     = 0;
    offload_data->max_z     = 0;
    offload_data->n_ppara   = 0;
    offload_data->min_ppara = 0;
    offload_data->max_ppara = 0;
    offload_data->n_pperp   = 0;
    offload_data->min_pperp = 0;
    offload_data->max_pperp = 0;
}

/**
 * @brief Initializes distribution from offload data
 *
 * @param dist_data pointer to data struct
 * @param offload_data pointer to offload data struct
 * @param offload_array offload array
 */
void dist_5D_init(dist_5D_data* dist_data, dist_5D_offload_data* offload_data,
                  real* offload_array) {
    dist_data->n_r       = offload_data->n_r;
    dist_data->min_r     = offload_data->min_r;
    dist_data->max_r     = offload_data->max_r;

    dist_data->n_phi     = offload_data->n_phi;
    dist_data->min_phi   = offload_data->min_phi;
    dist_data->max_phi   = offload_data->max_phi;

    dist_data->n_z       = offload_data->n_z;
    dist_data->min_z     = offload_data->min_z;
    dist_data->max_z     = offload_data->max_z;

    dist_data->n_ppara   = offload_data->n_ppara;
    dist_data->min_ppara = offload_data->min_ppara;
    dist_data->max_ppara = offload_data->max_ppara;

    dist_data->n_pperp   = offload_data->n_pperp;
    dist_data->min_pperp = offload_data->min_pperp;
    dist_data->max_pperp = offload_data->max_pperp;

    dist_data->n_time    = offload_data->n_time;
    dist_data->min_time  = offload_data->min_time;
    dist_data->max_time  = offload_data->max_time;

    dist_data->n_q       = offload_data->n_q;
    dist_data->min_q     = offload_data->min_q;
    dist_data->max_q     = offload_data->max_q;

    size_t n_q     = (size_t)(dist_data->n_q);
    size_t n_time  = (size_t)(dist_data->n_time);
    size_t n_pperp = (size_t)(dist_data->n_pperp);
    size_t n_ppara = (size_t)(dist_data->n_ppara);
    size_t n_z     = (size_t)(dist_data->n_z);
    size_t n_phi   = (size_t)(dist_data->n_phi);
    dist_data->step_6 = n_q * n_time * n_pperp * n_ppara * n_z * n_phi;
    dist_data->step_5 = n_q * n_time * n_pperp * n_ppara * n_z;
    dist_data->step_4 = n_q * n_time * n_pperp * n_ppara;
    dist_data->step_3 = n_q * n_time * n_pperp;
    dist_data->step_2 = n_q * n_time;
    dist_data->step_1 = n_q;

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
 * @param p_f pointer to SIMD particle struct at the end of current time step
 * @param p_i pointer to SIMD particle struct at the start of current time step
 * @param p_loc pre-allocated pointer to SIMD arrays used in diagnostics kernels
 *        to avoid dynamical allocation
 */
void dist_5D_update_fo(dist_5D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i, particle_loc* p_loc) {

    real* phi = p_loc->r_arr1;
    real* ppara = p_loc->r_arr2;
    real* pperp = p_loc->r_arr3;

    int* i_r = p_loc->i_arr1;
    int* i_phi = p_loc->i_arr2;
    int* i_z = p_loc->i_arr3;
    int* i_ppara = p_loc->i_arr4;
    int* i_pperp = p_loc->i_arr5;
    int* i_time = p_loc->i_arr6;
    int* i_q = p_loc->i_arr7;

    int* ok = p_loc->i_arr8;
    real* weight = p_loc->r_arr4;

#pragma acc data present(phi[0:NSIMD],ppara[0:NSIMD],pperp[0:NSIMD],i_r[0:NSIMD],i_phi[0:NSIMD],i_z[0:NSIMD],i_ppara[0:NSIMD],i_pperp[0:NSIMD],i_time[0:NSIMD],i_q[0:NSIMD],ok[0:NSIMD],weight[0:NSIMD])
    {
    GPU_PARALLEL_LOOP_ALL_LEVELS
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

            ppara[i] = (  p_f->p_r[i]   * p_f->B_r[i]
                        + p_f->p_phi[i] * p_f->B_phi[i]
                        + p_f->p_z[i]   * p_f->B_z[i])
                       / sqrt(  p_f->B_r[i]  * p_f->B_r[i]
                              + p_f->B_phi[i]* p_f->B_phi[i]
                              + p_f->B_z[i]  * p_f->B_z[i]);
            i_ppara[i] = floor((ppara[i] - dist->min_ppara)
                       / ((dist->max_ppara - dist->min_ppara) / dist->n_ppara));

            pperp[i] = sqrt(
                    p_f->p_r[i]   * p_f->p_r[i]
                  + p_f->p_phi[i] * p_f->p_phi[i]
                  + p_f->p_z[i]   * p_f->p_z[i]
                  - ppara[i] * ppara[i]);
            i_pperp[i] = floor((pperp[i] - dist->min_pperp)
                       / ((dist->max_pperp - dist->min_pperp) / dist->n_pperp));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i]/CONST_E - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_r[i]     >= 0  &&  i_r[i]     <= dist->n_r - 1      &&
               i_phi[i]   >= 0  &&  i_phi[i]   <= dist->n_phi - 1    &&
               i_z[i]     >= 0  &&  i_z[i]     <= dist->n_z - 1      &&
               i_ppara[i] >= 0  &&  i_ppara[i] <= dist->n_ppara - 1  &&
               i_pperp[i] >= 0  &&  i_pperp[i] <= dist->n_pperp - 1  &&
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

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i] && ok[i]) {
            unsigned long index = dist_5D_index(i_r[i], i_phi[i], i_z[i],
                                                i_ppara[i], i_pperp[i],
                                                i_time[i], i_q[i],
                                                dist->n_phi, dist->n_z,
                                                dist->n_ppara, dist->n_pperp,
                                                dist->n_time, dist->n_q);
	    GPU_ATOMIC
            dist->histogram[index] += weight[i];
        }
    }
}
}
/**
 * @brief Update the histogram from guiding center markers
 *
 * This function updates the histogram from the marker data. Bins are
 * calculated as vector op and histogram is updates as an atomic operation to
 * avoid race conditions.
 *
 * @param dist pointer to distribution parameter struct
 * @param p_f pointer to SIMD gc struct at the end of current time step
 * @param p_i pointer to SIMD gc struct at the start of current time step
 */
void dist_5D_update_gc(dist_5D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i) {
    real phi[NSIMD];
    real pperp[NSIMD];

    int i_r[NSIMD];
    int i_phi[NSIMD];
    int i_z[NSIMD];
    int i_ppara[NSIMD];
    int i_pperp[NSIMD];
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

            i_ppara[i] = floor((p_f->ppar[i] - dist->min_ppara)
                       / ((dist->max_ppara - dist->min_ppara) / dist->n_ppara));

            pperp[i] = sqrt(2 * sqrt(  p_f->B_r[i]   * p_f->B_r[i]
                                     + p_f->B_phi[i] * p_f->B_phi[i]
                                     + p_f->B_z[i]   * p_f->B_z[i] )
                            * p_f->mu[i] * p_f->mass[i]);
            i_pperp[i] = floor((pperp[i] - dist->min_pperp)
                       / ((dist->max_pperp - dist->min_pperp) / dist->n_pperp));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i]/CONST_E - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_r[i]     >= 0  &&  i_r[i]     <= dist->n_r - 1      &&
               i_phi[i]   >= 0  &&  i_phi[i]   <= dist->n_phi - 1    &&
               i_z[i]     >= 0  &&  i_z[i]     <= dist->n_z - 1      &&
               i_ppara[i] >= 0  &&  i_ppara[i] <= dist->n_ppara - 1  &&
               i_pperp[i] >= 0  &&  i_pperp[i] <= dist->n_pperp - 1  &&
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
            size_t index = dist_5D_index(
                i_r[i], i_phi[i], i_z[i], i_ppara[i], i_pperp[i], i_time[i],
                i_q[i], dist->step_6, dist->step_5, dist->step_4,
                dist->step_3, dist->step_2, dist->step_1);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}
