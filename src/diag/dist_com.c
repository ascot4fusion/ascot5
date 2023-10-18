/**
 * @file dist_com.c
 * @brief COM-distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../physlib.h"
#include "../particle.h"
#include "dist_com.h"

/**
 * @brief Internal function calculating the index in the histogram array
 */
unsigned long dist_COM_index(int i_mu, int i_Ekin, int i_Ptor,
                             int n_mu, int n_Ekin, int n_Ptor) {
    return i_mu    * (n_Ekin * n_Ptor)
        + i_Ekin   * (n_Ptor)
        + i_Ptor;
}

/**
 * @brief Frees the offload data
 *
 * @param offload_data pointer to offload data struct
 */
void dist_COM_free_offload(dist_COM_offload_data* offload_data) {
    offload_data->n_mu       = 0;
    offload_data->min_mu     = 0;
    offload_data->max_mu     = 0;
    offload_data->n_Ekin     = 0;
    offload_data->min_Ekin   = 0;
    offload_data->max_Ekin   = 0;
    offload_data->n_Ptor     = 0;
    offload_data->min_Ptor   = 0;
    offload_data->max_Ptor   = 0;
}

/**
 * @brief Initializes distribution from offload data
 *
 * @param dist_data pointer to data struct
 * @param offload_data pointer to offload data struct
 * @param offload_array offload array
 */
void dist_COM_init(dist_COM_data* dist_data,
                   dist_COM_offload_data* offload_data, real* offload_array) {
    dist_data->n_mu       = offload_data->n_mu;
    dist_data->min_mu     = offload_data->min_mu;
    dist_data->max_mu     = offload_data->max_mu;

    dist_data->n_Ekin     = offload_data->n_Ekin;
    dist_data->min_Ekin   = offload_data->min_Ekin;
    dist_data->max_Ekin   = offload_data->max_Ekin;

    dist_data->n_Ptor     = offload_data->n_Ptor;
    dist_data->min_Ptor   = offload_data->min_Ptor;
    dist_data->max_Ptor   = offload_data->max_Ptor;

    size_t n_Ptor = (size_t)(dist_data->n_Ptor);
    size_t n_Ekin = (size_t)(dist_data->n_Ekin);
    dist_data->step_2 = n_Ptor * n_Ekin;
    dist_data->step_1 = n_Ptor;

    dist_data->histogram = &offload_array[0];
}

/**
 * @brief Update the histogram from full-orbit markers
 * @param p_loc pre-allocated pointer to SIMD arrays used in diagnostics kernels
 *        to avoid dynamical allocation
 */
void dist_COM_update_fo(dist_COM_data* dist, B_field_data* Bdata,
                        particle_simd_fo* p_f, particle_simd_fo* p_i, particle_loc* p_loc) {

    real Ekin;
    real Ptor;
    real Bnorm;
    real psi;
    real mu;
    real xi;
    real pnorm;
    real ppar;

    real* weight  = p_loc->r_arr1;
    int*  i_mu    = p_loc->i_arr1;
    int*  i_Ekin  = p_loc->i_arr2;
    int*  i_Ptor  = p_loc->i_arr3;
    int*  ok      = p_loc->i_arr4;

//    #pragma omp simd
#pragma acc data present(weight[0:NSIMD],i_mu[0:NSIMD],i_Ekin[0:NSIMD],i_Ptor[0:NSIMD],ok[0:NSIMD] )
    {
    GPU_PARALLEL_LOOP_ALL_LEVELS  
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            real Ekin, Ptor, Bnorm, psi, mu, xi, pnorm, ppar;

            B_field_eval_psi(&psi, p_f->r[i], p_f->phi[i], p_f->z[i],
                             p_f->time[i], Bdata);

            real B [] = {p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]};
            Bnorm = math_normc(p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]);

            real p [] = {p_f->p_r[i], p_f->p_phi[i], p_f->p_z[i]};
            pnorm = math_normc(p_f->p_r[i], p_f->p_phi[i], p_f->p_z[i]);
            ppar = math_dot(p, B) / Bnorm;
            xi = ppar / pnorm;

            mu = physlib_gc_mu(p_f->mass[i], pnorm, xi, Bnorm);
            i_mu[i] = floor((mu - dist->min_mu)
                            / ((dist->max_mu - dist->min_mu)/dist->n_mu));
            Ekin = physlib_Ekin_pnorm(p_f->mass[i], pnorm);

            i_Ekin[i] = floor((Ekin - dist->min_Ekin)
                         / ((dist->max_Ekin - dist->min_Ekin) / dist->n_Ekin));
            Ptor = phys_ptoroid_fo(
                p_f->charge[i], p_f->r[i], p_f->p_phi[i], psi);

            i_Ptor[i] = floor((Ptor - dist->min_Ptor)
                         / ((dist->max_Ptor - dist->min_Ptor)/dist->n_Ptor));

            if(i_mu[i]   >= 0 && i_mu[i]   <= dist->n_mu - 1   &&
               i_Ekin[i] >= 0 && i_Ekin[i] <= dist->n_Ekin - 1 &&
               i_Ptor[i] >= 0 && i_Ptor[i] <= dist->n_Ptor - 1 ) {
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
            unsigned long index = dist_COM_index(i_mu[i], i_Ekin[i], i_Ptor[i],
                                                dist->n_mu,  dist->n_Ekin,
                                                dist->n_Ptor);
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
 * @param Bdata pointer to magnetic field data
 * @param p_f pointer to SIMD gc struct at the end of current time step
 * @param p_i pointer to SIMD gc struct at the start of current time step
 */
void dist_COM_update_gc(dist_COM_data* dist, B_field_data* Bdata,
                        particle_simd_gc* p_f, particle_simd_gc* p_i) {
    int i_mu[NSIMD];
    int i_Ekin[NSIMD];
    int i_Ptor[NSIMD];

    int ok[NSIMD];
    real weight[NSIMD];

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            real Ekin, Ptor, B, psi;

            B_field_eval_psi(&psi, p_f->r[i], p_f->phi[i], p_f->z[i],
                             p_f->time[i], Bdata);

            B = math_normc(p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]);
            i_mu[i] = floor((p_f->mu[i] - dist->min_mu)
                            / ((dist->max_mu - dist->min_mu)/dist->n_mu));

            Ekin = physlib_Ekin_ppar(p_f->mass[i], p_f->mu[i], p_f->ppar[i], B);

            i_Ekin[i] = floor((Ekin - dist->min_Ekin)
                         / ((dist->max_Ekin - dist->min_Ekin) / dist->n_Ekin));
            Ptor = phys_ptoroid_gc(p_f->charge[i], p_f->r[i], p_f->ppar[i], psi,
                                   B, p_f->B_phi[i]);

            i_Ptor[i] = floor((Ptor - dist->min_Ptor)
                         / ((dist->max_Ptor - dist->min_Ptor)/dist->n_Ptor));

            if(i_mu[i]   >= 0 && i_mu[i]   <= dist->n_mu - 1   &&
               i_Ekin[i] >= 0 && i_Ekin[i] <= dist->n_Ekin - 1 &&
               i_Ptor[i] >= 0 && i_Ptor[i] <= dist->n_Ptor - 1 ) {
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
            size_t index = dist_COM_index(i_mu[i], i_Ekin[i], i_Ptor[i],
                                          dist->step_2, dist->step_1);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}
