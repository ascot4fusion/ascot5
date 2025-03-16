/**
 * @file dist_6D.c
 * @brief 6D distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../physlib.h"
#include "dist_6D.h"
#include "../gctransform.h"

/**
 * @brief Internal function calculating the index in the histogram array
 */
size_t dist_6D_index(int i_r, int i_phi, int i_z, int i_pr, int i_pphi,
                     int i_pz, int i_time, int i_q, size_t step_7,
                     size_t step_6, size_t step_5, size_t step_4, size_t step_3,
                     size_t step_2, size_t step_1) {
    return (size_t)(i_r)    * step_7
         + (size_t)(i_phi)  * step_6
         + (size_t)(i_z)    * step_5
         + (size_t)(i_pr)   * step_4
         + (size_t)(i_pphi) * step_3
         + (size_t)(i_pz)   * step_2
         + (size_t)(i_time) * step_1
         + (size_t)(i_q);
}

/**
 * @brief Initializes distribution data
 */
int dist_6D_init(dist_6D_data* data) {

    size_t n_q    = (size_t)(data->n_q);
    size_t n_time = (size_t)(data->n_time);
    size_t n_pz   = (size_t)(data->n_pz);
    size_t n_pphi = (size_t)(data->n_pphi);
    size_t n_pr   = (size_t)(data->n_pr);
    size_t n_z    = (size_t)(data->n_z);
    size_t n_phi  = (size_t)(data->n_phi);
    data->step_7 = n_q * n_time * n_pz * n_pphi * n_pr * n_z * n_phi;
    data->step_6 = n_q * n_time * n_pz * n_pphi * n_pr * n_z;
    data->step_5 = n_q * n_time * n_pz * n_pphi * n_pr;
    data->step_4 = n_q * n_time * n_pz * n_pphi;
    data->step_3 = n_q * n_time * n_pz;
    data->step_2 = n_q * n_time;
    data->step_1 = n_q;

    data->histogram = calloc(data->step_7 * (size_t)data->n_r, sizeof(real));
    return data->histogram == NULL;
}

/**
 * @brief Free allocated resources
 */
void dist_6D_free(dist_6D_data* data) {
    free(data->histogram);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void dist_6D_offload(dist_6D_data* data) {
    GPU_MAP_TO_DEVICE(
        data->histogram[0:data->n_r*data->n_phi*data->n_z*data->n_pr*data->n_pphi*data->n_pz*data->n_time*data->n_q]
    )
}

/**
 * @brief Update the histogram from full-orbit particles
 *
 * This function updates the histogram from the particle data. Bins are
 * calculated as vector op and histogram is updates as an atomic operation to
 * avoid race conditions.
 *
 * @param dist pointer to distribution parameter struct
 * @param p_i pointer to SIMD particle struct at the beginning of time step
 * @param p_f pointer to SIMD particle struct at the end of time step
 */
void dist_6D_update_fo(dist_6D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i, int n_running_ref) {

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < n_running_ref; i++) {
        if(p_f->running[i]) {

            int i_r = floor((p_f->r[i] - dist->min_r)
                     / ((dist->max_r - dist->min_r)/dist->n_r));

            real phi = fmod(p_f->phi[i], 2*CONST_PI);
            if(phi < 0) {
                phi += 2*CONST_PI;
            }
            int i_phi = floor((phi - dist->min_phi)
                       / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            int i_z = floor((p_f->z[i] - dist->min_z)
                     / ((dist->max_z - dist->min_z) / dist->n_z));

            int i_pr = floor((p_f->p_r[i] - dist->min_pr)
                      / ((dist->max_pr - dist->min_pr) / dist->n_pr));

            int i_pphi = floor((p_f->p_phi[i] - dist->min_pphi)
                        / ((dist->max_pphi - dist->min_pphi) / dist->n_pphi));

            int i_pz = floor((p_f->p_z[i] - dist->min_pz)
                      / ((dist->max_pz - dist->min_pz) / dist->n_pz));

            int i_time = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            int i_q = floor((p_f->charge[i]/CONST_E - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_r    >= 0 && i_r    <= dist->n_r - 1    &&
               i_phi  >= 0 && i_phi  <= dist->n_phi - 1  &&
               i_z    >= 0 && i_z    <= dist->n_z - 1    &&
               i_pr   >= 0 && i_pr   <= dist->n_pr - 1   &&
               i_pphi >= 0 && i_pphi <= dist->n_pphi - 1 &&
               i_pz   >= 0 && i_pz   <= dist->n_pz - 1   &&
               i_time >= 0 && i_time <= dist->n_time - 1 &&
               i_q    >= 0 && i_q    <= dist->n_q - 1      ) {
                real weight = p_f->weight[i] * (p_f->time[i] - p_i->time[i]);
                size_t index = dist_6D_index(
                    i_r, i_phi, i_z, i_pr, i_pphi, i_pz,
                    i_time, i_q, dist->step_7, dist->step_6, dist->step_5,
                    dist->step_4, dist->step_3, dist->step_2, dist->step_1);
                GPU_ATOMIC
                dist->histogram[index] += weight;
            }
        }
    }
}

/**
 * @brief Update the histogram from guiding-center particles
 *
 * This function updates the histogram from the guiding center data. Bins are
 * calculated as vector op and histogram is updates as an atomic operation to
 * avoid race conditions.
 *
 * @param dist pointer to distribution parameter struct
 * @param p_i pointer to SIMD GC struct at the beginning of time step
 * @param p_f pointer to SIMD GC struct at the end of time step
 */
void dist_6D_update_gc(dist_6D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i) {
    real phi[NSIMD];

    int i_r[NSIMD];
    int i_phi[NSIMD];
    int i_z[NSIMD];
    int i_pr[NSIMD];
    int i_pphi[NSIMD];
    int i_pz[NSIMD];
    int i_time[NSIMD];
    int i_q[NSIMD];

    int ok[NSIMD];
    real weight[NSIMD];

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {

            real pr, pphi, pz;
            real B_dB[12] = {
                p_f->B_r[i], p_f->B_r_dr[i], p_f->B_r_dphi[i], p_f->B_r_dz[i],
                p_f->B_phi[i], p_f->B_phi_dr[i], p_f->B_phi_dphi[i],
                p_f->B_phi_dz[i],
                p_f->B_z[i], p_f->B_z_dr[i], p_f->B_z_dphi[i], p_f->B_z_dz[i]};
            gctransform_pparmuzeta2prpphipz(p_f->mass[i], p_f->charge[i], B_dB,
                                            p_f->phi[i], p_f->ppar[i],
                                            p_f->mu[i], p_f->zeta[i],
                                            &pr, &pphi, &pz);

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

            i_pr[i] = floor((pr - dist->min_pr)
                      / ((dist->max_pr - dist->min_pr) / dist->n_pr));

            i_pphi[i] = floor((pphi - dist->min_pphi)
                        / ((dist->max_pphi - dist->min_pphi) / dist->n_pphi));

            i_pz[i] = floor((pz - dist->min_pz)
                      / ((dist->max_pz - dist->min_pz) / dist->n_pz));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i]/CONST_E - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_r[i]    >= 0 && i_r[i]    <= dist->n_r - 1    &&
               i_phi[i]  >= 0 && i_phi[i]  <= dist->n_phi - 1  &&
               i_z[i]    >= 0 && i_z[i]    <= dist->n_z - 1    &&
               i_pr[i]   >= 0 && i_pr[i]   <= dist->n_pr - 1   &&
               i_pphi[i] >= 0 && i_pphi[i] <= dist->n_pphi - 1 &&
               i_pz[i]   >= 0 && i_pz[i]   <= dist->n_pz - 1   &&
               i_time[i] >= 0 && i_time[i] <= dist->n_time - 1 &&
               i_q[i]    >= 0 && i_q[i]    <= dist->n_q - 1      ) {
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
            size_t index = dist_6D_index(
                i_r[i], i_phi[i], i_z[i], i_pr[i], i_pphi[i], i_pz[i],
                i_time[i], i_q[i], dist->step_7, dist->step_6, dist->step_5,
                dist->step_4, dist->step_3, dist->step_2, dist->step_1);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}
