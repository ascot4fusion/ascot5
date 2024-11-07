/**
 * @file dist_rho6D.c
 * @brief rho 6D distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../physlib.h"
#include "../gctransform.h"
#include "dist_rho6D.h"

/**
 * @brief Internal function calculating the index in the histogram array
 */
size_t dist_rho6D_index(int i_rho, int i_theta, int i_phi, int i_pr, int i_pphi,
                        int i_pz, int i_time, int i_q, size_t step_7,
                        size_t step_6, size_t step_5, size_t step_4,
                        size_t step_3, size_t step_2, size_t step_1) {
    return (size_t)(i_rho)   * step_7
         + (size_t)(i_theta) * step_6
         + (size_t)(i_phi)   * step_5
         + (size_t)(i_pr)    * step_4
         + (size_t)(i_pphi)  * step_3
         + (size_t)(i_pz)    * step_2
         + (size_t)(i_time)  * step_1
         + (size_t)(i_q);
}

/**
 * @brief Frees the offload data
 */
void dist_rho6D_free_offload(dist_rho6D_offload_data* offload_data) {
    offload_data->n_rho     = 0;
    offload_data->min_rho   = 0;
    offload_data->max_rho   = 0;
    offload_data->n_theta   = 0;
    offload_data->min_theta = 0;
    offload_data->max_theta = 0;
    offload_data->n_phi     = 0;
    offload_data->min_phi   = 0;
    offload_data->max_phi   = 0;
    offload_data->n_pr      = 0;
    offload_data->min_pr    = 0;
    offload_data->max_pr    = 0;
    offload_data->n_pphi    = 0;
    offload_data->min_pphi  = 0;
    offload_data->max_pphi  = 0;
    offload_data->n_pz      = 0;
    offload_data->min_pz    = 0;
    offload_data->max_pz    = 0;
}

/**
 * @brief Initializes distribution from offload data
 */
void dist_rho6D_init(dist_rho6D_data* dist_data,
                     dist_rho6D_offload_data* offload_data,
                     real* offload_array) {
    dist_data->n_rho     = offload_data->n_rho;
    dist_data->min_rho   = offload_data->min_rho;
    dist_data->max_rho   = offload_data->max_rho;

    dist_data->n_theta   = offload_data->n_theta;
    dist_data->min_theta = offload_data->min_theta;
    dist_data->max_theta = offload_data->max_theta;

    dist_data->n_phi     = offload_data->n_phi;
    dist_data->min_phi   = offload_data->min_phi;
    dist_data->max_phi   = offload_data->max_phi;

    dist_data->n_pr      = offload_data->n_pr;
    dist_data->min_pr    = offload_data->min_pr;
    dist_data->max_pr    = offload_data->max_pr;

    dist_data->n_pphi    = offload_data->n_pphi;
    dist_data->min_pphi  = offload_data->min_pphi;
    dist_data->max_pphi  = offload_data->max_pphi;

    dist_data->n_pz      = offload_data->n_pz;
    dist_data->min_pz    = offload_data->min_pz;
    dist_data->max_pz    = offload_data->max_pz;

    dist_data->n_time    = offload_data->n_time;
    dist_data->min_time  = offload_data->min_time;
    dist_data->max_time  = offload_data->max_time;

    dist_data->n_q       = offload_data->n_q;
    dist_data->min_q     = offload_data->min_q;
    dist_data->max_q     = offload_data->max_q;

    size_t n_q     = (size_t)(dist_data->n_q);
    size_t n_time  = (size_t)(dist_data->n_time);
    size_t n_pz    = (size_t)(dist_data->n_pz);
    size_t n_pphi  = (size_t)(dist_data->n_pphi);
    size_t n_pr    = (size_t)(dist_data->n_pr);
    size_t n_phi   = (size_t)(dist_data->n_phi);
    size_t n_theta = (size_t)(dist_data->n_theta);
    dist_data->step_7 = n_q * n_time * n_pz * n_pphi * n_pr * n_phi * n_theta;
    dist_data->step_6 = n_q * n_time * n_pz * n_pphi * n_pr * n_phi;
    dist_data->step_5 = n_q * n_time * n_pz * n_pphi * n_pr;
    dist_data->step_4 = n_q * n_time * n_pz * n_pphi;
    dist_data->step_3 = n_q * n_time * n_pz;
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
 * @param p_i pointer to SIMD particle struct at the beginning of time step
 * @param p_f pointer to SIMD particle struct at the end of time step
 */
void dist_rho6D_update_fo(dist_rho6D_data* dist, particle_simd_fo* p_f,
                          particle_simd_fo* p_i) {

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < p_f->n_mrk_ref; i++) {
        if(p_f->running[i]) {

            int i_rho = floor((p_f->rho[i] - dist->min_rho)
                             / ((dist->max_rho - dist->min_rho)/dist->n_rho));

            real phi = fmod(p_f->phi[i], 2*CONST_PI);
            if(phi < 0) {
                phi = phi + 2*CONST_PI;
            }
            int i_phi = floor((phi - dist->min_phi)
                             / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            real theta = fmod(p_f->theta[i], 2*CONST_PI);
            if(theta < 0) {
                theta += 2*CONST_PI;
            }
            int i_theta = floor((theta - dist->min_theta)
                             / ((dist->max_theta - dist->min_theta)
                                / dist->n_theta));

            int i_pr = floor((p_f->p_r[i] - dist->min_pr)
                            / ((dist->max_pr - dist->min_pr) / dist->n_pr));

            int i_pphi = floor((p_f->p_phi[i] - dist->min_pphi)
                              / ((dist->max_pphi - dist->min_pphi)
                                 / dist->n_pphi));

            int i_pz = floor((p_f->p_z[i] - dist->min_pz)
                            / ((dist->max_pz - dist->min_pz) / dist->n_pz));

            int i_time = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            int i_q = floor((p_f->charge[i]/CONST_E - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_rho  >= 0 && i_rho  <= dist->n_rho - 1  &&
               i_theta  >=0  && i_theta  <= dist->n_theta -1   &&
               i_phi  >=0  && i_phi  <= dist->n_phi - 1  &&
               i_pr   >= 0 && i_pr   <= dist->n_pr - 1   &&
               i_pphi >= 0 && i_pphi <= dist->n_pphi - 1 &&
               i_pz   >= 0 && i_pz   <= dist->n_pz - 1   &&
               i_time >= 0 && i_time <= dist->n_time - 1 &&
               i_q    >= 0 && i_q    <= dist->n_q - 1      ) {
                real weight = p_f->weight[i] * (p_f->time[i] - p_i->time[i]);
                size_t index = dist_rho6D_index(
                    i_rho, i_theta, i_phi, i_pr, i_pphi, i_pz,
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
void dist_rho6D_update_gc(dist_rho6D_data* dist, particle_simd_gc* p_f,
                          particle_simd_gc* p_i) {
    real phi[NSIMD];
    real theta[NSIMD];

    int i_rho[NSIMD];
    int i_theta[NSIMD];
    int i_phi[NSIMD];
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
            real B_dB[12] = {p_f->B_r[i],
                             p_f->B_r_dr[i],
                             p_f->B_r_dphi[i],
                             p_f->B_r_dz[i],
                             p_f->B_phi[i],
                             p_f->B_phi_dr[i],
                             p_f->B_phi_dphi[i],
                             p_f->B_phi_dz[i],
                             p_f->B_z[i],
                             p_f->B_z_dr[i],
                             p_f->B_z_dphi[i],
                             p_f->B_z_dz[i]};
            gctransform_pparmuzeta2prpphipz(p_f->mass[i], p_f->charge[i], B_dB,
                                            p_f->phi[i], p_f->ppar[i],
                                            p_f->mu[i], p_f->zeta[i],
                                            &pr, &pphi, &pz);

            i_rho[i] = floor((p_f->rho[i] - dist->min_rho)
                             / ((dist->max_rho - dist->min_rho)/dist->n_rho));

            phi[i] = fmod(p_f->phi[i], 2*CONST_PI);
            if(phi[i] < 0) {
                phi[i] = phi[i] + 2*CONST_PI;
            }
            i_phi[i] = floor((phi[i] - dist->min_phi)
                             / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            theta[i] = fmod(p_f->theta[i], 2*CONST_PI);
            if(theta[i] < 0) {
                theta[i] = theta[i] + 2*CONST_PI;
            }
            i_theta[i] = floor((theta[i] - dist->min_theta)
                             / ((dist->max_theta - dist->min_theta)
                                / dist->n_theta));

            i_pr[i] = floor((pr - dist->min_pr)
                            / ((dist->max_pr - dist->min_pr) / dist->n_pr));

            i_pphi[i] = floor((pphi - dist->min_pphi)
                              / ((dist->max_pphi - dist->min_pphi)
                                 / dist->n_pphi));

            i_pz[i] = floor((pz - dist->min_pz)
                            / ((dist->max_pz - dist->min_pz) / dist->n_pz));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i]/CONST_E - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_rho[i]  >= 0 && i_rho[i]  <= dist->n_rho - 1  &&
               i_theta[i]  >= 0 && i_theta[i]  <= dist->n_theta -1   &&
               i_phi[i]  >= 0 && i_phi[i]  <= dist->n_phi - 1  &&
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
            size_t index = dist_rho6D_index(
                i_rho[i], i_theta[i], i_phi[i], i_pr[i], i_pphi[i], i_pz[i],
                i_time[i], i_q[i], dist->step_7, dist->step_6, dist->step_5,
                dist->step_4, dist->step_3, dist->step_2, dist->step_1);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}
