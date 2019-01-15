/**
 * @file dist_rho6D.c
 * @brief rho 6D distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "dist_rho6D.h"
#include "../particle.h"

/**
 * @brief Internal function calculating the index in the histogram array
 */
#pragma omp declare target
unsigned long dist_rho6D_index(int i_rho, int i_pol, int i_phi, int i_vr, 
                               int i_vphi, int i_vz, int i_time, int i_q,
                               int n_pol, int n_phi, int n_vr, int n_vphi,
                               int n_vz, int n_time, int n_q) {
    return i_rho  * (n_pol * n_phi * n_vr * n_vphi * n_vz * n_time * n_q)
        + i_pol  * (n_phi * n_vr * n_vphi * n_vz * n_time * n_q)
        + i_phi  * (n_vr * n_vphi * n_vz * n_time * n_q)
        + i_vr   * (n_vphi * n_vz * n_time * n_q)
        + i_vphi * (n_vz * n_time * n_q)
        + i_vz   * (n_time * n_q)
        + i_time * (n_q)
        + i_q;
}
#pragma omp end declare target

/**
 * @brief Frees the offload data
 */
void dist_rho6D_free_offload(dist_rho6D_offload_data* offload_data) {
    offload_data->n_rho = 0;
    offload_data->min_rho = 0;
    offload_data->max_rho = 0;
    offload_data->n_pol = 0;
    offload_data->min_pol = 0;
    offload_data->max_pol = 0;
    offload_data->n_phi = 0;
    offload_data->min_phi = 0;
    offload_data->max_phi = 0;
    offload_data->n_vr = 0;
    offload_data->min_vr = 0;
    offload_data->max_vr = 0;
    offload_data->n_vphi = 0;
    offload_data->min_vphi = 0;
    offload_data->max_vphi = 0;
    offload_data->n_vz = 0;
    offload_data->min_vz = 0;
    offload_data->max_vz = 0;
}

/**
 * @brief Initializes distribution from offload data
 */
void dist_rho6D_init(dist_rho6D_data* dist_data, dist_rho6D_offload_data* offload_data,
                     real* offload_array) {
    dist_data->n_rho = offload_data->n_rho;
    dist_data->min_rho = offload_data->min_rho;
    dist_data->max_rho = offload_data->max_rho;

    dist_data->n_pol = offload_data->n_pol;
    dist_data->min_pol = offload_data->min_pol;
    dist_data->max_pol = offload_data->max_pol;

    dist_data->n_phi = offload_data->n_phi;
    dist_data->min_phi = offload_data->min_phi;
    dist_data->max_phi = offload_data->max_phi;

    dist_data->n_vr   = offload_data->n_vr;
    dist_data->min_vr = offload_data->min_vr;
    dist_data->max_vr = offload_data->max_vr;

    dist_data->n_vphi   = offload_data->n_vphi;
    dist_data->min_vphi = offload_data->min_vphi;
    dist_data->max_vphi = offload_data->max_vphi;

    dist_data->n_vz   = offload_data->n_vz;
    dist_data->min_vz = offload_data->min_vz;
    dist_data->max_vz = offload_data->max_vz;

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
 * @param p_i pointer to SIMD particle struct at the beginning of time step
 * @param p_f pointer to SIMD particle struct at the end of time step
 */
void dist_rho6D_update_fo(dist_rho6D_data* dist, particle_simd_fo* p_f,
                          particle_simd_fo* p_i) {
    real phi[NSIMD];
    real pol[NSIMD];

    int i_rho[NSIMD];
    int i_pol[NSIMD];
    int i_phi[NSIMD];
    int i_vr[NSIMD];
    int i_vphi[NSIMD];
    int i_vz[NSIMD];
    int i_time[NSIMD];
    int i_q[NSIMD];

    int ok[NSIMD];
    real weight[NSIMD];

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            i_rho[i] = floor((p_f->rho[i] - dist->min_rho)
                             / ((dist->max_rho - dist->min_rho)/dist->n_rho));

            phi[i] = fmod(p_f->phi[i], 2*CONST_PI);
            if(phi[i] < 0) {
                phi[i] = phi[i] + 2*CONST_PI;
            }
            i_phi[i] = floor((phi[i] - dist->min_phi)
                             / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            pol[i] = fmod(p_f->pol[i], 2*CONST_PI);
            if(pol[i] < 0) {
                pol[i] = pol[i] + 2*CONST_PI;
            }
            i_pol[i] = floor((pol[i] - dist->min_pol)
                             / ((dist->max_pol - dist->min_pol) / dist->n_pol));

            i_vr[i] = floor((p_f->rdot[i] - dist->min_vr)
                            / ((dist->max_vr - dist->min_vr) / dist->n_vr));

            i_vphi[i] = floor((p_f->phidot[i]*p_f->r[i] - dist->min_vphi)
                              / ((dist->max_vphi - dist->min_vphi) / dist->n_vphi));

            i_vz[i] = floor((p_f->zdot[i] - dist->min_vz)
                            / ((dist->max_vz - dist->min_vz) / dist->n_vz));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i] - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_rho[i]  >= 0 && i_rho[i]  <= dist->n_rho - 1  &&
               i_pol[i]  >=0  && i_pol[i]  <= dist->n_pol -1   &&
               i_phi[i]  >=0  && i_phi[i]  <= dist->n_phi - 1  &&
               i_vr[i]   >= 0 && i_vr[i]   <= dist->n_vr - 1   &&
               i_vphi[i] >= 0 && i_vphi[i] <= dist->n_vphi - 1 &&
               i_vz[i]   >= 0 && i_vz[i]   <= dist->n_vz - 1   &&
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
            unsigned long index = dist_rho6D_index(i_rho[i], i_pol[i], i_phi[i],
                                                   i_vr[i], i_vphi[i], i_vz[i],
                                                   i_time[i], i_q[i], 
                                                   dist->n_pol, dist->n_phi,
                                                   dist->n_vr, dist->n_vphi,
                                                   dist->n_vz, dist->n_time,
                                                   dist->n_q);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
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
    real pol[NSIMD];

    int i_rho[NSIMD];
    int i_pol[NSIMD];
    int i_phi[NSIMD];
    int i_vr[NSIMD];
    int i_vphi[NSIMD];
    int i_vz[NSIMD];
    int i_time[NSIMD];
    int i_q[NSIMD];

    int ok[NSIMD];
    real weight[NSIMD];

    particle_state p_s;

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {
            particle_gc_to_state(p_f, i, &p_s, NULL);

            i_rho[i] = floor((p_s.rho - dist->min_rho)
                             / ((dist->max_rho - dist->min_rho)/dist->n_rho));

            phi[i] = fmod(p_s.phi, 2*CONST_PI);
            if(phi[i] < 0) {
                phi[i] = phi[i] + 2*CONST_PI;
            }
            i_phi[i] = floor((phi[i] - dist->min_phi)
                             / ((dist->max_phi - dist->min_phi)/dist->n_phi));

            pol[i] = fmod(p_s.pol, 2*CONST_PI);
            if(pol[i] < 0) {
                pol[i] = pol[i] + 2*CONST_PI;
            }
            i_pol[i] = floor((pol[i] - dist->min_pol)
                             / ((dist->max_pol - dist->min_pol) / dist->n_pol));

            i_vr[i] = floor((p_s.rdot - dist->min_vr)
                            / ((dist->max_vr - dist->min_vr) / dist->n_vr));

            i_vphi[i] = floor((p_s.phidot*p_s.r - dist->min_vphi)
                              / ((dist->max_vphi - dist->min_vphi) / dist->n_vphi));

            i_vz[i] = floor((p_s.zdot - dist->min_vz)
                            / ((dist->max_vz - dist->min_vz) / dist->n_vz));

            i_time[i] = floor((p_f->time[i] - dist->min_time)
                          / ((dist->max_time - dist->min_time) / dist->n_time));

            i_q[i] = floor((p_f->charge[i] - dist->min_q)
                           / ((dist->max_q - dist->min_q) / dist->n_q));

            if(i_rho[i]  >= 0 && i_rho[i]  <= dist->n_rho - 1  &&
               i_pol[i]  >= 0 && i_pol[i]  <= dist->n_pol -1   &&
               i_phi[i]  >= 0 && i_phi[i]  <= dist->n_phi - 1  &&
               i_vr[i]   >= 0 && i_vr[i]   <= dist->n_vr - 1   &&
               i_vphi[i] >= 0 && i_vphi[i] <= dist->n_vphi - 1 &&
               i_vz[i]   >= 0 && i_vz[i]   <= dist->n_vz - 1   &&
               i_time[i] >= 0 && i_time[i] <= dist->n_time - 1 &&
               i_q[i]    >= 0 && i_q[i]    <= dist->n_q - 1      ) {
                ok[i] = 1;
                weight[i] = p_s.weight * (p_s.time - p_i->time[i]);
            }
            else {
                ok[i] = 0;
            }
        }
    }

    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i] && ok[i]) {
            unsigned long index = dist_rho6D_index(i_rho[i], i_pol[i], i_phi[i],
                                                   i_vr[i], i_vphi[i], i_vz[i], 
                                                   i_time[i], i_q[i],
                                                   dist->n_pol, dist->n_phi,
                                                   dist->n_vr, dist->n_vphi,
                                                   dist->n_vz, dist->n_time,
                                                   dist->n_q);
            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}
