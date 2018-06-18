/**
 * @file dist_rho5D.c
 * @brief Distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "dist_rho5D.h"
#include "particle.h"

/** @brief Internal function calculating the index in the histogram array */
#pragma omp declare target
unsigned long dist_rho5D_index(int i_rho, int i_pol, int i_phi, int i_vpara,
                            int i_vperp, int n_phi, int n_pol, int n_vpara,
                            int n_vperp) {
    return   i_rho   * (n_pol * n_phi * n_vpara * n_vperp)
           + i_pol   * (n_phi * n_vpara * n_vperp)
           + i_phi   * (n_vpara * n_vperp)
           + i_vpara * (n_vperp)
           + i_vperp;
}
#pragma omp end declare target

/**
 * @brief Frees the offload data
 */
void dist_rho5D_free_offload(dist_rho5D_offload_data* offload_data) {
	offload_data->n_rho = 0;
	offload_data->min_rho = 0;
	offload_data->max_rho = 0;
	offload_data->n_pol = 0;
	offload_data->min_pol = 0;
	offload_data->max_pol = 0;
	offload_data->n_phi = 0;
	offload_data->min_phi = 0;
	offload_data->max_phi = 0;
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
void dist_rho5D_init(dist_rho5D_data* dist_data, dist_rho5D_offload_data* offload_data,
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
void dist_rho5D_update_fo(dist_rho5D_data* dist, particle_simd_fo* p_f,
                       particle_simd_fo* p_i) {
    real phi[NSIMD];
    real pol[NSIMD];
    real vpara[NSIMD];
    real vperp[NSIMD];

	int i_rho[NSIMD];
	int i_pol[NSIMD];
	int i_phi[NSIMD];
    int i_vpara[NSIMD];
    int i_vperp[NSIMD];

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

            if(i_rho[i] >= 0      && i_rho[i] <= dist->n_rho - 1
               && i_phi[i] >= 0   && i_phi[i] <= dist->n_phi - 1
               && i_pol[i] >= 0   && i_pol[i] <= dist->n_pol - 1
               && i_vpara[i] >= 0 && i_vpara[i] <= dist->n_vpara - 1
               && i_vperp[i] >= 0 && i_vperp[i] <= dist->n_vperp - 1) {
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
            unsigned long index = dist_rho5D_index(i_rho[i], i_pol[i], i_phi[i],
                                                i_vpara[i], i_vperp[i],
                                                dist->n_pol, dist->n_phi,
                                                dist->n_vpara, dist->n_vperp);
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
void dist_rho5D_update_gc(dist_rho5D_data* dist, particle_simd_gc* p_f,
                       particle_simd_gc* p_i) {
    real phi[NSIMD];
    real pol[NSIMD];
    real vperp[NSIMD];

    int i_rho[NSIMD];
    int i_phi[NSIMD];
    int i_pol[NSIMD];
    int i_vpara[NSIMD];
    int i_vperp[NSIMD];

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

            i_vpara[i] = floor((p_f->vpar[i] - dist->min_vpara)
                       / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));

            vperp[i] = sqrt(2 * sqrt(p_f->B_r[i]*p_f->B_r[i]
                                     +p_f->B_phi[i]*p_f->B_phi[i]
                                     +p_f->B_z[i]*p_f->B_z[i])
                            * p_f->mu[i] / p_f->mass[i]);
            i_vperp[i] = floor((vperp[i] - dist->min_vperp)
                       / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));

            if(i_rho[i] >= 0        && i_rho[i] <= dist->n_rho - 1
               && i_phi[i] >= 0   && i_phi[i] <= dist->n_phi - 1
               && i_pol[i] >= 0     && i_pol[i] <= dist->n_pol - 1
               && i_vpara[i] >= 0 && i_vpara[i] <= dist->n_vpara - 1
               && i_vperp[i] >= 0 && i_vperp[i] <= dist->n_vperp - 1) {
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
            unsigned long index = dist_rho5D_index(i_rho[i], i_pol[i], i_phi[i],
                                            i_vpara[i], i_vperp[i], dist->n_pol,
                                            dist->n_phi, dist->n_vpara,
                                            dist->n_vperp);

            #pragma omp atomic
            dist->histogram[index] += weight[i];
        }
    }
}


void dist_rho5D_sum(int start, int stop, real* array1, real* array2) {
    for(int i = start; i < stop; i++) {
        array1[i] += array2[i];
    }
}
