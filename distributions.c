/**
 * @file distributions.c
 * @brief Distribution collecting and processing functions.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "distributions.h"
#include "B_GS.h"
#include "particle.h"

/** @brief Internal function calculating the index in the histogram array */
#pragma omp declare target
unsigned long dist_rzvv_index(int i_r, int i_z,
                    int i_vpara, int i_vperp,
                    int n_z, int n_vpara, int n_vperp) {
    return   i_r     * (n_z * n_vpara * n_vperp)
           + i_z     * (n_vpara * n_vperp)
           + i_vpara * (n_vperp)
           + i_vperp;
}
#pragma omp end declare target

/**
 * @brief Initialize distribution
 *
 * Fill the offload struct with distribution parameters and 
 * allocates and empties the offload array.
 */
void dist_rzvv_init_offload(dist_rzvv_offload_data* offload_data,
                            real** offload_array) {
    offload_data->offload_array_length = offload_data->n_r * offload_data->n_z
        * offload_data->n_vpara * offload_data->n_vperp;

    *offload_array = (real*) malloc(sizeof(real)
                                    * offload_data->offload_array_length);
    if(*offload_array == 0) {
        printf("Histogram memory allocation failed\n");
        exit(1);
    }

    int i;
    for(i = 0; i < offload_data->offload_array_length; i++) {
        (*offload_array)[i] = 0.0;
    }
}

void dist_rzvv_free_offload(dist_rzvv_offload_data* offload_data,
                            real** offload_array) {
    offload_data->n_r = 0;
    offload_data->min_r = 0;
    offload_data->max_r = 0;
    offload_data->n_z = 0;
    offload_data->min_z = 0;
    offload_data->max_z = 0;
    offload_data->n_vpara = 0;
    offload_data->min_vpara = 0;
    offload_data->max_vpara = 0;
    offload_data->n_vperp = 0;
    offload_data->min_vperp = 0;
    offload_data->max_vperp = 0;
    free(*offload_array);
    *offload_array = NULL;
}

void dist_rzvv_init(dist_rzvv_data* dist_data,
                    dist_rzvv_offload_data* offload_data,
                    real* offload_array) {
    dist_data->n_r = offload_data->n_r;
    dist_data->min_r = offload_data->min_r;
    dist_data->max_r = offload_data->max_r;
    dist_data->n_z = offload_data->n_z;
    dist_data->min_z = offload_data->min_z;
    dist_data->max_z = offload_data->max_z;
    dist_data->n_vpara = offload_data->n_vpara;
    dist_data->min_vpara = offload_data->min_vpara;
    dist_data->max_vpara = offload_data->max_vpara;
    dist_data->n_vperp = offload_data->n_vperp;
    dist_data->min_vperp = offload_data->min_vperp;
    dist_data->max_vperp = offload_data->max_vperp;
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
void dist_rzvv_update_fo(dist_rzvv_data* dist, particle_simd_fo* p, real dt) {
    real r[NSIMD];
    real vpara[NSIMD];
    real vperp[NSIMD];
    real weight[NSIMD];
    int i_r[NSIMD];
    int i_z[NSIMD];
    int i_vpara[NSIMD];
    int i_vperp[NSIMD];
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        i_r[i] = floor((p->r[i] - dist->min_r)
            / ((dist->max_r - dist->min_r)/dist->n_r));
        if(i_r[i] < 0)
            i_r[i] = 0;
        if(i_r[i] > dist->n_r - 1)
            i_r[i] = dist->n_r - 1;

        i_z[i] = floor((p->z[i] - dist->min_z)
                / ((dist->max_z - dist->min_z) / dist->n_z));
        if(i_z[i] < 0)
            i_z[i] = 0;
        if(i_z[i] > dist->n_z - 1)
            i_z[i] = dist->n_z - 1;

	vpara[i] = (p->rdot[i] * p->B_r[i] + (p->phidot[i] * p->r[i]) * p->B_phi[i]
		    + p->zdot[i] * p->B_z[i])
	    / sqrt(p->B_r[i]*p->B_r[i]+p->B_phi[i]*p->B_phi[i]
		   + p->B_z[i]*p->B_z[i]);
        i_vpara[i] = floor((vpara[i] - dist->min_vpara)
                / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));
        if(i_vpara[i] < 0)
            i_vpara[i] = 0;
        if(i_vpara[i] > dist->n_vpara - 1)
            i_vpara[i] = dist->n_vpara - 1;

        vperp[i] = sqrt(p->rdot[i]*p->rdot[i] + (p->phidot[i]*p->phidot[i]*p->r[i]*p->r[i])
                        + p->zdot[i]*p->zdot[i] - vpara[i]*vpara[i]);
        i_vperp[i] = floor((vperp[i] - dist->min_vperp)
                / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));
        if(i_vperp[i] < 0)
            i_vperp[i] = 0;
        if(i_vperp[i] > dist->n_vperp - 1)
            i_vperp[i] = dist->n_vperp - 1;
        
        if(p->running[i]) {
            weight[i] = p->weight[i] * dt;
        }
        else {
            weight[i] = 0;
        }
    }

    for(i = 0; i < NSIMD; i++) {
        unsigned long index = dist_rzvv_index(i_r[i], i_z[i], i_vpara[i],
                                              i_vperp[i], dist->n_z,
                                              dist->n_vpara, dist->n_vperp);

        #pragma omp atomic 
        dist->histogram[index] += weight[i];
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
void dist_rzvv_update_gc(dist_rzvv_data* dist, particle_simd_gc* p, real dt) {
    real vperp[NSIMD];
    real weight[NSIMD];
    int i_r[NSIMD];
    int i_z[NSIMD];
    int i_vpara[NSIMD];
    int i_vperp[NSIMD];
    int i;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        i_r[i] = floor((p->r[i] - dist->min_r)
            / ((dist->max_r - dist->min_r)/dist->n_r));
        if(i_r[i] < 0)
            i_r[i] = 0;
        if(i_r[i] > dist->n_r - 1)
            i_r[i] = dist->n_r - 1;

        i_z[i] = floor((p->z[i] - dist->min_z)
                / ((dist->max_z - dist->min_z) / dist->n_z));
        if(i_z[i] < 0)
            i_z[i] = 0;
        if(i_z[i] > dist->n_z - 1)
            i_z[i] = dist->n_z - 1;

                   ;
        i_vpara[i] = floor((p->vpar[i] - dist->min_vpara)
                / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));
        if(i_vpara[i] < 0)
            i_vpara[i] = 0;
        if(i_vpara[i] > dist->n_vpara - 1)
            i_vpara[i] = dist->n_vpara - 1;

        vperp[i] = sqrt(2 * sqrt(p->B_r[i]*p->B_r[i]+p->B_phi[i]*p->B_phi[i]
                               +p->B_z[i]*p->B_z[i])
                        * p->mu[i] / p->mass[i]);
        i_vperp[i] = floor((vperp[i] - dist->min_vperp)
                / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));
        if(i_vperp[i] < 0)
            i_vperp[i] = 0;
        if(i_vperp[i] > dist->n_vperp - 1)
            i_vperp[i] = dist->n_vperp - 1;
        
        if(p->running[i]) {
            weight[i] = p->weight[i] * dt;
        }
        else {
            weight[i] = 0;
        }
    }

    for(i = 0; i < NSIMD; i++) {
        unsigned long index = dist_rzvv_index(i_r[i], i_z[i], i_vpara[i],
                                              i_vperp[i], dist->n_z,
                                              dist->n_vpara, dist->n_vperp);

        #pragma omp atomic 
        dist->histogram[index] += weight[i];
    }
}

/**
 * @brief Print rz-distribution
 *
 * This function prints the rz marginal distribution, summing over velocity
 * dimensions.
 *
 * @param dist pointer to distribution parameter struct
 * @param histogram pointer to histogram
 */
void dist_rzvv_print_rz(dist_rzvv_offload_data* dist, real* histogram) {
    int i, j, k, l;
    for(i = 0; i < dist->n_r; i++) {
        for(j = 0; j < dist->n_z; j++) {
            real sum = 0;
            for(k = 0; k < dist->n_vpara; k++) {
                for(l = 0; l < dist->n_vperp; l++) {
                    unsigned long index = dist_rzvv_index(i, j, k, l,
                                      dist->n_z, dist->n_vpara, dist->n_vperp);
                    sum += histogram[index];
                }
            }
            printf("%le ", sum);
        }
        printf("\n");
    } 
}

/**
 * @brief Print vv-distribution
 *
 * This function prints the vv marginal distribution, summing over spatial
 * dimensions.
 *
 * @param dist pointer to distribution parameter struct
 * @param histogram pointer to histogram
 */
void dist_rzvv_print_vv(dist_rzvv_offload_data* dist, real* histogram) {
    int i, j, k, l;
    for(i = 0; i < dist->n_vpara; i++) {
        for(j = 0; j < dist->n_vperp; j++) {
            real sum = 0;
            for(k = 0; k < dist->n_r; k++) {
                for(l = 0; l < dist->n_z; l++) {
                    unsigned long index = dist_rzvv_index(k, l, i, j,
                                      dist->n_z, dist->n_vpara, dist->n_vperp);
                    sum += histogram[index];
                }
            }
            printf("%le ", sum);
        }
        printf("\n");
    } 
}

void dist_rzvv_sum(dist_rzvv_offload_data* dist1, real* array1, real* array2) {
    int i;

    for(i=0; i < dist1->offload_array_length; i++) {
        array1[i] += array2[i];
    }
}


