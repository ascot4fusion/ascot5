/**
 * @file diag_transcoef.h
 * @brief Header file for diag_transcoef.c.
 *
 * Contains definitions for transport coefficient data structures.
 */
#ifndef DIAG_TRANSCOEF_H
#define DIAG_TRANSCOEF_H

#include "../ascot5.h"
#include "../particle.h"

/**
 * @brief Simple linked list link for storing data points.
 */
typedef struct diag_transcoef_link{
    real rho;      /**< Current radial coordinate (R or rho)                */
    real time;     /**< Current time                                        */
    int pitchsign; /**< Sign of the pitch at this point                     */
    struct diag_transcoef_link* prevlink; /**< Pointer to the previous link */
}diag_transcoef_link;

/**
 * @brief Transport coefficient diagnostics offload data struct.
 */
typedef struct{
    integer Nmrk;  /**< Number of markers in this simulation                  */
    int Navg;      /**< Data is divided into subarrays of length Navg, before
                        taking average value and evaluating K and D           */
    int recordrho; /**< Flag for whether the spatial unit is rho or R.        */
    real interval; /**< Interval at which markers are recorded.               */

}diag_transcoef_offload_data;

/**
 * @brief Transport coefficient diagnostics offload data struct.
 */
typedef struct{
    int Navg;      /**< Data is divided into subarrays of length Navg, before
                        taking average value and evaluating K and D           */
    int recordrho; /**< Flag for whether the spatial unit is rho or R.        */
    real interval; /**< Interval at which markers are recorded.               */
    diag_transcoef_link** datapoints; /**< Temporary data storage             */

    real* id;    /**< Marker ID whose data is stored at this index            */
    real* Kcoef; /**< Calculated drift coefficients                           */
    real* Dcoef; /**< Calculated diffusion coefficients where negative value
                    means coefficients are/were not calculated                */

}diag_transcoef_data;

#pragma omp declare target
void diag_transcoef_init(diag_transcoef_data* data,
                         diag_transcoef_offload_data* offload_data,
                         real* offload_array);

void diag_transcoef_free(diag_transcoef_data* data);

void diag_transcoef_update_fo(diag_transcoef_data* data,
                              particle_simd_fo* p_f, particle_simd_fo* p_i);

void diag_transcoef_update_gc(diag_transcoef_data* data,
                              particle_simd_gc* p_f, particle_simd_gc* p_i);

void diag_transcoef_update_ml(diag_transcoef_data* data,
                              particle_simd_ml* p_f, particle_simd_ml* p_i);

#pragma omp end declare target
#endif
