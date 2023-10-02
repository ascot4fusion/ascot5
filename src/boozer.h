/**
 * @file boozer.h
 * @brief Header file for boozer.c
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "spline/interp.h"

/**
 * @brief offload data for maps between boozer and cylindrical coordinates
 */
typedef struct {
    int  npsi;     /**< Number of psi grid points for other fields            */
    real psi_min;  /**< Minimum psi in other fields                           */
    real psi_max;  /**< Maximum psi in other fields                           */
    int  ntheta;   /**< number of boozer theta grid points                    */
    int  nthetag;  /**< number of geometric theta grid points                 */
    int  nrzs;     /**< number of elements in rs and zs                       */
    int  offload_array_length; /**< Number of elements in offload_array       */
} boozer_offload_data;

/**
 * @brief Boozer parameters on the target
 */
typedef struct {
    real psi_min;  /**< Minimum psi in other fields                           */
    real psi_max;  /**< Maximum psi in other fields                           */
    real* rs;  /**< R points of outermost poloidal psi-surface contour,
                    nrzs elements, the first and last points are the same     */
    real* zs;  /**< z points of outermost poloidal psi-surface contour,
                    nrzs elements, the first and last points are the same     */
    int  nrzs; /**< number of elements in rs and zs                           */
    interp2D_data nu_psitheta; /**< the nu-function, phi=zeta+nu(psi,theta),
                                    with phi the cylindrical angle            */
    interp2D_data theta_psithetageom; /**< boozer_theta(psi,thetag)           */
} boozer_data;

int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array);
void boozer_free_offload(boozer_offload_data* offload_data,
                         real** offload_array);

#pragma omp declare target
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array);

#pragma omp declare simd uniform(Bdata, boozerdata)
a5err boozer_eval_psithetazeta(real psithetazeta[12], int* isinside, real r,
                               real phi, real z, B_field_data* Bdata,
                               boozer_data* boozerdata);

#pragma omp end declare target

#endif
