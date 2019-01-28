/**
 * @file boozer.h
 * @brief Header file for boozer.c
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "ascot5.h"
#include "error.h"

/**
 * @brief Boozer parameters that will be offloaded to target
 */
typedef struct {
    int nR;         /**< Number of major radius points for the psi_pol(R,z)
                         map                                                  */
    int nz;         /**< Number of major radius points for the psi_pol(R,z)
                         map                                                  */
    int npsi_pol;   /**< Number of psi_pol values for the coordinate grids    */
    int ntheta_geo; /**< Number of geometric theta angle values               */
    int ntheta;     /**< Number of the Boozer theta values                    */

    int offload_array_length; /**< Number of elements in offload_array        */
} boozer_offload_data;

/**
 * @brief Boozer parameters on the target
 */
typedef struct {
    int nR;          /**< Number of major radius points for the psi_pol(R,z)
                          map                                                 */
    int nz;          /**< Number of major radius points for the psi_pol(R,z)
                          map                                                 */
    int npsi_pol;    /**< Number of psi_pol values for the coordinate grids   */
    int ntheta_geo;  /**< Number of geometric theta angle values              */
    int ntheta;      /**< Number of the Boozer theta values                   */

    real* theta_geo; /**< Geometric theta angle                               */
    real* g;         /**< Toroidal covariant component of the magnetic field  */
    real* q;         /**< q-profile                                           */
    real* I;         /**< poloidal covariant component of the magnetic field  */
    real* delta;     /**< radial covariant component of the magnetic field    */
    real* theta;     /**< Boozer theta angle                                  */
    real* nu;        /**< the nu-function, phi=zeta+nu(psi,theta),
                          phi the cylindrical angle                           */
} boozer_data;

int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array);

#pragma omp declare target
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array);

#pragma omp declare simd uniform(boozerdata)
a5err boozer_eval(boozer_data* boozerdata);

#pragma omp end declare target

#endif
