/**
 * @file boozer.h
 * @brief Header file for boozer.c
 */
#ifndef BOOZER_H
#define BOOZER_H

#include "ascot5.h"
#include "error.h"
#include "spline/interp1D.h"
#include "spline/interp1Dcomp.h"
#include "spline/interp2Dcomp.h"

/**
 * @brief Boozer parameters that will be offloaded to target
 *
 * The actual data consists of g(psi), q(psi), I(psi), delta(psi, theta_bzr),
 * nu(psi, theta_bzr), theta_bzr(psi, theta_geo), and theta_geo(R, z)values.
 * This struct holds enough information to construct the grids where data is
 * given.
 */
typedef struct {
    int nR;         /**< Number R grid points               */
    real R_min;     /**< Minimum R grid point               */
    real R_max;     /**< Maximum R grid point               */
    int nz;         /**< Number z grid points               */
    real z_min;     /**< Minimum R grid point               */
    real z_max;     /**< Minimum R grid point               */
    int npsi;       /**< Number psi grid points             */
    real psi_min;   /**< Minimum psi grid point             */
    real psi_max;   /**< Maximum psi grid point             */
    int ntheta_geo; /**< Number of geometric theta points   */
    int ntheta_bzr; /**< Number of boozer theta grid points */

    int offload_array_length; /**< Number of elements in offload_array        */
} boozer_offload_data;

/**
 * @brief Boozer parameters on the target
 */
typedef struct {
    interp2D_data theta_geo; /**< Geometric theta angle                       */
    interp2D_data theta_bzr; /**< Boozer theta angle                          */
    interp1D_data g; /**< Toroidal covariant component of the magnetic field  */
    interp1D_data q; /**< q-profile                                           */
    interp1D_data I; /**< poloidal covariant component of the magnetic field  */
    interp2D_data delta; /**< radial covariant component of the magnetic field*/
    interp2D_data nu;    /**< the nu-function, phi=zeta+nu(psi,theta),
                              phi the cylindrical angle                       */
} boozer_data;

int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array);

#pragma omp declare target
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array);

#pragma omp declare simd uniform(boozerdata)
a5err boozer_cyl2booz(real ptz[3], real r, real phi, real z,
                      boozer_data* boozerdata);
#pragma omp declare simd uniform(boozerdata)
a5err boozer_booz2cyl(real rz[2], real psi, real theta, real zeta,
                      boozer_data* boozerdata);
#pragma omp declare simd uniform(boozerdata)
a5err boozer_eval_gradients(real ptz_dptz[12], real r, real phi, real z,
                            boozer_data* boozerdata);

#pragma omp end declare target

#endif
