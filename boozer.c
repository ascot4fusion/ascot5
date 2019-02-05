/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "error.h"
#include "boozer.h"
#include "spline/interp.h"

/**
 * @brief Load Boozer data and prepare parameters for offload.
 *
 * This function fills the boozer offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * The offload data struct should be fully initialized before calling this
 * function and offload array should hold the input data in order
 * [g, q, I, delta, nu, theta_bzr, theta_geo]. This function fits splines to
 * input data, reallocates the offload array and stores spline coefficients
 * there.
 *
 * Multidimensional arrays must be stored as
 * - delta(psi_i, thetabzr_j)     = array[j*npsi + i]
 * - nu(psi_i, thetabzr_j)        = array[j*npsi + i]
 * - theta_bzr(psi_i, thetageo_j) = array[j*npsi + i]
 * - theta_geo(R_i, z_j)          = array[j*nR + i]
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 */
int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array) {

    int err = 0;

    /* Size of 1D and 2D input data arrays */
    int psisize   = offload_data->npsi;
    int thgeosize = offload_data->nr   * offload_data->nz;
    int thbzrsize = offload_data->npsi * offload_data->ntheta_geo;
    int psithsize = offload_data->npsi * offload_data->ntheta_bzr;

    /* Grid limits for theta_bzr and theta_geo grids*/
    int THETAMIN = 0;
    int THETAMAX = CONST_2PI;

    /* Allocate array for storing coefficients (which later replaces the
       offload array) */
    real* coeff_array = (real*)malloc( ( 3*psisize*NSIZE_COMP1D
                                         + 2*psithsize*NSIZE_COMP2D
                                         + thbzrsize*NSIZE_COMP2D
                                         + thgeosize*NSIZE_COMP2D)
                                       * sizeof(real) );

    /* Evaluate and store coefficients */

    /* g */
    err += interp1Dcomp_init_coeff(
            &coeff_array[0],
            &(*offload_array)[0],
            offload_data->npsi, NATURALBC,
            offload_data->psi_min, offload_data->psi_max);

    /* q */
    err += interp1Dcomp_init_coeff(
            &coeff_array[NSIZE_COMP1D * psisize],
            &(*offload_array)[psisize],
            offload_data->npsi, NATURALBC,
            offload_data->psi_min, offload_data->psi_max);

    /* I */
    err += interp1Dcomp_init_coeff(
            &coeff_array[2 * NSIZE_COMP1D * psisize],
            &(*offload_array)[2 * psisize],
            offload_data->npsi, NATURALBC,
            offload_data->psi_min, offload_data->psi_max);

    /* delta */
    err += interp2Dcomp_init_coeff(
            &coeff_array[3 * NSIZE_COMP1D * psisize],
            &(*offload_array)[3 * psisize],
            offload_data->npsi, offload_data->ntheta_bzr,
            NATURALBC, PERIODICBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN, THETAMAX);

    /* nu */
    err += interp2Dcomp_init_coeff(
            &coeff_array[3 * NSIZE_COMP1D * psisize
                         + psithsize * NSIZE_COMP2D],
            &(*offload_array)[3 * psisize + psithsize],
            offload_data->npsi, offload_data->ntheta_bzr,
            NATURALBC, PERIODICBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN, THETAMAX);

    /* theta_bzr */
    err += interp2Dcomp_init_coeff(
            &coeff_array[3 * NSIZE_COMP1D * psisize
                         + 2 * psithsize * NSIZE_COMP2D],
            &(*offload_array)[3 * psisize + 2 * psithsize],
            offload_data->npsi, offload_data->ntheta_geo,
            NATURALBC, PERIODICBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN, THETAMAX);

    /* theta_geo */
    err += interp2Dcomp_init_coeff(
            &coeff_array[3 * NSIZE_COMP1D * psisize
                         + 2 * psithsize * NSIZE_COMP2D
                         + thbzrsize * NSIZE_COMP2D],
            &(*offload_array)[3 * psisize + 2 * psithsize + thbzrsize],
            offload_data->nr, offload_data->nz,
            NATURALBC, NATURALBC,
            offload_data->r_min, offload_data->r_max,
            offload_data->z_min, offload_data->z_max);

    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = 3*psisize*NSIZE_COMP1D
                                         + 2*psithsize*NSIZE_COMP2D
                                         + thbzrsize*NSIZE_COMP2D
                                         + thgeosize*NSIZE_COMP2D;

    return err;
}

/**
 * @brief Initialize boozer data struct on target
 *
 * @param boozerdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array) {

    /* Grid limits for theta_bzr and theta_geo grids*/
    int THETAMIN = 0;
    int THETAMAX = CONST_2PI;

    /* Size of 1D and 2D input data arrays */
    int psisize   = offload_data->npsi * NSIZE_COMP1D;
    int psithsize = offload_data->npsi * offload_data->ntheta_bzr*NSIZE_COMP2D;
    int thbzrsize = offload_data->npsi * offload_data->ntheta_geo*NSIZE_COMP2D;

    /* Initialize splines */

    interp1Dcomp_init_spline(&boozerdata->g,
                             &(offload_array[0 * psisize]),
                             offload_data->npsi, NATURALBC,
                             offload_data->psi_min, offload_data->psi_max);

    interp1Dcomp_init_spline(&boozerdata->q,
                             &(offload_array[1 * psisize]),
                             offload_data->npsi, NATURALBC,
                             offload_data->psi_min, offload_data->psi_max);

    interp1Dcomp_init_spline(&boozerdata->I,
                             &(offload_array[2 * psisize]),
                             offload_data->npsi, NATURALBC,
                             offload_data->psi_min, offload_data->psi_max);

    interp2Dcomp_init_spline(&boozerdata->delta,
                             &(offload_array[3 * psisize + 0 * psithsize]),
                             offload_data->npsi,
                             offload_data->ntheta_bzr,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN, THETAMAX);

    interp2Dcomp_init_spline(&boozerdata->nu,
                             &(offload_array[3 * psisize + 1 * psithsize]),
                             offload_data->npsi,
                             offload_data->ntheta_bzr,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN, THETAMAX);

    interp2Dcomp_init_spline(&boozerdata->theta_bzr,
                             &(offload_array[3 * psisize + 2 * psithsize
                                             + 0 * thbzrsize]),
                             offload_data->npsi,
                             offload_data->ntheta_geo,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN, THETAMAX);

    interp2Dcomp_init_spline(&boozerdata->theta_geo,
                             &(offload_array[3 * psisize + 2 * psithsize
                                             + 1 * thbzrsize]),
                             offload_data->nr,
                             offload_data->nz,
                             NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->z_min,
                             offload_data->z_max);
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void boozer_free_offload(boozer_offload_data* offload_data,
                         real** offload_array) {
    free(*offload_array);
}

/**
 * @brief Evaluate Boozer angular coordinates and partial derivatives.
 *
 * The values are stored in the given array as:
 * - ptz_dptz[0] = theta
 * - ptz_dptz[1] = dtheta/dR
 * - ptz_dptz[2] = dtheta/dphi
 * - ptz_dptz[3] = dtheta/dz
 * - ptz_dptz[4] = zeta
 * - ptz_dptz[5] = dzeta/dR
 * - ptz_dptz[6] = dzeta/dphi
 * - ptz_dptz[7] = dzeta/dz
 *
 * @param thetazeta evaluated Boozer angular coordinates and their gradients.
 *
 * @return zero on success
 */
a5err boozer_eval_thetazeta(real thetazeta[8], real r, real phi, real z,
                            real psi_dpsi[4], boozer_data* boozerdata) {
    a5err err = 0;
    int interperr = 0;

    /* geometrical theta and derivatives as helper variables */
    real thgeo[6];
    interperr += interp2Dcomp_eval_df(thgeo, &boozerdata->theta_geo, r, z);

    /* theta and derivatives */
    real theta[6];
    interperr += interp2Dcomp_eval_df(theta, &boozerdata->theta_bzr,
                                      psi_dpsi[0], thgeo[0]);
    thetazeta[0] = theta[0];
    thetazeta[1] = psi_dpsi[1] * theta[1] + thgeo[1] * theta[2];
    thetazeta[2] = 0;
    thetazeta[3] = psi_dpsi[3] * theta[1] + thgeo[3] * theta[2];

    /* zeta and derivatives */
    real q[3];
    interperr += interp1Dcomp_eval_df(q, &boozerdata->q, psi_dpsi[0]);

    thetazeta[4] = fmod(phi + theta[0]*q[0], CONST_2PI);
    thetazeta[5] = thetazeta[1]*q[0] + theta[0]*q[1]*psi_dpsi[1];
    thetazeta[6] = 1;
    thetazeta[7] = thetazeta[3]*q[0] + theta[0]*q[1]*psi_dpsi[3];

    if(interperr) {
        err = 1;
    }

    return err;
}
