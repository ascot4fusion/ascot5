/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include <stdlib.h>
#include "ascot5.h"
#include "consts.h"
#include "error.h"
#include "boozer.h"
#include "spline/interp1D.h"
#include "spline/interp1Dcomp.h"

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

    real psi_grid = (offload_data->psi_max - offload_data->psi_min)
        / (offload_data->npsi - 1);
    int npsi       = offload_data->npsi;
    int nR         = offload_data->nR;
    int nz         = offload_data->nz;
    int ntheta_bzr = offload_data->ntheta_bzr;
    int ntheta_geo = offload_data->ntheta_geo;

    /* Allocate array for storing coefficients (which later replaces the
       offload array) */
    offload_data->offload_array_length =
        npsi * (3*2 + 4*ntheta_geo + 4*2*ntheta_bzr) + 4*nR*nz;
    real* coeff_array = (real*)malloc(offload_data->offload_array_length
                                      * sizeof(real));

    int cpernode = 2; // Number of coefficients in one node.
    /* Initialize spline for g and copy the coefficients */
    interp1D_data spline1D;
    err += interp1Dcomp_init(&spline1D,
                             &(*offload_array)[0],
                             offload_data->npsi, offload_data->psi_min,
                             offload_data->psi_max,
                             psi_grid);

    for(int i=0; i < cpernode*npsi; i++) {
        coeff_array[i] = spline1D.c[i];
    }
    interp1Dcomp_free(&spline1D);

    /* Initialize spline for q and copy the coefficients */
    err += interp1Dcomp_init(&spline1D,
                             &(*offload_array)[npsi],
                             offload_data->npsi, offload_data->psi_min,
                             offload_data->psi_max,
                             psi_grid);

    for(int i=0; i < cpernode*npsi; i++) {
        coeff_array[cpernode*npsi + i] = spline1D.c[i];
    }
    interp1Dcomp_free(&spline1D);

    /* Initialize spline for I and copy the coefficients */
    err += interp1Dcomp_init(&spline1D,
                             &(*offload_array)[2*npsi],
                             offload_data->npsi, offload_data->psi_min,
                             offload_data->psi_max,
                             psi_grid);

    for(int i=0; i < cpernode*npsi; i++) {
        coeff_array[2*cpernode*npsi + i] = spline1D.c[i];
    }
    interp1Dcomp_free(&spline1D);

    /* Initialize spline for delta and copy the coefficients */
    err += interp2Dcomp_init_coeff(&(coeff_array[(3 + 2*ntheta_bzr)*cpernode*npsi]),
                                   &(*offload_array)[3*npsi],
                                   offload_data->npsi, offload_data->ntheta_bzr,
                                   NATURALBC, PERIODICBC,
                                   offload_data->psi_min,
                                   offload_data->psi_max,
                                   0, CONST_2PI);

    /* Initialize spline for nu and copy the coefficients */
    err += interp2Dcomp_init_coeff(&(coeff_array[(3 + 2*ntheta_bzr)*cpernode*npsi]),
                                   &(*offload_array)[(3 + ntheta_bzr)*npsi],
                                   offload_data->npsi, offload_data->ntheta_bzr,
                                   NATURALBC, PERIODICBC,
                                   offload_data->psi_min,
                                   offload_data->psi_max,
                                   0, CONST_2PI);

    /* Initialize spline for theta_bzr and copy the coefficients */
    err += interp2Dcomp_init_coeff(&(coeff_array[(3 + 2*ntheta_bzr)*cpernode*npsi]),
                                   &(*offload_array)[(3 + 2*ntheta_bzr)*npsi],
                                   offload_data->npsi, offload_data->ntheta_geo,
                                   NATURALBC, NATURALBC,
                                   offload_data->R_min,
                                   offload_data->R_max,
                                   offload_data->z_min,
                                   offload_data->z_max);

    /* Initialize spline for theta_geo and copy the coefficients */

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

    interp2Dcomp_init_spline(&boozerdata->delta,
                             &(offload_array[0]),
                             offload_data->npsi,
                             offload_data->ntheta_bzr,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             0, CONST_2PI);

    interp2Dcomp_init_spline(&boozerdata->nu,
                             &(offload_array[0]),
                             offload_data->npsi,
                             offload_data->ntheta_bzr,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             0, CONST_2PI);

    interp2Dcomp_init_spline(&boozerdata->theta_bzr,
                             &(offload_array[0]),
                             offload_data->npsi,
                             offload_data->ntheta_geo,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             0, CONST_2PI);

    interp2Dcomp_init_spline(&boozerdata->theta_geo,
                             &(offload_array[0]),
                             offload_data->nR,
                             offload_data->nz,
                             NATURALBC, NATURALBC,
                             offload_data->R_min,
                             offload_data->R_max,
                             offload_data->z_min,
                             offload_data->z_max);
}

/**
 * @brief Transform cylindrical coordinates to Boozer coordinates.
 *
 * @todo This is just a dummy.
 *
 * @param ptz Boozer coordinates as [psi, theta, zeta].
 *
 * @return zero on success
 */
a5err boozer_cyl2booz(real ptz[3], real r, real phi, real z,
                      boozer_data* boozerdata) {
    return 0;
}

/**
 * @brief Transform Boozer coordinates to cylindrical coordinates.
 *
 * @todo This is just a dummy.
 *
 * @param rz cylindrical coordinates as [R, z].
 *
 * @return zero on success
 */
a5err boozer_booz2cyl(real rz[2], real psi, real theta, real zeta,
                      boozer_data* boozerdata) {
    return 0;
}

/**
 * @brief Evaluate Boozer coordinates and gradients on a given location.
 *
 * @todo This is just a dummy.
 *
 * The values are stored in the given array as:
 * - ptz_dptz[0]  = psi
 * - ptz_dptz[1]  = dpsi/dR
 * - ptz_dptz[2]  = dpsi/dphi
 * - ptz_dptz[3]  = dpsi/dz
 * - ptz_dptz[4]  = theta
 * - ptz_dptz[5]  = dtheta/dR
 * - ptz_dptz[6]  = dtheta/dphi
 * - ptz_dptz[7]  = dtheta/dz
 * - ptz_dptz[8]  = zeta
 * - ptz_dptz[9]  = dzeta/dR
 * - ptz_dptz[10] = dzeta/dphi
 * - ptz_dptz[11] = dzeta/dz
 *
 * @param ptz_dptz evaluated Boozer coordinates and their gradients.
 *
 * @return zero on success
 */
a5err boozer_eval_gradients(real ptz_dptz[12], real r, real phi, real z,
                            boozer_data* boozerdata) {
    return 0;
}
