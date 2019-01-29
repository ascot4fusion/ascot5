/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include "ascot5.h"
#include "error.h"
#include "boozer.h"

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
 *
 * @todo Konsta will write this.
 */
int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array) {

    /* Initialize spline for g and copy the coefficients */
    interp1D_data spline1D;
    err += interp1Dcomp_init(&spline,
                             &(*offload_array)[0],
                             offload_data->npsi, offload_data->psimin,
                             offload_data->psimax,
                             offload_data->psigrid);

    for(int i=0; i < cpernode*npsi; i++) {
        coeff_array[i] = spline1D.c[i];
    }
    interp1Dcomp_free(&spline1D);

    /* Initialize spline for q and copy the coefficients */
    err += interp1Dcomp_init(&spline1D,
                             &(*offload_array)[npsi],
                             offload_data->npsi, offload_data->psimin,
                             offload_data->psimax,
                             offload_data->psigrid);

    for(int i=0; i < cpernode*npsi; i++) {
        coeff_array[cpernode*npsi + i] = spline1D.c[i];
    }
    interp1Dcomp_free(&spline1D);

    /* Initialize spline for I and copy the coefficients */
    err += interp1Dcomp_init(&spline1D,
                             &(*offload_array)[2*npsi],
                             offload_data->npsi, offload_data->psimin,
                             offload_data->psimax,
                             offload_data->psigrid);

    for(int i=0; i < cpernode*npsi; i++) {
        coeff_array[2*cpernode*npsi + i] = spline1D.c[i];
    }
    interp1Dcomp_free(&spline1D);

    /* Initialize spline for delta and copy the coefficients */
    interp1D_data spline2D;
    cpernode = 4;
    err += interp2Dcomp_init(&spline2D,
                             &(*offload_array)[3*npsi],
                             offload_data->npsi, offload_data->ntheta_bzr,
                             offload_data->psimin,
                             offload_data->psimax,
                             offload_data->psigrid,
                             offload_data->thetabzrmin,
                             offload_data->thetabzrmax,
                             offload_data->thetabzrgrid);

    for(int i=0; i < cpernode*npsi*ntheta_bzr; i++) {
        coeff_array[3*cpernode*npsi + i] = spline2D.c[i];
    }
    interp2Dcomp_free(&spline2D);

    /* Initialize spline for nu and copy the coefficients */
    cpernode = 4;
    err += interp2Dcomp_init(&spline2D,
                             &(*offload_array)[(3 + ntheta_bzr)*npsi],
                             offload_data->npsi, offload_data->ntheta_bzr,
                             offload_data->psimin,
                             offload_data->psimax,
                             offload_data->psigrid,
                             offload_data->thetabzrmin,
                             offload_data->thetabzrmax,
                             offload_data->thetabzrgrid);

    for(int i=0; i < cpernode*npsi*ntheta_bzr; i++) {
        coeff_array[(3 + ntheta_bzr)*cpernode*npsi + i] = spline2D.c[i];
    }
    interp2Dcomp_free(&spline2D);

    /* Initialize spline for theta_bzr and copy the coefficients */
    cpernode = 4;
    err += interp2Dcomp_init(&spline2D,
                             &(*offload_array)[(3 + 2*ntheta_bzr)*npsi],
                             offload_data->npsi, offload_data->ntheta_geo,
                             offload_data->psimin,
                             offload_data->psimax,
                             offload_data->psigrid,
                             offload_data->thetageomin,
                             offload_data->thetageomax,
                             offload_data->thetageogrid);

    for(int i=0; i < cpernode*npsi*ntheta_geo; i++) {
        coeff_array[(3 + 2*ntheta_bzr)*cpernode*npsi + i] = spline2D.c[i];
    }
    interp2Dcomp_free(&spline2D);

    /* Initialize spline for theta_geo and copy the coefficients */

    return err;
}

/**
 * @brief Initialize boozer data struct on target
 *
 * @param boozerdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @todo Konsta will write this.
 */
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                 real* offload_array) {
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
