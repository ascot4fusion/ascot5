/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include <stdlib.h>
#include <math.h>
#include "print.h"
#include "ascot5.h"
#include "consts.h"
#include "math.h"
#include "error.h"
#include "B_field.h"
#include "boozer.h"
#include "spline/interp.h"

/**
 * @brief Initialize boozer coordinate transformation
 *
 * Multidimensional arrays must be stored as
 * - nu(psi_i, theta_j)     = array[j*npsi + i]
 * - theta(psi_i, thetag_j) = array[j*npsi + i]
 *
 * @param data pointer to the data struct
 * @param npsi Number of psi grid points in `nu` and `theta` data
 * @param psi_min minimum value in the psi grid
 * @param psi_max maximum value in the psi grid
 * @param ntheta number of boozer theta grid points in `nu` data
 * @param nthetag number of geometric theta grid points in `theta` data
 * @param nu the difference between cylindrical angle phi and toroidal boozer
 *           coordinate zeta, phi = zeta + nu [rad]
 * @param theta the boozer poloidal angle [rad]
 * @param nrzs the number of elements in `rs` and `zs`
 * @param rs separatrix contour R coordinates [m]
 * @param zs separatrix contour z coordinates [m]
 *
 * @return zero if initialization succeeded.
 */
int boozer_init(boozer_data* data, int npsi, real psi_min, real psi_max,
                int ntheta, int nthetag, real* nu, real* theta,
                int nrzs, real* rs, real* zs) {

    int err = 0;
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = ( 4.0*CONST_2PI ) / ( nthetag - 2*4.0 - 1 );
    data->psi_min = psi_min;
    data->psi_max = psi_max;

    err = interp2Dcomp_setup(&data->nu_psitheta, nu, npsi, ntheta,
                             NATURALBC, PERIODICBC, psi_min, psi_max,
                             THETAMIN, THETAMAX);
    err = interp2Dcomp_setup(&data->theta_psithetageom, theta, npsi, nthetag,
                             NATURALBC, NATURALBC, psi_min, psi_max,
                             THETAMIN-padding, THETAMAX+padding);

    data->nrzs = nrzs;
    for(int i = 0; i < nrzs; i++) {
        data->rs[i] = rs[i];
        data->zs[i] = zs[i];
    }

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nBoozer input\n");
    print_out(VERBOSE_IO, "psi grid: n = %4.d min = %3.3f max = %3.3f\n",
              npsi, psi_min, psi_max);
    print_out(VERBOSE_IO, "thetageo grid: n = %4.d\n", nthetag);
    print_out(VERBOSE_IO, "thetabzr grid: n = %4.d\n", ntheta);

    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void boozer_free(boozer_data* data) {
    free(data->rs);
    free(data->zs);
    free(data->nu_psitheta.c);
    free(data->theta_psithetageom.c);
}

/**
 * @brief Evaluate Boozer coordinates and partial derivatives
 *
 * The output vector has the following elements:
 *
 * - psithetazeta[0]  = psi
 * - psithetazeta[1]  = dpsi/dR
 * - psithetazeta[2]  = dpsi/dphi
 * - psithetazeta[3]  = dpsi/dz
 * - psithetazeta[4]  = theta
 * - psithetazeta[5]  = dtheta/dR
 * - psithetazeta[6]  = dtheta/dphi
 * - psithetazeta[7]  = dtheta/dz
 * - psithetazeta[8]  = zeta
 * - psithetazeta[9]  = dzeta/dR
 * - psithetazeta[10] = dzeta/dphi
 * - psithetazeta[11] = dzeta/dz
 *
 * @param psithetazeta evaluated Boozer coordinates and their gradients.
 * @param isinside a flag indicating whether the queried point was inside
 *        boozer grid
 * @param r R coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data
 * @param boozerdata pointer to boozerdata
 *
 * @return zero on success
 */
a5err boozer_eval_psithetazeta(real psithetazeta[12], int* isinside,
                               real r, real phi, real z, B_field_data* Bdata,
                               boozer_data* boozerdata) {
    a5err err = 0;
    int interperr = 0;

    /* Winding number to test whether we are inside the plasma (and not in the
       private plasma region) */
    isinside[0]=0;
    if(math_point_in_polygon(r, z, boozerdata->rs, boozerdata->zs,
                             boozerdata->nrzs)) {
        /* Get the psi value and check that it is within the psi grid (the grid
           does not extend all the way to the axis) Use t = 0.0 s */
        real psi[4], rho[2], r0, z0;
        err = B_field_eval_psi_dpsi(psi, r, phi, z, 0.0, Bdata);
        if(!err) {
            err = B_field_eval_rho(rho, psi[0], Bdata);
        }
        if(!err) {
            real rz[2];
            err = B_field_get_axis_rz(rz, Bdata, phi);
            r0 = rz[0];
            z0 = rz[1];
        }

        if(!err && psi[0] < boozerdata->psi_max &&
           psi[0] > boozerdata->psi_min) {

            /* Update the flag, and we are good to go */
            isinside[0]=1;

            /* Geometrical theta */
            real thgeo;
            thgeo = fmod( atan2(z-z0,r-r0) + CONST_2PI,
                          CONST_2PI);

            /* Boozer theta and derivatives */
            real theta[6];
            interperr += interp2Dcomp_eval_df(
                theta, &boozerdata->theta_psithetageom, psi[0], thgeo);

            /* Boozer nu function and derivatives */
            real nu[6];
            interperr += interp2Dcomp_eval_df(
                nu, &boozerdata->nu_psitheta, psi[0], theta[0]);

            /* Set up data for returning the requested values */

            /* Psi and derivatives */
            psithetazeta[0]=psi[0]; /* psi       */
            psithetazeta[1]=psi[1]; /* dpsi_dr   */
            psithetazeta[2]=0;      /* dpsi_dphi */
            psithetazeta[3]=psi[3]; /* dpsi_dz   */

            /* Helpers */
            real asq;
            asq = (r - r0) * (r - r0)
                + (z - z0) * (z - z0);
            real dthgeo_dr;
            dthgeo_dr=-(z-z0)/asq;
            real dthgeo_dz;
            dthgeo_dz=(r-r0)/asq;

            /* Theta and derivatives */
            psithetazeta[4]=theta[0];                          /* theta       */
            psithetazeta[5]=theta[1]*psi[1]+theta[2]*dthgeo_dr;/* dtheta_dr   */
            psithetazeta[6]=0;                                 /* dtheta_dphi */
            psithetazeta[7]=theta[1]*psi[3]+theta[2]*dthgeo_dz;/* dtheta_dz   */

            /* Zeta and derivatives */
            psithetazeta[8]=fmod(phi+nu[0], CONST_2PI);          /* zeta      */
            psithetazeta[9]=nu[1]*psi[1]+nu[2]*psithetazeta[5];  /* dzeta_dR  */
            psithetazeta[10]=1.0;                                /* dzeta_dphi*/
            psithetazeta[11]=nu[1]*psi[3]+nu[2]*psithetazeta[7]; /* dzeta_dz  */

            /* Make sure zeta is between [0, 2pi]*/
            psithetazeta[8]=fmod(psithetazeta[8] + CONST_2PI, CONST_2PI);
        }
    }

    if(!err && interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_BOOZER );
    }

    return err;
}
