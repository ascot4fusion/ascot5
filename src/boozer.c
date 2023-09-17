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
 * @brief Load Boozer data and prepare parameters for offload.
 *
 * This function fills the boozer offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * The offload data struct should be fully initialized before calling this
 * function and offload array should hold the input data in order
 * [psi, nu, theta_bzr]. This function fits splines to input data, reallocates
 * the offload array and stores spline coefficients there.
 *
 * Multidimensional arrays must be stored as
 * - nu(psi_i, thetabzr_j)        = array[j*npsi + i]
 * - theta_bzr(psi_i, thetageo_j) = array[j*npsi + i]
 * - psi(R_i, z_j)                = array[j*nR + i]
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 */
int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array) {

    int err = 0;

    /* Size of 2D grids and 1D data */
    int nusize      = offload_data->npsi * offload_data->ntheta;
    int thetasize   = offload_data->npsi * offload_data->nthetag;
    int contoursize = offload_data->nrzs;

    /* Grid limits for theta_bzr and theta_geo grids */
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = (4.0*CONST_2PI)/(offload_data->nthetag-2*4.0 -1);

    /* Allocate array for storing coefficients (which later replaces the
       offload array) and contour points */
    real* coeff_array = (real*)malloc( ( ( nusize + thetasize)
                                         * NSIZE_COMP2D + 2*contoursize )
                                       * sizeof(real) );

    /* Evaluate and store coefficients */

    /* nu */
    err += interp2Dcomp_init_coeff(
            &coeff_array[0],
            &(*offload_array)[0],
            offload_data->npsi, offload_data->ntheta,
            NATURALBC, PERIODICBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN, THETAMAX);

    /* theta_bzr */
    err += interp2Dcomp_init_coeff(
            &coeff_array[nusize * NSIZE_COMP2D],
            &(*offload_array)[nusize],
            offload_data->npsi, offload_data->nthetag,
            NATURALBC, NATURALBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN-padding, THETAMAX+padding);

    for(int i = 0; i < contoursize; i++) {
        coeff_array[(nusize + thetasize)*NSIZE_COMP2D + i] =
            (*offload_array)[nusize + thetasize + i];
        coeff_array[(nusize + thetasize)*NSIZE_COMP2D + contoursize + i] =
            (*offload_array)[nusize + thetasize + contoursize + i];
    }

    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = (nusize + thetasize)
                                         * NSIZE_COMP2D + 2 * contoursize;

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nBoozer input\n");
    print_out(VERBOSE_IO, "psi grid: n = %4.d min = %3.3f max = %3.3f\n",
              offload_data->npsi,
              offload_data->psi_min, offload_data->psi_max);
    print_out(VERBOSE_IO, "thetageo grid: n = %4.d\n", offload_data->nthetag);
    print_out(VERBOSE_IO, "thetabzr grid: n = %4.d\n", offload_data->ntheta);

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

    boozerdata->psi_min  = offload_data->psi_min;
    boozerdata->psi_max  = offload_data->psi_max;

    /* Grid limits for theta_bzr and theta_geo grids*/
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = (4.0*CONST_2PI)/(offload_data->nthetag-2*4.0-1);

    /* Size of 1D and 2D input data arrays */
    int nusize      = offload_data->npsi * offload_data->ntheta * NSIZE_COMP2D;
    int thetasize   = offload_data->npsi * offload_data->nthetag * NSIZE_COMP2D;
    int contoursize = offload_data->nrzs;

    /* Initialize splines */
    interp2Dcomp_init_spline(&boozerdata->nu_psitheta,
                             &(offload_array[0]),
                             offload_data->npsi,
                             offload_data->ntheta,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN, THETAMAX);

    interp2Dcomp_init_spline(&boozerdata->theta_psithetageom,
                             &(offload_array[nusize]),
                             offload_data->npsi,
                             offload_data->nthetag,
                             NATURALBC, NATURALBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN-padding, THETAMAX+padding);

    boozerdata->rs   = &(offload_array[nusize + thetasize]);
    boozerdata->zs   = &(offload_array[nusize + thetasize + contoursize]);
    boozerdata->nrzs = offload_data->nrzs;
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
