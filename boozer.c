/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "consts.h"
#include "math.h"
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
    int rzsize      = offload_data->nr   * offload_data->nz;
    int nusize      = offload_data->npsi * offload_data->ntheta;
    int thetasize   = offload_data->npsi * offload_data->nthetag;
    int contoursize = offload_data->nrzs;

    /* Grid limits for theta_bzr and theta_geo grids */
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = (4.0*CONST_2PI)/(offload_data->nthetag-2*4.0 -1);

    /* Allocate array for storing coefficients (which later replaces the
       offload array) and contour points */
    real* coeff_array = (real*)malloc( ( ( rzsize + nusize + thetasize)
                                         * NSIZE_COMP2D + 2*contoursize )
                                       * sizeof(real) );

    /* Evaluate and store coefficients */

    /* psi */
    err += interp2Dcomp_init_coeff(
            &coeff_array[0],
            &(*offload_array)[0],
            offload_data->nr, offload_data->nz,
            NATURALBC, NATURALBC,
            offload_data->r_min, offload_data->r_max,
            offload_data->z_min, offload_data->z_max);

    /* nu */
    err += interp2Dcomp_init_coeff(
            &coeff_array[rzsize * NSIZE_COMP2D],
            &(*offload_array)[rzsize],
            offload_data->npsi, offload_data->ntheta,
            NATURALBC, PERIODICBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN, THETAMAX);

    /* theta_bzr */
    err += interp2Dcomp_init_coeff(
            &coeff_array[(rzsize + nusize) * NSIZE_COMP2D],
            &(*offload_array)[rzsize + nusize],
            offload_data->npsi, offload_data->nthetag,
            NATURALBC, NATURALBC,
            offload_data->psi_min, offload_data->psi_max,
            THETAMIN-padding, THETAMAX+padding);

    for(int i = 0; i < contoursize; i++) {
        coeff_array[(rzsize + nusize + thetasize)*NSIZE_COMP2D + i] =
            (*offload_array)[rzsize + nusize + thetasize + i];
        coeff_array[(rzsize + nusize + thetasize)*NSIZE_COMP2D + contoursize + i] =
            (*offload_array)[rzsize + nusize + thetasize + contoursize + i];
    }

    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = (rzsize + nusize + thetasize)
                                         * NSIZE_COMP2D + 2* contoursize;

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

    boozerdata->psi_inner = offload_data->psi_inner;
    boozerdata->psi_outer = offload_data->psi_outer;
    boozerdata->r0        = offload_data->r0;
    boozerdata->z0        = offload_data->z0;
    boozerdata->nrzs      = offload_data->nrzs;

    /* Grid limits for theta_bzr and theta_geo grids*/
    real THETAMIN = 0;
    real THETAMAX = CONST_2PI;
    real padding = (4.0*CONST_2PI)/(offload_data->nthetag-2*4.0-1);

    /* Size of 1D and 2D input data arrays */
    int rzsize      = offload_data->nr * offload_data->nz * NSIZE_COMP2D;
    int nusize      = offload_data->npsi * offload_data->ntheta * NSIZE_COMP2D;
    int thetasize   = offload_data->npsi * offload_data->nthetag * NSIZE_COMP2D;
    int contoursize = offload_data->nrzs;

    /* Initialize splines */

    interp2Dcomp_init_spline(&boozerdata->psi_rz,
                             &(offload_array[0]),
                             offload_data->nr,
                             offload_data->nz,
                             NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->z_min,
                             offload_data->z_max);

    interp2Dcomp_init_spline(&boozerdata->nu_psitheta,
                             &(offload_array[rzsize]),
                             offload_data->npsi,
                             offload_data->ntheta,
                             NATURALBC, PERIODICBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN, THETAMAX);

    interp2Dcomp_init_spline(&boozerdata->theta_psithetageom,
                             &(offload_array[rzsize + nusize]),
                             offload_data->npsi,
                             offload_data->nthetag,
                             NATURALBC, NATURALBC,
                             offload_data->psi_min,
                             offload_data->psi_max,
                             THETAMIN-padding, THETAMAX+padding);

    boozerdata->rs = &(offload_array[rzsize + nusize + thetasize]);
    boozerdata->zs = &(offload_array[rzsize + nusize + thetasize + contoursize]);
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
 * @brief Evaluates the normalized psi and derivatives from
 * given psi and its derivatives
 *
 * @param psi a vector [psi,dpsi/dR, dpsi/dphi, dpsi/dz]
 * @param psin a normalized psi
 */
a5err boozer_eval_psinormalized(real psi[4], real psin[4],
                                boozer_data* boozerdata){
    a5err err = 0;

    psin[0] = ( boozerdata-> psi_outer - psi[0] )
        / ( boozerdata->psi_outer - boozerdata->psi_inner );
    psin[1] = psi[1] / ( boozerdata->psi_outer - boozerdata->psi_inner );
    psin[2] = psi[2] / ( boozerdata->psi_outer - boozerdata->psi_inner );
    psin[3] = psi[3] / ( boozerdata->psi_outer - boozerdata->psi_inner );

    return err;
}

/**
 * @brief Evaluate Boozer coordinates and partial derivatives
 *
 * the modified vectors
 *
 * - isinside /= 0 , the point (r,phi,z) is inside the grid
 *            == 0 , the point (r,phi,z) is outside the grid
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
 * @param thetazeta evaluated Boozer angular coordinates and their gradients.
 * @param isinside a flag indicating whether the queried point was inside
 *        boozer grid
 * @param r R coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param boozerdata pointer to boozerdata
 *
 * @return zero on success
 */
a5err boozer_eval_psithetazeta(real psithetazeta[12], int* isinside,
                               real r, real phi, real z,
                               boozer_data* boozerdata) {
    a5err err = 0;
    int interperr = 0;


    /* winding number to test whether we are inside the plasma */
    if(math_point_in_polygon(r, z, boozerdata->rs, boozerdata->zs,
                             boozerdata->nrzs)) {

        /* get the psi value and check that it is within the psi grid (the grid
           does not extend all the way to the axis) */
        real psi[6];
        interperr += interp2Dcomp_eval_df(psi, &boozerdata->psi_rz, r, z);

        if(psi[0] >= boozerdata->psi_inner && psi[0] <= boozerdata->psi_outer) {

            /* update the flag, and we are good to go */
            isinside[0]=1;

            /* geometrical theta */
            real thgeo;
            thgeo = fmod( atan2(z-boozerdata->z0,r-boozerdata->r0) + CONST_2PI,
                          CONST_2PI);

            /* boozer theta and derivatives */
            real theta[6];
            interperr += interp2Dcomp_eval_df(
                theta, &boozerdata->theta_psithetageom, psi[0], thgeo);

            /* boozer nu function and derivatives */
            real nu[6];
            interperr += interp2Dcomp_eval_df(
                nu, &boozerdata->nu_psitheta, psi[0], theta[0]);

            /* set up data for returning the requested values */

            /* psi and derivatives */
            psithetazeta[0]=psi[0]; /* psi       */
            psithetazeta[1]=psi[1]; /* dpsi_dr   */
            psithetazeta[2]=0;      /* dpsi_dphi */
            psithetazeta[3]=psi[2]; /* dpsi_dz   */

            /* helpers */
            real asq;
            asq = (r - boozerdata->r0) * (r - boozerdata->r0)
                + (z - boozerdata->z0) * (z - boozerdata->z0);
            real dthgeo_dr;
            dthgeo_dr=-(z-boozerdata->z0)/asq;
            real dthgeo_dz;
            dthgeo_dz=(r-boozerdata->r0)/asq;

            /* theta and derivatives */
            psithetazeta[4]=theta[0];                          /* theta       */
            psithetazeta[5]=theta[1]*psi[1]+theta[2]*dthgeo_dr;/* dtheta_dr   */
            psithetazeta[6]=0;                                 /* dtheta_dphi */
            psithetazeta[7]=theta[1]*psi[2]+theta[2]*dthgeo_dz;/* dtheta_dz   */

            /* zeta and derivatives */
            psithetazeta[8]=fmod(phi-nu[0], CONST_2PI);          /* zeta      */
            psithetazeta[9]=-nu[1]*psi[1]-nu[2]*psithetazeta[5]; /* dzeta_dR  */
            psithetazeta[10]=1.0;                                /* dzeta_dphi*/
            psithetazeta[11]=-nu[1]*psi[2]-nu[2]*psithetazeta[7];/* dzeta_dz  */
        }
        else {
            /* This would mean that (r,z) is in the very center of the plasma */
            isinside[0]=0;
        }
    }
    else {
        /* The winding number test indicates that (r,z) is outside the outermost
           boozer grid psi contour. */
        isinside[0]=0;
    }


    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_BOOZER );
    }

    return err;
}


/**
 * @brief Evaluate normalized poloidal flux rho
 *
 * @param rho pointer where rho value will be stored
 * @param psi poloidal flux value from which rho is evaluated
 * @param boozerdata pointer to boozer data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err boozer_eval_rho(real* rho, real psi, boozer_data* boozerdata) {

    /* Check that the values seem valid */
    real delta = (boozerdata->psi_outer - boozerdata->psi_inner);
    if( (boozerdata->psi_outer-psi) / delta < 0 ) {
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_BOOZER );
    }

    rho[0] = sqrt( (boozerdata->psi_outer - psi ) / delta );
    return 0;
}
