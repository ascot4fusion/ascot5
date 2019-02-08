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
 * @brief Evaluates the normalized psi and derivatives from 
 * given psi and its derivatives
 * 
 * @param psi a vector [psi,dpsi/dR, dpsi/dphi, dpsi/dz]
 * @param psin a normalized psi
 */
a5err boozer_eval_psinormalized(real psi[4], real psin[4],boozer_data* boozerdata){
    a5err err = 0;
    psin[0]=(psi[0]-boozerdata->psi_inner)/(boozerdata->psi_outer-boozerdata->psi_inner);
    psin[1]=psi[1]/(boozerdata->psi_outer-boozerdata->psi_inner);
    psin[2]=psi[2]/(boozerdata->psi_outer-boozerdata->psi_inner);
    psin[3]=psi[3]/(boozerdata->psi_outer-boozerdata->psi_inner);
}

/**
 * @brief Evaluate Boozer coordinates and partial derivatives 
 *
 *
 * the modified vectors
 *
 * - isinside /= 0 , the point (r,phi,z) is inside the grid
 *            == 0 , the point (r,phi,z) is outside the grid
 *
 * - psithetazeta[0] = psi
 * - psithetazeta[1] = dpsi/dR
 * - psithetazeta[2] = dpsi/dphi
 * - psithetazeta[3] = dpsi/dz
 * - psithetazeta[4] = theta
 * - psithetazeta[5] = dtheta/dR
 * - psithetazeta[6] = dtheta/dphi
 * - psithetazeta[7] = dtheta/dz
 * - psithetazeta[8] = zeta
 * - psithetazeta[9] = dzeta/dR
 * - psithetazeta[10] = dzeta/dphi
 * - psithetazeta[11] = dzeta/dz
 *
 * @param thetazeta evaluated Boozer angular coordinates and their gradients.
 *
 * @return zero on success
 */
a5err boozer_eval_psithetazeta(real psithetazeta[12], int* isinside, real r,
			       real phi, real z, boozer_data* boozerdata) {
    a5err err = 0;
    int interperr = 0;


    /* winding number to test whether we are inside the plasma */
    if(math_point_in_polygon(r, z, boozerdata->rs, boozerdata->zs,
			     boozerdata->nrzs)) {
	
	/* get the psi value and check that it is within the psi grid (the grid 
	   does not extend all the way to the axis) */ 
	real psi[6]; 
	interperr += interp2Dcomp_eval_df(psi, &boozerdata->psi_rz,r,z);
	if(psi[0] >= boozerdata->psimin && psi[0] <= boozerdata->psimax) {

	    /* update the flag, and we are good to go */
	    isinside[0]=1;
	    
	    /* geometrical theta */
	    real thgeo;
	    thgeo = fmod(atan2(z-boozerdata->z0,r-boozerdata->r0),CONSTANT_2PI);
	    
	    /* boozer theta and derivatives */
	    real theta[6]; 
	    interperr += interp2Dcomp_eval_df(
		theta,&boozerdata->theta_psithetageom,psi[0],thgeo);
	    
	    /* boozer nu function and derivatives */
	    real nu[6];
	    interperr += interp2dcomp_eval_df(
		nu,&boozerdata->nu_psitheta,psi[0],theta[0]);
	    
	    /* set up data for returning the requested values */
	    
	    /* psi and derivatives */
	    psithetazeta[0]=psi[0]; /* psi */
	    psithetazeta[1]=psi[1]; /* dpsi_dr */
	    psithetazeta[2]=0; /* dpsi_dphi */
	    psithetazeta[3]=psi[2]; /* dpsi_dz */
	    
	    /* helpers */
	    real asq; 
	    asq=(r-boozerdata->r0)*(r-boozerdata->r0)+(z-boozerdata->z0)*(z-boozerdata->z0);
	    real dthgeo_dr;
	    dthgeo_dr=-(z-boozerdata->z0)/asq;
	    real dthgeo_dz;
	    dthgeo_dz=(r-boozerdata->r0)/asq;
	    
	    /* theta and derivatives */
	    psithetazeta[4]=theta[0]; /* theta */
	    psithetazeta[5]=theta[1]*psi[1]+theta[2]*dthgeo_dr; /* dtheta_dr */
	    psithetazeta[6]=0; /* dtheta_dphi */
	    psithetazeta[7]=theta[1]*psi[2]+theta[2]*dthgeo_dz; /* dtheta_dz */

	    /* zeta and derivatives */
	    psithetazeta[8]=phi-nu[0]; /* zeta */
	    psithetazeta[9]=-nu[1]*psi[1]-nu[2]*psithetazeta[5]; /* dzeta_dR */
	    psithetazeta[10]=1.0; /* dzeta_dphi */
	    psithetazeta[11]=-nu[1]*psi[2]-nu[2]*psithetazeta[7]; /* dzeta_dz */
	}
	else {
	    /* This whould mean that (r,z) is in the very center of the plasma */
	    isinside[0]=0;
	}
    }
    else {
	/* The winding number test indicates that (r,z) is outside the outermost 
	   boozer grid psi contour. */
	isinside[0]=0;
    }


    if(interperr) {
        err = 1;
    }

    return err;
}
