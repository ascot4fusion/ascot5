/**
 * @file interp3Dcomp.c
 * @brief Tricubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <stdio.h> /* Needed for printf debugging purposes */
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "interp3Dcomp.h"
#include "spline1Dcomp.h"

/**
 * @brief Calculate tricubic spline interpolation coefficients for scalar 3D data
 *
 * This function calculates the tricubic spline interpolation coefficients for
 * the given data and stores them in the data struct. Compact  cofficients are
 * calculated directly.
 *
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 * @param f 3D data to be interpolated
 * @param n_r number of data points in the r direction
 * @param n_phi number of data points in the phi direction
 * @param n_z number of data points in the z direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 * @param phi_min minimum value of the phi axis
 * @param phi_max maximum value of the phi axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 */
int interp3Dcomp_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
		       real r_min, real r_max, real r_grid,
		       real phi_min, real phi_max, real phi_grid,
		       real z_min, real z_max, real z_grid) {

    int err = 0;

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_phi = n_phi;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    str->phi_min = phi_min;
    str->phi_max = phi_max;
    str->phi_grid = phi_grid;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;
    str->c = malloc(n_phi*n_z*n_r*8*sizeof(real));

    /* Declare and allocate the needed variables */
    int i_r;                                    /**< index for r variable */
    int i_phi;                                  /**< index for phi variable */
    int i_z;                                    /**< index for z variable */
    real* f_r = malloc(n_r*sizeof(real));       /**< Temporary array for data along r */
    real* f_phi = malloc(n_phi*sizeof(real));   /**< Temp array for data along phi */
    real* f_z = malloc(n_z*sizeof(real));       /**< Temp array for data along z */
    real* c_r = malloc(n_r*2*sizeof(real));     /**< Temp array for coefficients along r */
    real* c_phi = malloc(n_phi*2*sizeof(real)); /**< Temp array for coefs along phi */
    real* c_z = malloc(n_z*2*sizeof(real));     /**< Temp array for coefficients along z */

    if(f_r == NULL || f_z == NULL || c_r == NULL || c_z == NULL) {
	err = 1;
    }
    else {

	/* Tricubic spline volume coefficients: For i_r, i_phi and i_z the 8 coefficients
	   are [f, frr, fphiphi, fzz, frrphiphi, frrzz, fphiphizz, frrphiphizz].
	   Note how we account for normalized grid. */
	/* Bicubic spline surface over rz-grid for each phi */
	for(i_phi=0; i_phi<n_phi; i_phi++) {
	    /* Cubic spline along r for each z to get frr */
	    for(i_z=0; i_z<n_z; i_z++) {
		for(i_r=0; i_r<n_r; i_r++) {
		    f_r[i_r] = f[i_phi*n_z*n_r+i_z*n_r+i_r];
		}
		spline1Dcomp(f_r,n_r,0,c_r);
		for(i_r=0; i_r<n_r; i_r++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8] = c_r[i_r*2];
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+1] = c_r[i_r*2+1]/(r_grid*r_grid);
		}
	    }
	
	    /* Two cubic splines along z for each r using f and frr */
	    for(i_r=0; i_r<n_r; i_r++) {
		/* fzz */
		for(i_z=0; i_z<n_z; i_z++) {
		    f_z[i_z] =  f[i_phi*n_z*n_r+i_z*n_r+i_r];
		}
		spline1Dcomp(f_z,n_z,0,c_z);
		for(i_z=0; i_z<n_z; i_z++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+3] = c_z[i_z*2+1]/(z_grid*z_grid);
		}
		/* frrzz */
		for(i_z=0; i_z<n_z; i_z++) {
		    f_z[i_z] = str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+1];
		}
		spline1Dcomp(f_z,n_z,0,c_z);
		for(i_z=0; i_z<n_z; i_z++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+5] = c_z[i_z*2+1]/(z_grid*z_grid);
		}
	    }

	}

	/* Cubic spline along phi for each rz-pair to find the compact coefficients
	   of the tricubic spline volume */
	for(i_z=0; i_z<n_z; i_z++) {
	    for(i_r=0; i_r<n_r; i_r++) {
		/* fphiphi */
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    f_phi[i_phi] = f[i_phi*n_z*n_r+i_z*n_r+i_r];
		}
		spline1Dcomp(f_phi,n_phi,1,c_phi);
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+2] = c_phi[i_phi*2+1]/(phi_grid*phi_grid);
		}
		/* frrphiphi */
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    f_phi[i_phi] = str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+1];
		}
		spline1Dcomp(f_phi,n_phi,1,c_phi);
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+4] = c_phi[i_phi*2+1]/(phi_grid*phi_grid);
		}
		/* fphiphizz */
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    f_phi[i_phi] = str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+3];
		}
		spline1Dcomp(f_phi,n_phi,1,c_phi);
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+6] = c_phi[i_phi*2+1]/(phi_grid*phi_grid);
		}
		/* frrphiphizz */
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    f_phi[i_phi] = str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+5];
		}
		spline1Dcomp(f_phi,n_phi,1,c_phi);
		for(i_phi=0; i_phi<n_phi; i_phi++) {
		    str->c[i_phi*n_z*n_r*8+i_z*n_r*8+i_r*8+7] = c_phi[i_phi*2+1]/(phi_grid*phi_grid);
		}
	    }
	}
    }

    /* Free allocated memory */
    free(f_r);
    free(f_phi);
    free(f_z);
    free(c_r);
    free(c_phi);
    free(c_z);

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * tricubic spline interpolation coefficients of the compact form.
 * 
 * @todo Seg fault if in last phi-sector becaause of compact eval +1-coefs
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
integer interp3Dcomp_eval_B(real* B, interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI + phi;}

    int i_r = (r-str->r_min)/str->r_grid;     /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dri = 1.0-dr;
    real dr3 = dr*dr*dr-dr;
    real dri3 = (1.0-dr)*(1.0-dr)*(1.0-dr)-(1.0-dr);
    real rg2 = str->r_grid*str->r_grid;       /**< Square of cell length in r direction */
    int i_phi = (phi-str->phi_min)/str->phi_grid; /**< index for phi variable */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
									   coordinate in
									   current cell */
    real dphii = 1.0-dphi;
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphii3 = (1.0-dphi)*(1.0-dphi)*(1.0-dphi)-(1.0-dphi);
    real phig2 = str->phi_grid*str->phi_grid; /**< Square of cell length in phi dir */
    int i_z = (z-str->z_min)/str->z_grid;     /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */
    real dzi = 1.0-dz;
    real dz3 = dz*dz*dz-dz;
    real dzi3 = (1.0-dz)*(1.0-dz)*(1.0-dz)-(1.0-dz);
    real zg2 = str->z_grid*str->z_grid;       /**< Square of cell length in z direction */
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8; /**< Index jump to cell */
    int r1 = 8;                               /**< Index jump one r forward */
    int phi1 = str->n_z*str->n_r*8;           /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
	phi1 = -(str->n_phi-1)*phi1;          /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r*8;                      /**< Index jump one z forward */
	   
    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {

	*B = (
	    dzi*(
		dri*(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])+
		dr*(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
	    +dz*(
		dri*(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])+
		dr*(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
	    +rg2/6*(
		dzi*(
		    dri3*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
		    dr3*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
		+dz*(
		    dri3*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
		    dr3*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
	    +phig2/6*(
		dzi*(
		    dri*(dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2])+
		    dr*(dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
		+dz*(
		    dri*(dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2])+
		    dr*(dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
	    +zg2/6*(
		dzi3*(
		    dri*(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])+
		    dr*(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
		+dz3*(
		    dri*(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])+
		    dr*(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
	    +rg2*phig2/36*(
		dzi*(
		    dri3*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
		    dr3*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
		+dz*(
		    dri3*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
		    dr3*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
	    +rg2*zg2/36*(
		dzi3*(
		    dri3*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
		    dr3*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
		+dz3*(
		    dri3*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
		    dr3*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
	    +phig2*zg2/36*(
		dzi3*(
		    dri*(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])+
		    dr*(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
		+dz3*(
		    dri*(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])+
		    dr*(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
	    +rg2*phig2*zg2/216*(
		dzi3*(
		    dri3*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
		    dr3*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
		+dz3*(
		    dri3*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
		    dr3*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * tricubic spline interpolation coefficients of the compact form.
 * 
 * @todo Seg fault if in last phi-sector becaause of compact eval +1-coefs
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
integer interp3Dcomp_eval_B_SIMD(int i, real B[NSIMD], interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI + phi;}

    int i_r = (r-str->r_min)/str->r_grid;     /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dri = 1.0-dr;
    real dr3 = dr*dr*dr-dr;
    real dri3 = (1.0-dr)*(1.0-dr)*(1.0-dr)-(1.0-dr);
    real rg2 = str->r_grid*str->r_grid;       /**< Square of cell length in r direction */
    int i_phi = (phi-str->phi_min)/str->phi_grid; /**< index for phi variable */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
									   coordinate in
									   current cell */
    real dphii = 1.0-dphi;
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphii3 = (1.0-dphi)*(1.0-dphi)*(1.0-dphi)-(1.0-dphi);
    real phig2 = str->phi_grid*str->phi_grid; /**< Square of cell length in phi dir */
    int i_z = (z-str->z_min)/str->z_grid;     /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */
    real dzi = 1.0-dz;
    real dz3 = dz*dz*dz-dz;
    real dzi3 = (1.0-dz)*(1.0-dz)*(1.0-dz)-(1.0-dz);
    real zg2 = str->z_grid*str->z_grid;       /**< Square of cell length in z direction */
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8; /**< Index jump to cell */
    int r1 = 8;                               /**< Index jump one r forward */
    int phi1 = str->n_z*str->n_r*8;           /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
	phi1 = -(str->n_phi-1)*phi1;          /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r*8;                      /**< Index jump one z forward */
	   
    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {

	B[i] = (
	    dzi*(
		dri*(dphii*str->c[n+0]+dphi*str->c[n+phi1+0])+
		dr*(dphii*str->c[n+r1+0]+dphi*str->c[n+phi1+r1+0]))
	    +dz*(
		dri*(dphii*str->c[n+z1+0]+dphi*str->c[n+phi1+z1+0])+
		dr*(dphii*str->c[n+r1+z1+0]+dphi*str->c[n+phi1+z1+r1+0])))
	    +rg2/6*(
		dzi*(
		    dri3*(dphii*str->c[n+1]+dphi*str->c[n+phi1+1])+
		    dr3*(dphii*str->c[n+r1+1]+dphi*str->c[n+phi1+r1+1]))
		+dz*(
		    dri3*(dphii*str->c[n+z1+1]+dphi*str->c[n+phi1+z1+1])+
		    dr3*(dphii*str->c[n+r1+z1+1]+dphi*str->c[n+phi1+z1+r1+1])))
	    +phig2/6*(
		dzi*(
		    dri*(dphii3*str->c[n+2]+dphi3*str->c[n+phi1+2])+
		    dr*(dphii3*str->c[n+r1+2]+dphi3*str->c[n+phi1+r1+2]))
		+dz*(
		    dri*(dphii3*str->c[n+z1+2]+dphi3*str->c[n+phi1+z1+2])+
		    dr*(dphii3*str->c[n+r1+z1+2]+dphi3*str->c[n+phi1+z1+r1+2])))
	    +zg2/6*(
		dzi3*(
		    dri*(dphii*str->c[n+3]+dphi*str->c[n+phi1+3])+
		    dr*(dphii*str->c[n+r1+3]+dphi*str->c[n+phi1+r1+3]))
		+dz3*(
		    dri*(dphii*str->c[n+z1+3]+dphi*str->c[n+phi1+z1+3])+
		    dr*(dphii*str->c[n+r1+z1+3]+dphi*str->c[n+phi1+z1+r1+3])))
	    +rg2*phig2/36*(
		dzi*(
		    dri3*(dphii3*str->c[n+4]+dphi3*str->c[n+phi1+4])+
		    dr3*(dphii3*str->c[n+r1+4]+dphi3*str->c[n+phi1+r1+4]))
		+dz*(
		    dri3*(dphii3*str->c[n+z1+4]+dphi3*str->c[n+phi1+z1+4])+
		    dr3*(dphii3*str->c[n+r1+z1+4]+dphi3*str->c[n+phi1+z1+r1+4])))
	    +rg2*zg2/36*(
		dzi3*(
		    dri3*(dphii*str->c[n+5]+dphi*str->c[n+phi1+5])+
		    dr3*(dphii*str->c[n+r1+5]+dphi*str->c[n+phi1+r1+5]))
		+dz3*(
		    dri3*(dphii*str->c[n+z1+5]+dphi*str->c[n+phi1+z1+5])+
		    dr3*(dphii*str->c[n+r1+z1+5]+dphi*str->c[n+phi1+z1+r1+5])))
	    +phig2*zg2/36*(
		dzi3*(
		    dri*(dphii3*str->c[n+6]+dphi3*str->c[n+phi1+6])+
		    dr*(dphii3*str->c[n+r1+6]+dphi3*str->c[n+phi1+r1+6]))
		+dz3*(
		    dri*(dphii3*str->c[n+z1+6]+dphi3*str->c[n+phi1+z1+6])+
		    dr*(dphii3*str->c[n+r1+z1+6]+dphi3*str->c[n+phi1+z1+r1+6])))
	    +rg2*phig2*zg2/216*(
		dzi3*(
		    dri3*(dphii3*str->c[n+7]+dphi3*str->c[n+phi1+7])+
		    dr3*(dphii3*str->c[n+r1+7]+dphi3*str->c[n+phi1+r1+7]))
		+dz3*(
		    dri3*(dphii3*str->c[n+z1+7]+dphi3*str->c[n+phi1+z1+7])+
		    dr3*(dphii3*str->c[n+r1+z1+7]+dphi3*str->c[n+phi1+z1+r1+7])));

    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 3D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 * 
 * @todo Seg fault if in last phi-sector becaause of compact eval +1-coefs
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
integer interp3Dcomp_eval_dB(real* B_dB, interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI + phi;}

    int i_r = (r-str->r_min)/str->r_grid;       /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dr3 = dr*dr*dr-dr;
    real dr3dr = 3*dr*dr-1.0;           /**< r-derivative of dr3, not including 1/r_grid */
    real dri = 1.0-dr;
    real dri3 = dri*dri*dri-dri;
    real dri3dr = -3*dri*dri+1.0;      /**< r-derivative of dri3, not including 1/r_grid */
    real rg = str->r_grid;                      /**< Cell length in r direction */
    real rg2 = rg*rg;
    real rgi = 1.0/rg;
    int i_phi = (phi-str->phi_min)/str->phi_grid; /**< index for phi variable */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
									   coordinate in
									   current cell */
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphi3dphi = 3*dphi*dphi-1.0; /**< phi-derivative of dphi3,
					 not including 1/phi_grid */
    real dphii = 1.0-dphi;
    real dphii3 = dphii*dphii*dphii-dphii;
    real dphii3dphi = -3*dphii*dphii+1.0; /**< phi-derivative of dphii3,
					     not including 1/r_grid */
    real phig = str->phi_grid;                  /**< Cell length in phi direction */
    real phig2 = phig*phig;
    real phigi = 1.0/phig;
    int i_z = (z-str->z_min)/str->z_grid;       /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */
    real dz3 = dz*dz*dz-dz;
    real dz3dz = 3*dz*dz-1.0;          /**< z-derivative of dz3, not including 1/z_grid */
    real dzi = 1.0-dz;
    real dzi3 = dzi*dzi*dzi-dzi;
    real dzi3dz = -3*dzi*dzi+1.0;      /**< z-derivative of dzi3, not including 1/z_grid */
    real zg = str->z_grid;                      /**< Cell length in z direction */
    real zg2 = zg*zg;
    real zgi = 1.0/zg;
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8; /**< Index jump to cell */
    int r1 = 8;                                 /**< Index jump one r forward */
    int phi1 = str->n_z*str->n_r*8;             /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
	phi1 = -(str->n_phi-1)*phi1;            /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r*8;                        /**< Index jump one z forward */


    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {

	real c0000 = str->c[n+0];
	real c0001 = str->c[n+1];
	real c0002 = str->c[n+2];
	real c0003 = str->c[n+3];
	real c0004 = str->c[n+4];
	real c0005 = str->c[n+5];
	real c0006 = str->c[n+6];
	real c0007 = str->c[n+7];

	real c0010 = str->c[n+r1+0];
	real c0011 = str->c[n+r1+1];
	real c0012 = str->c[n+r1+2];
	real c0013 = str->c[n+r1+3];
	real c0014 = str->c[n+r1+4];
	real c0015 = str->c[n+r1+5];
	real c0016 = str->c[n+r1+6];
	real c0017 = str->c[n+r1+7];

	real c0100 = str->c[n+phi1+0];
	real c0101 = str->c[n+phi1+1];
	real c0102 = str->c[n+phi1+2];
	real c0103 = str->c[n+phi1+3];
	real c0104 = str->c[n+phi1+4];
	real c0105 = str->c[n+phi1+5];
	real c0106 = str->c[n+phi1+6];
	real c0107 = str->c[n+phi1+7];

	real c1000 = str->c[n+z1+0];
	real c1001 = str->c[n+z1+1];
	real c1002 = str->c[n+z1+2];
	real c1003 = str->c[n+z1+3];
	real c1004 = str->c[n+z1+4];
	real c1005 = str->c[n+z1+5];
	real c1006 = str->c[n+z1+6];
	real c1007 = str->c[n+z1+7];

	real c0110 = str->c[n+phi1+r1+0];
	real c0111 = str->c[n+phi1+r1+1];
	real c0112 = str->c[n+phi1+r1+2];
	real c0113 = str->c[n+phi1+r1+3];
	real c0114 = str->c[n+phi1+r1+4];
	real c0115 = str->c[n+phi1+r1+5];
	real c0116 = str->c[n+phi1+r1+6];
	real c0117 = str->c[n+phi1+r1+7];

	real c1100 = str->c[n+phi1+z1+0];
	real c1101 = str->c[n+phi1+z1+1];
	real c1102 = str->c[n+phi1+z1+2];
	real c1103 = str->c[n+phi1+z1+3];
	real c1104 = str->c[n+phi1+z1+4];
	real c1105 = str->c[n+phi1+z1+5];
	real c1106 = str->c[n+phi1+z1+6];
	real c1107 = str->c[n+phi1+z1+7];

	real c1010 = str->c[n+r1+z1+0];
	real c1011 = str->c[n+r1+z1+1];
	real c1012 = str->c[n+r1+z1+2];
	real c1013 = str->c[n+r1+z1+3];
	real c1014 = str->c[n+r1+z1+4];
	real c1015 = str->c[n+r1+z1+5];
	real c1016 = str->c[n+r1+z1+6];
	real c1017 = str->c[n+r1+z1+7];

	real c1110 = str->c[n+r1+phi1+z1+0];
	real c1111 = str->c[n+r1+phi1+z1+1];
	real c1112 = str->c[n+r1+phi1+z1+2];
	real c1113 = str->c[n+r1+phi1+z1+3];
	real c1114 = str->c[n+r1+phi1+z1+4];
	real c1115 = str->c[n+r1+phi1+z1+5];
	real c1116 = str->c[n+r1+phi1+z1+6];
	real c1117 = str->c[n+r1+phi1+z1+7];


	/* f */
	B_dB[0] = (
	       dzi*(
		   dri*(dphii*c0000+dphi*c0100)+
		   dr*(dphii*c0010+dphi*c0110))
	       +dz*(
		   dri*(dphii*c1000+dphi*c1100)+
		   dr*(dphii*c1010+dphi*c1110)))
	+rg2/6*(
	    dzi*(
		dri3*(dphii*c0001+dphi*c0101)+
		dr3*(dphii*c0011+dphi*c0111))
	    +dz*(
		dri3*(dphii*c1001+dphi*c1101)+
		dr3*(dphii*c1011+dphi*c1111)))
	+phig2/6*(
	    dzi*(
		dri*(dphii3*c0002+dphi3*c0102)+
		dr*(dphii3*c0012+dphi3*c0112))
	    +dz*(
		dri*(dphii3*c1002+dphi3*c1102)+
		dr*(dphii3*c1012+dphi3*c1112)))
	+zg2/6*(
	    dzi3*(
		dri*(dphii*c0003+dphi*c0103)+
		dr*(dphii*c0013+dphi*c0113))
	    +dz3*(
		dri*(dphii*c1003+dphi*c1103)+
		dr*(dphii*c1013+dphi*c1113)))
	+rg2*phig2/36*(
	    dzi*(
		dri3*(dphii3*c0004+dphi3*c0104)+
		dr3*(dphii3*c0014+dphi3*c0114))
	    +dz*(
		dri3*(dphii3*c1004+dphi3*c1104)+
		dr3*(dphii3*c1014+dphi3*c1114)))
	+rg2*zg2/36*(
	    dzi3*(
		dri3*(dphii*c0005+dphi*c0105)+
		dr3*(dphii*c0015+dphi*c0115))
	    +dz3*(
		dri3*(dphii*c1005+dphi*c1105)+
		dr3*(dphii*c1015+dphi*c1115)))
	+phig2*zg2/36*(
	    dzi3*(
		dri*(dphii3*c0006+dphi3*c0106)+
		dr*(dphii3*c0016+dphi3*c0116))
	    +dz3*(
		dri*(dphii3*c1006+dphi3*c1106)+
		dr*(dphii3*c1016+dphi3*c1116)))
	+rg2*phig2*zg2/216*(
	    dzi3*(
		dri3*(dphii3*c0007+dphi3*c0107)+
		dr3*(dphii3*c0017+dphi3*c0117))
	    +dz3*(
		dri3*(dphii3*c1007+dphi3*c1107)+
		dr3*(dphii3*c1017+dphi3*c1117)));
    
    /* df/dr */
    B_dB[1] = rgi*(
	dzi*(
	    -(dphii*c0000+dphi*c0100)
	    +(dphii*c0010+dphi*c0110))
	+dz*(
	    -(dphii*c1000+dphi*c1100)
	    +(dphii*c1010+dphi*c1110)))
	+rg/6*(
	    dzi*(
		dri3dr*(dphii*c0001+dphi*c0101)+
		dr3dr*(dphii*c0011+dphi*c0111))
	    +dz*(
		dri3dr*(dphii*c1001  +dphi*c1101)+
		dr3dr*(dphii*c1011+dphi*c1111)))
	+rgi*phig2/6*(
	    dzi*(
		-(dphii3*c0002+dphi3*c0102)
		+(dphii3*c0012+dphi3*c0112))
	    +dz*(
		-(dphii3*c1002+dphi3*c1102)
		+(dphii3*c1012+dphi3*c1112)))
	+rgi*zg2/6*(
	    dzi3*(
		-(dphii*c0003+dphi*c0103)
		+(dphii*c0013+dphi*c0113))
	    +dz3*(
		-(dphii*c1003+dphi*c1103)
		+(dphii*c1013+dphi*c1113)))
	+rg*phig2/36*(
	    dzi*(
		dri3dr*(dphii3*c0004+dphi3*c0104)+
		dr3dr*(dphii3*c0014+dphi3*c0114))
	    +dz*(
		dri3dr*(dphii3*c1004+dphi3*c1104)+
		dr3dr*(dphii3*c1014+dphi3*c1114)))
	+rg*zg2/36*(
	    dzi3*(
		dri3dr*(dphii*c0005+dphi*c0105)+
		dr3dr*(dphii*c0015+dphi*c0115))
	    +dz3*(
		dri3dr*(dphii*c1005+dphi*c1105)+
		dr3dr*(dphii*c1015+dphi*c1115)))
	+rgi*phig2*zg2/36*(
	    dzi3*(
		-(dphii3*c0006+dphi3*c0106)
		+(dphii3*c0016+dphi3*c0116))
	    +dz3*(
		-(dphii3*c1006+dphi3*c1106)
		+(dphii3*c1016+dphi3*c1116)))
	+rg*phig2*zg2/216*(
	    dzi3*(
		dri3dr*(dphii3*c0007+dphi3*c0107)+
		dr3dr*(dphii3*c0017+dphi3*c0117))
	    +dz3*(
		dri3dr*(dphii3*c1007+dphi3*c1107)+
		dr3dr*(dphii3*c1017+dphi3*c1117)));
    
    /* df/dphi */
    B_dB[2] = phigi*(
	dzi*(
	    dri*(-c0000+c0100)+
	    dr*(-c0010+c0110))
	+dz*(
	    dri*(-c1000+c1100)+
	    dr*(-c1010+c1110)))
	+phigi*rg2/6*(
	    dzi*(
		dri3*(-c0001+c0101)+
		dr3*(-c0011+c0111))
	    +dz*(
		dri3*(-c1001+c1101)+
		dr3*(-c1011+c1111)))
	+phig/6*(
	    dzi*(
		dri*(dphii3dphi*c0002+dphi3dphi*c0102)+
		dr*(dphii3dphi*c0012+dphi3dphi*c0112))
	    +dz*(
		dri*(dphii3dphi*c1002+dphi3dphi*c1102)+
		dr*(dphii3dphi*c1012+dphi3dphi*c1112)))
	+phigi*zg2/6*(
	    dzi3*(
		dri*(-c0003+c0103)+
		dr*(-c0013+c0113))
	    +dz3*(
		dri*(-c1003+c1103)+
		dr*(-c1013+c1113)))
	+rg2*phig/36*(
	    dzi*(
		dri3*(dphii3dphi*c0004+dphi3dphi*c0104)+
		dr3*(dphii3dphi*c0014+dphi3dphi*c0114))
	    +dz*(
		dri3*(dphii3dphi*c1004+dphi3dphi*c1104)+
		dr3*(dphii3dphi*c1014+dphi3dphi*c1114)))
	+phigi*rg2*zg2/36*(
	    dzi3*(
		dri3*(-c0005+c0105)+
		dr3*(-c0015+c0115))
	    +dz3*(
		dri3*(-c1005+c1105)+
		dr3*(-c1015+c1115)))
	+phig*zg2/36*(
	    dzi3*(
		dri*(dphii3dphi*c0006+dphi3dphi*c0106)+
		dr*(dphii3dphi*c0016+dphi3dphi*c0116))
	    +dz3*(
		dri*(dphii3dphi*c1006+dphi3dphi*c1106)+
		dr*(dphii3dphi*c1016+dphi3dphi*c1116)))
	+rg2*phig*zg2/216*(
	    dzi3*(
		dri3*(dphii3dphi*c0007+dphi3dphi*c0107)+
		dr3*(dphii3dphi*c0017+dphi3dphi*c0117))
	    +dz3*(
		dri3*(dphii3dphi*c1007+dphi3dphi*c1107)+
		dr3*(dphii3dphi*c1017+dphi3dphi*c1117)));
    
    /* df/dz */
    B_dB[3] = zgi*(
	-(
	    dri*(dphii*c0000+dphi*c0100)+
	    dr*(dphii*c0010+dphi*c0110))
	+(
	    dri*(dphii*c1000+dphi*c1100)+
	    dr*(dphii*c1010+dphi*c1110)))
	+rg2*zgi/6*(
	    -(
		dri3*(dphii*c0001+dphi*c0101)+
		dr3*(dphii*c0011+dphi*c0111))
	    +(
		dri3*(dphii*c1001+dphi*c1101)+
		dr3*(dphii*c1011+dphi*c1111)))
	+phig2*zgi/6*(
	    -(
		dri*(dphii3*c0002+dphi3*c0102)+
		dr*(dphii3*c0012+dphi3*c0112))
	    +(
		dri*(dphii3*c1002+dphi3*c1102)+
		dr*(dphii3*c1012+dphi3*c1112)))
	+zg/6*(
	    dzi3dz*(
		dri*(dphii*c0003+dphi*c0103)+
		dr*(dphii*c0013+dphi*c0113))
	    +dz3dz*(
		dri*(dphii*c1003+dphi*c1103)+
		dr*(dphii*c1013+dphi*c1113)))
	+rg2*phig2*zgi/36*(
	    -(
		dri3*(dphii3*c0004+dphi3*c0104)+
		dr3*(dphii3*c0014+dphi3*c0114))
	    +(
		dri3*(dphii3*c1004+dphi3*c1104)+
		dr3*(dphii3*c1014+dphi3*c1114)))
	+rg2*zg/36*(
	    dzi3dz*(
		dri3*(dphii*c0005+dphi*c0105)+
		dr3*(dphii*c0015+dphi*c0115))
	    +dz3dz*(
		dri3*(dphii*c1005+dphi*c1105)+
		dr3*(dphii*c1015+dphi*c1115)))
	+phig2*zg/36*(
	    dzi3dz*(
		dri*(dphii3*c0006+dphi3*c0106)+
		dr*(dphii3*c0016+dphi3*c0116))
	    +dz3dz*(
		dri*(dphii3*c1006+dphi3*c1106)+
		dr*(dphii3*c1016+dphi3*c1116)))
	+rg2*phig2*zg/216*(
	    dzi3dz*(
		dri3*(dphii3*c0007+dphi3*c0107)+
		dr3*(dphii3*c0017+dphi3*c0117))
	    +dz3dz*(
		dri3*(dphii3*c1007+dphi3*c1107)+
		dr3*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/dr2 */
    B_dB[4] = (
	dzi*(
	    dri*(dphii*c0001+dphi*c0101)+
	    dr*(dphii*c0011+dphi*c0111))
	+dz*(
	    dri*(dphii*c1001+dphi*c1101)+
	    dr*(dphii*c1011+dphi*c1111)))
	+phig2/6*(
	    dzi*(
		dri*(dphii3*c0004+dphi3*c0104)+
		dr*(dphii3*c0014+dphi3*c0114))
	    +dz*(
		dri*(dphii3*c1004+dphi3*c1104)+
		dr*(dphii3*c1014+dphi3*c1114)))
	+zg2/6*(
	    dzi3*(
		dri*(dphii*c0005+dphi*c0105)+
		dr*(dphii*c0015+dphi*c0115))
	    +dz3*(
		dri*(dphii*c1005+dphi*c1105)+
		dr*(dphii*c1015+dphi*c1115)))
	+phig2*zg2/36*(
	    dzi3*(
		dri*(dphii3*c0007+dphi3*c0107)+
		dr*(dphii3*c0017+dphi3*c0117))
	    +dz3*(
		dri*(dphii3*c1007+dphi3*c1107)+
		dr*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/dphi2 */
    B_dB[5] = (
	dzi*(
	    dri*(dphii*c0002+dphi*c0102)+
	    dr*(dphii*c0012+dphi*c0112))
	+dz*(
	    dri*(dphii*c1002+dphi*c1102)+
	    dr*(dphii*c1012+dphi*c1112)))
	+rg2/6*(
	    dzi*(
		dri3*(dphii*c0004+dphi*c0104)+
		dr3*(dphii*c0014+dphi*c0114))
	    +dz*(
		dri3*(dphii*c1004+dphi*c1104)+
		dr3*(dphii*c1014+dphi*c1114)))
	+zg2/6*(
	    dzi3*(
		dri*(dphii*c0006+dphi*c0106)+
		dr*(dphii*c0016+dphi*c0116))
	    +dz3*(
		dri*(dphii*c1006+dphi*c1106)+
		dr*(dphii*c1016+dphi*c1116)))
	+rg2*zg2/36*(
	    dzi3*(
		dri3*(dphii*c0007+dphi*c0107)+
		dr3*(dphii*c0017+dphi*c0117))
	    +dz3*(
		dri3*(dphii*c1007+dphi*c1107)+
		dr3*(dphii*c1017+dphi*c1117)));
    
    /* d2f/dz2 */
    B_dB[6] = (
	dzi*(
	    dri*(dphii*c0003+dphi*c0103)+
	    dr*(dphii*c0013+dphi*c0113))
	+dz*(
	    dri*(dphii*c1003+dphi*c1103)+
	    dr*(dphii*c1013+dphi*c1113)))
	+rg2/6*(
	    dzi*(
		dri3*(dphii*c0005+dphi*c0105)+
		dr3*(dphii*c0015+dphi*c0115))
	    +dz*(
		dri3*(dphii*c1005+dphi*c1105)+
		dr3*(dphii*c1015+dphi*c1115)))
	+phig2/6*(
	    dzi*(
		dri*(dphii3*c0006+dphi3*c0106)+
		dr*(dphii3*c0016+dphi3*c0116))
	    +dz*(
		dri*(dphii3*c1006+dphi3*c1106)+
		dr*(dphii3*c1016+dphi3*c1116)))
	+rg2*phig2/36*(
	    dzi*(
		dri3*(dphii3*c0007+dphi3*c0107)+
		dr3*(dphii3*c0017+dphi3*c0117))
	    +dz*(
		dri3*(dphii3*c1007+dphi3*c1107)+
		dr3*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/drdphi */
    B_dB[7] = rgi*phigi*(
	dzi*(
	    (c0000  -c0100)-
	    (c0010-c0110))
	+dz*(
	    (c1000  -c1100)-
	    (c1010-c1110)))
	+phigi*rg/6*(
	    dzi*(
		dri3dr*(-c0001+c0101)+
		dr3dr*(-c0011+c0111))
	    +dz*(
		dri3dr*(-c1001+c1101)+
		dr3dr*(-c1011+c1111)))
	+rgi*phig/6*(
	    dzi*(
		-(dphii3dphi*c0002+dphi3dphi*c0102)
		+(dphii3dphi*c0012+dphi3dphi*c0112))
	    +dz*(
		-(dphii3dphi*c1002+dphi3dphi*c1102)
		+(dphii3dphi*c1012+dphi3dphi*c1112)))
	+rgi*phigi*zg2/6*(
	    dzi3*(
		(c0003  -c0103)-
		(c0013-c0113))
	    +dz3*(
		(c1003  -c1103)-
		(c1013-c1113)))
	+rg*phig/36*(
	    dzi*(
		dri3dr*(dphii3dphi*c0004+dphi3dphi*c0104)+
		dr3dr*(dphii3dphi*c0014+dphi3dphi*c0114))
	    +dz*(
		dri3dr*(dphii3dphi*c1004+dphi3dphi*c1104)+
		dr3dr*(dphii3dphi*c1014+dphi3dphi*c1114)))
	+phigi*rg*zg2/36*(
	    dzi3*(
		dri3dr*(-c0005+c0105)+
		dr3dr*(-c0015+c0115))
	    +dz3*(
		dri3dr*(-c1005+c1105)+
		dr3dr*(-c1015+c1115)))
	+rgi*phig*zg2/36*(
	    dzi3*(
		-(dphii3dphi*c0006+dphi3dphi*c0106)
		+(dphii3dphi*c0016+dphi3dphi*c0116))
	    +dz3*(
		-(dphii3dphi*c1006+dphi3dphi*c1106)
		+(dphii3dphi*c1016+dphi3dphi*c1116)))
	+rg*phig*zg2/216*(
	    dzi3*(
		dri3dr*(dphii3dphi*c0007+dphi3dphi*c0107)+
		dr3dr*(dphii3dphi*c0017+dphi3dphi*c0117))
	    +dz3*(
		dri3dr*(dphii3dphi*c1007+dphi3dphi*c1107)+
		dr3dr*(dphii3dphi*c1017+dphi3dphi*c1117)));
    
    /* d2f/drdz */
    B_dB[8] = rgi*zgi*(
	(
	    (dphii*c0000+dphi*c0100) -
	    (dphii*c0010+dphi*c0110))
	-(
	    (dphii*c1000+dphi*c1100) -
	    (dphii*c1010+dphi*c1110)))
	+rg*zgi/6*(
	    -(
		dri3dr*(dphii*c0001+dphi*c0101)+
		dr3dr*(dphii*c0011+dphi*c0111))
	    +(
		dri3dr*(dphii*c1001+dphi*c1101)+
		dr3dr*(dphii*c1011+dphi*c1111)))
	+rgi*phig2*zgi/6*(
	    (
		(dphii3*c0002+dphi3*c0102) -
		(dphii3*c0012+dphi3*c0112))
	    -(
		(dphii3*c1002+dphi3*c1102) -
		(dphii3*c1012+dphi3*c1112)))
	+rgi*zg/6*(
	    dzi3dz*(
		-(dphii*c0003+dphi*c0103)
		+(dphii*c0013+dphi*c0113))
	    +dz3dz*(
		-(dphii*c1003+dphi*c1103)
		+(dphii*c1013+dphi*c1113)))
	+rg*phig2*zgi/36*(
	    -(
		dri3dr*(dphii3*c0004+dphi3*c0104)+
		dr3dr*(dphii3*c0014+dphi3*c0114))
	    +(
		dri3dr*(dphii3*c1004+dphi3*c1104)+
		dr3dr*(dphii3*c1014+dphi3*c1114)))
	+rg*zg/36*(
	    dzi3dz*(
		dri3dr*(dphii*c0005+dphi*c0105)+
		dr3dr*(dphii*c0015+dphi*c0115))
	    +dz3dz*(
		dri3dr*(dphii*c1005+dphi*c1105)+
		dr3dr*(dphii*c1015+dphi*c1115)))
	+rgi*phig2*zg/36*(
	    dzi3dz*(
		-(dphii3*c0006+dphi3*c0106)
		+(dphii3*c0016+dphi3*c0116))
	    +dz3dz*(
		-(dphii3*c1006+dphi3*c1106)
		+(dphii3*c1016+dphi3*c1116)))
	+rg*phig2*zg/216*(
	    dzi3dz*(
		dri3dr*(dphii3*c0007+dphi3*c0107)+
		dr3dr*(dphii3*c0017+dphi3*c0117))
	    +dz3dz*(
		dri3dr*(dphii3*c1007+dphi3*c1107)+
		dr3dr*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/dphidz */
    B_dB[9] = phigi*zgi*(
	(
	    dri*(c0000  -c0100)+
	    dr*(c0010-c0110))
	-(
	    dri*(c1000  -c1100)+
	    dr*(c1010-c1110)))
	+phigi*rg2*zgi/6*(
	    (
		dri3*(c0001  -c0101)+
		dr3*(c0011-c0111))
	    -(
		dri3*(c1001  -c1101)+
		dr3*(c1011-c1111)))
	+phig*zgi/6*(
	    -(
		dri*(dphii3dphi*c0002+dphi3dphi*c0102)+
		dr*(dphii3dphi*c0012+dphi3dphi*c0112))
	    +(
		dri*(dphii3dphi*c1002+dphi3dphi*c1102)+
		dr*(dphii3dphi*c1012+dphi3dphi*c1112)))
	+phigi*zg/6*(
	    dzi3dz*(
		dri*(-c0003+c0103)+
		dr*(-c0013+c0113))
	    +dz3dz*(
		dri*(-c1003+c1103)+
		dr*(-c1013+c1113)))
	+rg2*phig*zgi/36*(
	    -(
		dri3*(dphii3dphi*c0004+dphi3dphi*c0104)+
		dr3*(dphii3dphi*c0014+dphi3dphi*c0114))
	    +(
		dri3*(dphii3dphi*c1004+dphi3dphi*c1104)+
		dr3*(dphii3dphi*c1014+dphi3dphi*c1114)))
	+phigi*rg2*zg/36*(
	    dzi3dz*(
		dri3*(-c0005+c0105)+
		dr3*(-c0015+c0115))
	    +dz3dz*(
		dri3*(-c1005+c1105)+
		dr3*(-c1015+c1115)))
	+phig*zg/36*(
	    dzi3dz*(
		dri*(dphii3dphi*c0006+dphi3dphi*c0106)+
		dr*(dphii3dphi*c0016+dphi3dphi*c0116))
	    +dz3dz*(
		dri*(dphii3dphi*c1006+dphi3dphi*c1106)+
		dr*(dphii3dphi*c1016+dphi3dphi*c1116)))
	+rg2*phig*zg/216*(
	    dzi3dz*(
		dri3*(dphii3dphi*c0007+dphi3dphi*c0107)+
		dr3*(dphii3dphi*c0017+dphi3dphi*c0117))
	    +dz3dz*(
		dri3*(dphii3dphi*c1007+dphi3dphi*c1107)+
		dr3*(dphii3dphi*c1017+dphi3dphi*c1117)));

    }

    return err;
}

integer interp3Dcomp_eval_dB_SIMD(int i, real B_dB[10][NSIMD], interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI + phi;}

    int i_r = (r-str->r_min)/str->r_grid;       /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dr3 = dr*dr*dr-dr;
    real dr3dr = 3*dr*dr-1.0;           /**< r-derivative of dr3, not including 1/r_grid */
    real dri = 1.0-dr;
    real dri3 = dri*dri*dri-dri;
    real dri3dr = -3*dri*dri+1.0;      /**< r-derivative of dri3, not including 1/r_grid */
    real rg = str->r_grid;                      /**< Cell length in r direction */
    real rg2 = rg*rg;
    real rgi = 1.0/rg;
    int i_phi = (phi-str->phi_min)/str->phi_grid; /**< index for phi variable */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
									   coordinate in
									   current cell */
    real dphi3 = dphi*dphi*dphi-dphi;
    real dphi3dphi = 3*dphi*dphi-1.0; /**< phi-derivative of dphi3,
					 not including 1/phi_grid */
    real dphii = 1.0-dphi;
    real dphii3 = dphii*dphii*dphii-dphii;
    real dphii3dphi = -3*dphii*dphii+1.0; /**< phi-derivative of dphii3,
					     not including 1/r_grid */
    real phig = str->phi_grid;                  /**< Cell length in phi direction */
    real phig2 = phig*phig;
    real phigi = 1.0/phig;
    int i_z = (z-str->z_min)/str->z_grid;       /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */
    real dz3 = dz*dz*dz-dz;
    real dz3dz = 3*dz*dz-1.0;          /**< z-derivative of dz3, not including 1/z_grid */
    real dzi = 1.0-dz;
    real dzi3 = dzi*dzi*dzi-dzi;
    real dzi3dz = -3*dzi*dzi+1.0;      /**< z-derivative of dzi3, not including 1/z_grid */
    real zg = str->z_grid;                      /**< Cell length in z direction */
    real zg2 = zg*zg;
    real zgi = 1.0/zg;
    int n = i_phi*str->n_z*str->n_r*8+i_z*str->n_r*8+i_r*8; /**< Index jump to cell */
    int r1 = 8;                                 /**< Index jump one r forward */
    int phi1 = str->n_z*str->n_r*8;             /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
	phi1 = -(str->n_phi-1)*phi1;            /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r*8;                        /**< Index jump one z forward */


    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {

	real c0000 = str->c[n+0];
	real c0001 = str->c[n+1];
	real c0002 = str->c[n+2];
	real c0003 = str->c[n+3];
	real c0004 = str->c[n+4];
	real c0005 = str->c[n+5];
	real c0006 = str->c[n+6];
	real c0007 = str->c[n+7];

	real c0010 = str->c[n+r1+0];
	real c0011 = str->c[n+r1+1];
	real c0012 = str->c[n+r1+2];
	real c0013 = str->c[n+r1+3];
	real c0014 = str->c[n+r1+4];
	real c0015 = str->c[n+r1+5];
	real c0016 = str->c[n+r1+6];
	real c0017 = str->c[n+r1+7];

	real c0100 = str->c[n+phi1+0];
	real c0101 = str->c[n+phi1+1];
	real c0102 = str->c[n+phi1+2];
	real c0103 = str->c[n+phi1+3];
	real c0104 = str->c[n+phi1+4];
	real c0105 = str->c[n+phi1+5];
	real c0106 = str->c[n+phi1+6];
	real c0107 = str->c[n+phi1+7];

	real c1000 = str->c[n+z1+0];
	real c1001 = str->c[n+z1+1];
	real c1002 = str->c[n+z1+2];
	real c1003 = str->c[n+z1+3];
	real c1004 = str->c[n+z1+4];
	real c1005 = str->c[n+z1+5];
	real c1006 = str->c[n+z1+6];
	real c1007 = str->c[n+z1+7];

	real c0110 = str->c[n+phi1+r1+0];
	real c0111 = str->c[n+phi1+r1+1];
	real c0112 = str->c[n+phi1+r1+2];
	real c0113 = str->c[n+phi1+r1+3];
	real c0114 = str->c[n+phi1+r1+4];
	real c0115 = str->c[n+phi1+r1+5];
	real c0116 = str->c[n+phi1+r1+6];
	real c0117 = str->c[n+phi1+r1+7];

	real c1100 = str->c[n+phi1+z1+0];
	real c1101 = str->c[n+phi1+z1+1];
	real c1102 = str->c[n+phi1+z1+2];
	real c1103 = str->c[n+phi1+z1+3];
	real c1104 = str->c[n+phi1+z1+4];
	real c1105 = str->c[n+phi1+z1+5];
	real c1106 = str->c[n+phi1+z1+6];
	real c1107 = str->c[n+phi1+z1+7];

	real c1010 = str->c[n+r1+z1+0];
	real c1011 = str->c[n+r1+z1+1];
	real c1012 = str->c[n+r1+z1+2];
	real c1013 = str->c[n+r1+z1+3];
	real c1014 = str->c[n+r1+z1+4];
	real c1015 = str->c[n+r1+z1+5];
	real c1016 = str->c[n+r1+z1+6];
	real c1017 = str->c[n+r1+z1+7];

	real c1110 = str->c[n+r1+phi1+z1+0];
	real c1111 = str->c[n+r1+phi1+z1+1];
	real c1112 = str->c[n+r1+phi1+z1+2];
	real c1113 = str->c[n+r1+phi1+z1+3];
	real c1114 = str->c[n+r1+phi1+z1+4];
	real c1115 = str->c[n+r1+phi1+z1+5];
	real c1116 = str->c[n+r1+phi1+z1+6];
	real c1117 = str->c[n+r1+phi1+z1+7];


	/* f */
	B_dB[0][i] = (
	       dzi*(
		   dri*(dphii*c0000+dphi*c0100)+
		   dr*(dphii*c0010+dphi*c0110))
	       +dz*(
		   dri*(dphii*c1000+dphi*c1100)+
		   dr*(dphii*c1010+dphi*c1110)))
	+rg2/6*(
	    dzi*(
		dri3*(dphii*c0001+dphi*c0101)+
		dr3*(dphii*c0011+dphi*c0111))
	    +dz*(
		dri3*(dphii*c1001+dphi*c1101)+
		dr3*(dphii*c1011+dphi*c1111)))
	+phig2/6*(
	    dzi*(
		dri*(dphii3*c0002+dphi3*c0102)+
		dr*(dphii3*c0012+dphi3*c0112))
	    +dz*(
		dri*(dphii3*c1002+dphi3*c1102)+
		dr*(dphii3*c1012+dphi3*c1112)))
	+zg2/6*(
	    dzi3*(
		dri*(dphii*c0003+dphi*c0103)+
		dr*(dphii*c0013+dphi*c0113))
	    +dz3*(
		dri*(dphii*c1003+dphi*c1103)+
		dr*(dphii*c1013+dphi*c1113)))
	+rg2*phig2/36*(
	    dzi*(
		dri3*(dphii3*c0004+dphi3*c0104)+
		dr3*(dphii3*c0014+dphi3*c0114))
	    +dz*(
		dri3*(dphii3*c1004+dphi3*c1104)+
		dr3*(dphii3*c1014+dphi3*c1114)))
	+rg2*zg2/36*(
	    dzi3*(
		dri3*(dphii*c0005+dphi*c0105)+
		dr3*(dphii*c0015+dphi*c0115))
	    +dz3*(
		dri3*(dphii*c1005+dphi*c1105)+
		dr3*(dphii*c1015+dphi*c1115)))
	+phig2*zg2/36*(
	    dzi3*(
		dri*(dphii3*c0006+dphi3*c0106)+
		dr*(dphii3*c0016+dphi3*c0116))
	    +dz3*(
		dri*(dphii3*c1006+dphi3*c1106)+
		dr*(dphii3*c1016+dphi3*c1116)))
	+rg2*phig2*zg2/216*(
	    dzi3*(
		dri3*(dphii3*c0007+dphi3*c0107)+
		dr3*(dphii3*c0017+dphi3*c0117))
	    +dz3*(
		dri3*(dphii3*c1007+dphi3*c1107)+
		dr3*(dphii3*c1017+dphi3*c1117)));
    
    /* df/dr */
    B_dB[1][i] = rgi*(
	dzi*(
	    -(dphii*c0000+dphi*c0100)
	    +(dphii*c0010+dphi*c0110))
	+dz*(
	    -(dphii*c1000+dphi*c1100)
	    +(dphii*c1010+dphi*c1110)))
	+rg/6*(
	    dzi*(
		dri3dr*(dphii*c0001+dphi*c0101)+
		dr3dr*(dphii*c0011+dphi*c0111))
	    +dz*(
		dri3dr*(dphii*c1001  +dphi*c1101)+
		dr3dr*(dphii*c1011+dphi*c1111)))
	+rgi*phig2/6*(
	    dzi*(
		-(dphii3*c0002+dphi3*c0102)
		+(dphii3*c0012+dphi3*c0112))
	    +dz*(
		-(dphii3*c1002+dphi3*c1102)
		+(dphii3*c1012+dphi3*c1112)))
	+rgi*zg2/6*(
	    dzi3*(
		-(dphii*c0003+dphi*c0103)
		+(dphii*c0013+dphi*c0113))
	    +dz3*(
		-(dphii*c1003+dphi*c1103)
		+(dphii*c1013+dphi*c1113)))
	+rg*phig2/36*(
	    dzi*(
		dri3dr*(dphii3*c0004+dphi3*c0104)+
		dr3dr*(dphii3*c0014+dphi3*c0114))
	    +dz*(
		dri3dr*(dphii3*c1004+dphi3*c1104)+
		dr3dr*(dphii3*c1014+dphi3*c1114)))
	+rg*zg2/36*(
	    dzi3*(
		dri3dr*(dphii*c0005+dphi*c0105)+
		dr3dr*(dphii*c0015+dphi*c0115))
	    +dz3*(
		dri3dr*(dphii*c1005+dphi*c1105)+
		dr3dr*(dphii*c1015+dphi*c1115)))
	+rgi*phig2*zg2/36*(
	    dzi3*(
		-(dphii3*c0006+dphi3*c0106)
		+(dphii3*c0016+dphi3*c0116))
	    +dz3*(
		-(dphii3*c1006+dphi3*c1106)
		+(dphii3*c1016+dphi3*c1116)))
	+rg*phig2*zg2/216*(
	    dzi3*(
		dri3dr*(dphii3*c0007+dphi3*c0107)+
		dr3dr*(dphii3*c0017+dphi3*c0117))
	    +dz3*(
		dri3dr*(dphii3*c1007+dphi3*c1107)+
		dr3dr*(dphii3*c1017+dphi3*c1117)));
    
    /* df/dphi */
    B_dB[2][i] = phigi*(
	dzi*(
	    dri*(-c0000+c0100)+
	    dr*(-c0010+c0110))
	+dz*(
	    dri*(-c1000+c1100)+
	    dr*(-c1010+c1110)))
	+phigi*rg2/6*(
	    dzi*(
		dri3*(-c0001+c0101)+
		dr3*(-c0011+c0111))
	    +dz*(
		dri3*(-c1001+c1101)+
		dr3*(-c1011+c1111)))
	+phig/6*(
	    dzi*(
		dri*(dphii3dphi*c0002+dphi3dphi*c0102)+
		dr*(dphii3dphi*c0012+dphi3dphi*c0112))
	    +dz*(
		dri*(dphii3dphi*c1002+dphi3dphi*c1102)+
		dr*(dphii3dphi*c1012+dphi3dphi*c1112)))
	+phigi*zg2/6*(
	    dzi3*(
		dri*(-c0003+c0103)+
		dr*(-c0013+c0113))
	    +dz3*(
		dri*(-c1003+c1103)+
		dr*(-c1013+c1113)))
	+rg2*phig/36*(
	    dzi*(
		dri3*(dphii3dphi*c0004+dphi3dphi*c0104)+
		dr3*(dphii3dphi*c0014+dphi3dphi*c0114))
	    +dz*(
		dri3*(dphii3dphi*c1004+dphi3dphi*c1104)+
		dr3*(dphii3dphi*c1014+dphi3dphi*c1114)))
	+phigi*rg2*zg2/36*(
	    dzi3*(
		dri3*(-c0005+c0105)+
		dr3*(-c0015+c0115))
	    +dz3*(
		dri3*(-c1005+c1105)+
		dr3*(-c1015+c1115)))
	+phig*zg2/36*(
	    dzi3*(
		dri*(dphii3dphi*c0006+dphi3dphi*c0106)+
		dr*(dphii3dphi*c0016+dphi3dphi*c0116))
	    +dz3*(
		dri*(dphii3dphi*c1006+dphi3dphi*c1106)+
		dr*(dphii3dphi*c1016+dphi3dphi*c1116)))
	+rg2*phig*zg2/216*(
	    dzi3*(
		dri3*(dphii3dphi*c0007+dphi3dphi*c0107)+
		dr3*(dphii3dphi*c0017+dphi3dphi*c0117))
	    +dz3*(
		dri3*(dphii3dphi*c1007+dphi3dphi*c1107)+
		dr3*(dphii3dphi*c1017+dphi3dphi*c1117)));
    
    /* df/dz */
    B_dB[3][i] = zgi*(
	-(
	    dri*(dphii*c0000+dphi*c0100)+
	    dr*(dphii*c0010+dphi*c0110))
	+(
	    dri*(dphii*c1000+dphi*c1100)+
	    dr*(dphii*c1010+dphi*c1110)))
	+rg2*zgi/6*(
	    -(
		dri3*(dphii*c0001+dphi*c0101)+
		dr3*(dphii*c0011+dphi*c0111))
	    +(
		dri3*(dphii*c1001+dphi*c1101)+
		dr3*(dphii*c1011+dphi*c1111)))
	+phig2*zgi/6*(
	    -(
		dri*(dphii3*c0002+dphi3*c0102)+
		dr*(dphii3*c0012+dphi3*c0112))
	    +(
		dri*(dphii3*c1002+dphi3*c1102)+
		dr*(dphii3*c1012+dphi3*c1112)))
	+zg/6*(
	    dzi3dz*(
		dri*(dphii*c0003+dphi*c0103)+
		dr*(dphii*c0013+dphi*c0113))
	    +dz3dz*(
		dri*(dphii*c1003+dphi*c1103)+
		dr*(dphii*c1013+dphi*c1113)))
	+rg2*phig2*zgi/36*(
	    -(
		dri3*(dphii3*c0004+dphi3*c0104)+
		dr3*(dphii3*c0014+dphi3*c0114))
	    +(
		dri3*(dphii3*c1004+dphi3*c1104)+
		dr3*(dphii3*c1014+dphi3*c1114)))
	+rg2*zg/36*(
	    dzi3dz*(
		dri3*(dphii*c0005+dphi*c0105)+
		dr3*(dphii*c0015+dphi*c0115))
	    +dz3dz*(
		dri3*(dphii*c1005+dphi*c1105)+
		dr3*(dphii*c1015+dphi*c1115)))
	+phig2*zg/36*(
	    dzi3dz*(
		dri*(dphii3*c0006+dphi3*c0106)+
		dr*(dphii3*c0016+dphi3*c0116))
	    +dz3dz*(
		dri*(dphii3*c1006+dphi3*c1106)+
		dr*(dphii3*c1016+dphi3*c1116)))
	+rg2*phig2*zg/216*(
	    dzi3dz*(
		dri3*(dphii3*c0007+dphi3*c0107)+
		dr3*(dphii3*c0017+dphi3*c0117))
	    +dz3dz*(
		dri3*(dphii3*c1007+dphi3*c1107)+
		dr3*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/dr2 */
    B_dB[4][i] = (
	dzi*(
	    dri*(dphii*c0001+dphi*c0101)+
	    dr*(dphii*c0011+dphi*c0111))
	+dz*(
	    dri*(dphii*c1001+dphi*c1101)+
	    dr*(dphii*c1011+dphi*c1111)))
	+phig2/6*(
	    dzi*(
		dri*(dphii3*c0004+dphi3*c0104)+
		dr*(dphii3*c0014+dphi3*c0114))
	    +dz*(
		dri*(dphii3*c1004+dphi3*c1104)+
		dr*(dphii3*c1014+dphi3*c1114)))
	+zg2/6*(
	    dzi3*(
		dri*(dphii*c0005+dphi*c0105)+
		dr*(dphii*c0015+dphi*c0115))
	    +dz3*(
		dri*(dphii*c1005+dphi*c1105)+
		dr*(dphii*c1015+dphi*c1115)))
	+phig2*zg2/36*(
	    dzi3*(
		dri*(dphii3*c0007+dphi3*c0107)+
		dr*(dphii3*c0017+dphi3*c0117))
	    +dz3*(
		dri*(dphii3*c1007+dphi3*c1107)+
		dr*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/dphi2 */
    B_dB[5][i] = (
	dzi*(
	    dri*(dphii*c0002+dphi*c0102)+
	    dr*(dphii*c0012+dphi*c0112))
	+dz*(
	    dri*(dphii*c1002+dphi*c1102)+
	    dr*(dphii*c1012+dphi*c1112)))
	+rg2/6*(
	    dzi*(
		dri3*(dphii*c0004+dphi*c0104)+
		dr3*(dphii*c0014+dphi*c0114))
	    +dz*(
		dri3*(dphii*c1004+dphi*c1104)+
		dr3*(dphii*c1014+dphi*c1114)))
	+zg2/6*(
	    dzi3*(
		dri*(dphii*c0006+dphi*c0106)+
		dr*(dphii*c0016+dphi*c0116))
	    +dz3*(
		dri*(dphii*c1006+dphi*c1106)+
		dr*(dphii*c1016+dphi*c1116)))
	+rg2*zg2/36*(
	    dzi3*(
		dri3*(dphii*c0007+dphi*c0107)+
		dr3*(dphii*c0017+dphi*c0117))
	    +dz3*(
		dri3*(dphii*c1007+dphi*c1107)+
		dr3*(dphii*c1017+dphi*c1117)));
    
    /* d2f/dz2 */
    B_dB[6][i] = (
	dzi*(
	    dri*(dphii*c0003+dphi*c0103)+
	    dr*(dphii*c0013+dphi*c0113))
	+dz*(
	    dri*(dphii*c1003+dphi*c1103)+
	    dr*(dphii*c1013+dphi*c1113)))
	+rg2/6*(
	    dzi*(
		dri3*(dphii*c0005+dphi*c0105)+
		dr3*(dphii*c0015+dphi*c0115))
	    +dz*(
		dri3*(dphii*c1005+dphi*c1105)+
		dr3*(dphii*c1015+dphi*c1115)))
	+phig2/6*(
	    dzi*(
		dri*(dphii3*c0006+dphi3*c0106)+
		dr*(dphii3*c0016+dphi3*c0116))
	    +dz*(
		dri*(dphii3*c1006+dphi3*c1106)+
		dr*(dphii3*c1016+dphi3*c1116)))
	+rg2*phig2/36*(
	    dzi*(
		dri3*(dphii3*c0007+dphi3*c0107)+
		dr3*(dphii3*c0017+dphi3*c0117))
	    +dz*(
		dri3*(dphii3*c1007+dphi3*c1107)+
		dr3*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/drdphi */
    B_dB[7][i] = rgi*phigi*(
	dzi*(
	    (c0000  -c0100)-
	    (c0010-c0110))
	+dz*(
	    (c1000  -c1100)-
	    (c1010-c1110)))
	+phigi*rg/6*(
	    dzi*(
		dri3dr*(-c0001+c0101)+
		dr3dr*(-c0011+c0111))
	    +dz*(
		dri3dr*(-c1001+c1101)+
		dr3dr*(-c1011+c1111)))
	+rgi*phig/6*(
	    dzi*(
		-(dphii3dphi*c0002+dphi3dphi*c0102)
		+(dphii3dphi*c0012+dphi3dphi*c0112))
	    +dz*(
		-(dphii3dphi*c1002+dphi3dphi*c1102)
		+(dphii3dphi*c1012+dphi3dphi*c1112)))
	+rgi*phigi*zg2/6*(
	    dzi3*(
		(c0003  -c0103)-
		(c0013-c0113))
	    +dz3*(
		(c1003  -c1103)-
		(c1013-c1113)))
	+rg*phig/36*(
	    dzi*(
		dri3dr*(dphii3dphi*c0004+dphi3dphi*c0104)+
		dr3dr*(dphii3dphi*c0014+dphi3dphi*c0114))
	    +dz*(
		dri3dr*(dphii3dphi*c1004+dphi3dphi*c1104)+
		dr3dr*(dphii3dphi*c1014+dphi3dphi*c1114)))
	+phigi*rg*zg2/36*(
	    dzi3*(
		dri3dr*(-c0005+c0105)+
		dr3dr*(-c0015+c0115))
	    +dz3*(
		dri3dr*(-c1005+c1105)+
		dr3dr*(-c1015+c1115)))
	+rgi*phig*zg2/36*(
	    dzi3*(
		-(dphii3dphi*c0006+dphi3dphi*c0106)
		+(dphii3dphi*c0016+dphi3dphi*c0116))
	    +dz3*(
		-(dphii3dphi*c1006+dphi3dphi*c1106)
		+(dphii3dphi*c1016+dphi3dphi*c1116)))
	+rg*phig*zg2/216*(
	    dzi3*(
		dri3dr*(dphii3dphi*c0007+dphi3dphi*c0107)+
		dr3dr*(dphii3dphi*c0017+dphi3dphi*c0117))
	    +dz3*(
		dri3dr*(dphii3dphi*c1007+dphi3dphi*c1107)+
		dr3dr*(dphii3dphi*c1017+dphi3dphi*c1117)));
    
    /* d2f/drdz */
    B_dB[8][i] = rgi*zgi*(
	(
	    (dphii*c0000+dphi*c0100) -
	    (dphii*c0010+dphi*c0110))
	-(
	    (dphii*c1000+dphi*c1100) -
	    (dphii*c1010+dphi*c1110)))
	+rg*zgi/6*(
	    -(
		dri3dr*(dphii*c0001+dphi*c0101)+
		dr3dr*(dphii*c0011+dphi*c0111))
	    +(
		dri3dr*(dphii*c1001+dphi*c1101)+
		dr3dr*(dphii*c1011+dphi*c1111)))
	+rgi*phig2*zgi/6*(
	    (
		(dphii3*c0002+dphi3*c0102) -
		(dphii3*c0012+dphi3*c0112))
	    -(
		(dphii3*c1002+dphi3*c1102) -
		(dphii3*c1012+dphi3*c1112)))
	+rgi*zg/6*(
	    dzi3dz*(
		-(dphii*c0003+dphi*c0103)
		+(dphii*c0013+dphi*c0113))
	    +dz3dz*(
		-(dphii*c1003+dphi*c1103)
		+(dphii*c1013+dphi*c1113)))
	+rg*phig2*zgi/36*(
	    -(
		dri3dr*(dphii3*c0004+dphi3*c0104)+
		dr3dr*(dphii3*c0014+dphi3*c0114))
	    +(
		dri3dr*(dphii3*c1004+dphi3*c1104)+
		dr3dr*(dphii3*c1014+dphi3*c1114)))
	+rg*zg/36*(
	    dzi3dz*(
		dri3dr*(dphii*c0005+dphi*c0105)+
		dr3dr*(dphii*c0015+dphi*c0115))
	    +dz3dz*(
		dri3dr*(dphii*c1005+dphi*c1105)+
		dr3dr*(dphii*c1015+dphi*c1115)))
	+rgi*phig2*zg/36*(
	    dzi3dz*(
		-(dphii3*c0006+dphi3*c0106)
		+(dphii3*c0016+dphi3*c0116))
	    +dz3dz*(
		-(dphii3*c1006+dphi3*c1106)
		+(dphii3*c1016+dphi3*c1116)))
	+rg*phig2*zg/216*(
	    dzi3dz*(
		dri3dr*(dphii3*c0007+dphi3*c0107)+
		dr3dr*(dphii3*c0017+dphi3*c0117))
	    +dz3dz*(
		dri3dr*(dphii3*c1007+dphi3*c1107)+
		dr3dr*(dphii3*c1017+dphi3*c1117)));
    
    /* d2f/dphidz */
    B_dB[9][i] = phigi*zgi*(
	(
	    dri*(c0000  -c0100)+
	    dr*(c0010-c0110))
	-(
	    dri*(c1000  -c1100)+
	    dr*(c1010-c1110)))
	+phigi*rg2*zgi/6*(
	    (
		dri3*(c0001  -c0101)+
		dr3*(c0011-c0111))
	    -(
		dri3*(c1001  -c1101)+
		dr3*(c1011-c1111)))
	+phig*zgi/6*(
	    -(
		dri*(dphii3dphi*c0002+dphi3dphi*c0102)+
		dr*(dphii3dphi*c0012+dphi3dphi*c0112))
	    +(
		dri*(dphii3dphi*c1002+dphi3dphi*c1102)+
		dr*(dphii3dphi*c1012+dphi3dphi*c1112)))
	+phigi*zg/6*(
	    dzi3dz*(
		dri*(-c0003+c0103)+
		dr*(-c0013+c0113))
	    +dz3dz*(
		dri*(-c1003+c1103)+
		dr*(-c1013+c1113)))
	+rg2*phig*zgi/36*(
	    -(
		dri3*(dphii3dphi*c0004+dphi3dphi*c0104)+
		dr3*(dphii3dphi*c0014+dphi3dphi*c0114))
	    +(
		dri3*(dphii3dphi*c1004+dphi3dphi*c1104)+
		dr3*(dphii3dphi*c1014+dphi3dphi*c1114)))
	+phigi*rg2*zg/36*(
	    dzi3dz*(
		dri3*(-c0005+c0105)+
		dr3*(-c0015+c0115))
	    +dz3dz*(
		dri3*(-c1005+c1105)+
		dr3*(-c1015+c1115)))
	+phig*zg/36*(
	    dzi3dz*(
		dri*(dphii3dphi*c0006+dphi3dphi*c0106)+
		dr*(dphii3dphi*c0016+dphi3dphi*c0116))
	    +dz3dz*(
		dri*(dphii3dphi*c1006+dphi3dphi*c1106)+
		dr*(dphii3dphi*c1016+dphi3dphi*c1116)))
	+rg2*phig*zg/216*(
	    dzi3dz*(
		dri3*(dphii3dphi*c0007+dphi3dphi*c0107)+
		dr3*(dphii3dphi*c0017+dphi3dphi*c0117))
	    +dz3dz*(
		dri3*(dphii3dphi*c1007+dphi3dphi*c1107)+
		dr3*(dphii3dphi*c1017+dphi3dphi*c1117)));

    }

    return err;
}

/**
 * @brief Free allocated memory in interpolation data struct
 *
 * This function frees the memory allocated for interpolation coefficients
 * in the interpolation data struct
 * 
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 */
void interp3Dcomp_free(interp3D_data* str) {
    free(str->c);
}
