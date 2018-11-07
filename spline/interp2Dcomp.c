/**
 * @file interp2Dcomp.c
 * @brief Bicubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <stdio.h> /* Needed for printf debugging purposes */
#include "../ascot5.h"
#include "interp2Dcomp.h"
#include "spline1Dcomp.h"

/**
 * @brief Calculate bicubic spline interpolation coefficients for scalar 2D data
 *
 * This function calculates the bicubic spline interpolation coefficients for
 * the given data and stores them in the data struct. Compact  cofficients are
 * calculated directly.
 * 
 * @todo Error checking
 *
 * @param str data struct for data interpolation
 * @param f 2D data to be interpolated
 * @param n_r number of data points in the r direction
 * @param n_z number of data points in the z direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 */
int interp2Dcomp_init(interp2D_data* str, real* f, int n_r, int n_z,
		      real r_min, real r_max, real r_grid,
		      real z_min, real z_max, real z_grid) {
    int err = 0;

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;
    str->c = malloc(n_z*n_r*4*sizeof(real));

    /* Declare and allocate the needed variables */
    int i_r;                                 /**< index for r variable */
    int i_z;                                 /**< index for z variable */
    real* f_r = malloc(n_r*sizeof(real));    /**< Temporary array for data along r */
    real* f_z = malloc(n_z*sizeof(real));    /**< Temporary array for data along z */
    real* c_r = malloc(n_r*2*sizeof(real));  /**< Temp array for coefficients along r */
    real* c_z = malloc(n_z*2*sizeof(real));  /**< Temp array for coefficients along z */

    if(f_r == NULL || f_z == NULL || c_r == NULL || c_z == NULL) {
	err = 1;
    }
    else {
	/* Bicubic spline surface over rz-grid. Note how we account for normalized grid. */
	/* Cubic spline along r for each z to get frr */
	for(i_z=0; i_z<n_z; i_z++) {
	    for(i_r=0; i_r<n_r; i_r++) {
		f_r[i_r] = f[i_z*n_r+i_r];
	    }
	    spline1Dcomp(f_r,n_r,0,c_r);
	    for(i_r=0; i_r<n_r; i_r++) {
		str->c[i_z*n_r*4+i_r*4] = c_r[i_r*2];
		str->c[i_z*n_r*4+i_r*4+1] = c_r[i_r*2+1]/(r_grid*r_grid);
	    }
	}

	/* Two cubic splines along z for each r using f and frr */
	for(i_r=0; i_r<n_r; i_r++) {
	    /* fzz */
	    for(i_z=0; i_z<n_z; i_z++) {
		f_z[i_z] =  f[i_z*n_r+i_r];
	    }
	    spline1Dcomp(f_z,n_z,0,c_z);
	    for(i_z=0; i_z<n_z; i_z++) {
		str->c[i_z*n_r*4+i_r*4+2] = c_z[i_z*2+1]/(z_grid*z_grid);
	    }
	    /* frrzz */
	    for(i_z=0; i_z<n_z; i_z++) {
		f_z[i_z] =  str->c[i_z*n_r*4+i_r*4+1];
	    }
	    spline1Dcomp(f_z,n_z,0,c_z);
	    for(i_z=0; i_z<n_z; i_z++) {
		str->c[i_z*n_r*4+i_r*4+3] = c_z[i_z*2+1]/(z_grid*z_grid);
	    }
	}
    }

    /* Free allocated memory */
    free(f_r);
    free(f_z);
    free(c_r);
    free(c_z);

    return err;
}

/**
 * @brief Evaluate interpolated value of 2D scalar field
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bicubic spline interpolation coefficients of the compact form.
 * 
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param z z-coordinate
 */
integer interp2Dcomp_eval_B(real* B, interp2D_data* str, real r, real z) {
    int i_r = (r-str->r_min)/str->r_grid; /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dr3 = dr*(dr*dr-1.0);
    real dri = 1.0-dr;
    real dri3 = dri*(dri*dri-1.0);
    real rg2 = str->r_grid*str->r_grid;        /**< Square of cell length in r direction */
    int i_z = (z-str->z_min)/str->z_grid;                   /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */
    real dz3 = dz*(dz*dz-1.0);
    real dzi = 1.0-dz;
    real dzi3 = dzi*(dzi*dzi-1.0);
    real zg2 = str->z_grid*str->z_grid; /**< Square of cell length in z direction */
    int n = i_z*str->n_r*4+i_r*4;       /**< Index jump to cell */
    int r1 = 4;                         /**< Index jump one r forward */
    int z1 = str->n_r*4;                /**< Index jump one z forward */

    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {
	*B = (
	    dri*(dzi*str->c[n]   +dz*str->c[n+z1])+
	    dr*(dzi*str->c[n+r1]+dz*str->c[n+z1+r1]))
	    +rg2/6.0*(
		dri3*(dzi*str->c[n+1]   +dz*str->c[n+z1+1])+
		dr3*(dzi*str->c[n+r1+1]+dz*str->c[n+z1+r1+1]))
	    +zg2/6.0*(
		dri*(dzi3*str->c[n+2]   +dz3*str->c[n+z1+2])+
		dr*(dzi3*str->c[n+r1+2]+dz3*str->c[n+z1+r1+2]))
	    +rg2*zg2/36.0*(
		dri3*(dzi3*str->c[n+3]   +dz3*str->c[n+z1+3])+
		dr3*(dzi3*str->c[n+r1+3]+dz3*str->c[n+z1+r1+3]));
    }
    return err;
}

/**
 * @brief Evaluate interpolated value of 2D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 2D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 * 
 * @todo Check discrepency to ascot4 and explicit version
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param z z-coordinate
 */
integer interp2Dcomp_eval_dB(real* B_dB, interp2D_data* str, real r, real z) {
    int i_r = (r-str->r_min)/str->r_grid;                   /**< index for r variable */
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dr3 = dr*(dr*dr-1);
    real dr3dr = 3*dr*dr-1;           /**< r-derivative of dr3, not including 1/r_grid */
    real dri = 1.0-dr;
    real dri3 = dri*(dri*dri-1);
    real dri3dr = -3*dri*dri+1;       /**< r-derivative of dri3, not including 1/r_grid */
    real rg = str->r_grid;            /**< Cell length in r direction */
    real rg2 = rg*rg;
    real rgi = 1.0/rg;
    int i_z = (z-str->z_min)/str->z_grid; /**< index for z variable */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */
    real dz3 = dz*(dz*dz-1);
    real dz3dz = 3*dz*dz-1;           /**< z-derivative of dz3, not including 1/z_grid */
    real dzi = 1.0-dz;
    real dzi3 = dzi*(dzi*dzi-1);
    real dzi3dz = -3*dzi*dzi+1;       /**< z-derivative of dzi3, not including 1/z_grid */
    real zg = str->z_grid;            /**< Cell length in z direction */
    real zg2 = zg*zg;
    real zgi = 1.0/zg;
    int n = i_z*str->n_r*4+i_r*4;     /**< Index jump to cell */
    int r1 = 4;                       /**< Index jump one r forward */
    int z1 = str->n_r*4;              /**< Index jump one z forward */

    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {
	/* f */
	B_dB[0] = (
	    dri*(dzi*str->c[n]+dz*str->c[n+z1])+
	    dr*(dzi*str->c[n+r1]+dz*str->c[n+z1+r1]))
	    +(rg2/6)*(
		dri3*(dzi*str->c[n+1] + dz*str->c[n+z1+1])+
		dr3*(dzi*str->c[n+r1+1] + dz*str->c[n+z1+r1+1]))
	    +(zg2/6)*(
		dri*(dzi3*str->c[n+2]+dz3*str->c[n+z1+2])+
		dr*(dzi3*str->c[n+r1+2]+dz3*str->c[n+z1+r1+2]))
	    +(rg2*zg2/36)*(
		dri3*(dzi3*str->c[n+3]+dz3*str->c[n+z1+3])+
		dr3*(dzi3*str->c[n+r1+3]+dz3*str->c[n+z1+r1+3]));

	/* df/dr */
	B_dB[1] = rgi*(
	    -(dzi*str->c[n]  +dz*str->c[n+z1])
	    +(dzi*str->c[n+r1]+dz*str->c[n+z1+r1]))
	    +(rg/6)*(
		dri3dr*(dzi*str->c[n+1]  +dz*str->c[n+z1+1])+
		dr3dr*(dzi*str->c[n+r1+1]+dz*str->c[n+z1+r1+1]))
	    +(rgi*zg2/6)*(
		-(dzi3*str->c[n+2]  +dz3*str->c[n+z1+2])
		+(dzi3*str->c[n+r1+2]+dz3*str->c[n+z1+r1+2]))
	    +(rg*zg2/36)*(
		dri3dr*(dzi3*str->c[n+3]  +dz3*str->c[n+z1+3])+
		dr3dr*(dzi3*str->c[n+r1+3]+dz3*str->c[n+z1+r1+3]));

	/* df/dz */
	B_dB[2] = zgi*(
	    dri*(-str->c[n]  +str->c[n+z1])+
	    dr*(-str->c[n+r1]+str->c[n+z1+r1]))
	    +(rg2*zgi/6)*(
		dri3*(-str->c[n+1]  +str->c[n+z1+1])+
		dr3*(-str->c[n+r1+1]+str->c[n+z1+r1+1]))
	    +(zg/6)*(
		dri*(dzi3dz*str->c[n+2]  +dz3dz*str->c[n+z1+2])+
		dr*(dzi3dz*str->c[n+r1+2]+dz3dz*str->c[n+z1+r1+2]))
	    +(rg2*zg/36)*(
		dri3*(dzi3dz*str->c[n+3]  +dz3dz*str->c[n+z1+3])+
		dr3*(dzi3dz*str->c[n+r1+3]+dz3dz*str->c[n+z1+r1+3]));

	/* d2f/dr2 */
	B_dB[3] = (
	    dri*(dzi*str->c[n+1]  +dz*str->c[n+z1+1])+
	    dr*(dzi*str->c[n+r1+1]+dz*str->c[n+z1+r1+1]))
	    +(zg2/6)*(
		dri*(dzi3*str->c[n+3]  +dz3*str->c[n+z1+3])+
		dr*(dzi3*str->c[n+r1+3]+dz3*str->c[n+z1+r1+3]));

	/* d2f/dz2 */
	B_dB[4] = (
	      dri*(dzi*str->c[n+2]  +dz*str->c[n+z1+2])+
	      dr*(dzi*str->c[n+r1+2]+dz*str->c[n+z1+r1+2]))
	+rg2/6*(
	    dri3*(dzi*str->c[n+3]  +dz*str->c[n+z1+3])+
	    dr3*(dzi*str->c[n+r1+3]+dz*str->c[n+z1+r1+3]));

	/* d2f/dzdr */
	B_dB[5] = rgi*zgi*(
	    str->c[n]  -str->c[n+z1]
	    -str->c[n+r1]+str->c[n+z1+r1])
	    +(rg/6*zgi)*(
		dri3dr*(-str->c[n+1]  +str->c[n+z1+1])+
		dr3dr*(-str->c[n+r1+1]+str->c[n+z1+r1+1]))
	    +(rgi/6*zg)*(
		-(dzi3dz*str->c[n+2]  +dz3dz*str->c[n+z1+2])
		+(dzi3dz*str->c[n+r1+2]+dz3dz*str->c[n+z1+r1+2]))
	    +(rg*zg/36)*(
		dri3dr*(dzi3dz*str->c[n+3]  +dz3dz*str->c[n+z1+3])+
		dr3dr*(dzi3dz*str->c[n+r1+3]+dz3dz*str->c[n+z1+r1+3]));
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
void interp2Dcomp_free(interp2D_data* str) {
    free(str->c);
}
