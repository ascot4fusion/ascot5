/**
 * @file interp2Dexpl.c
 * @brief Bicubic spline interpolation, i.e. cubic spline interpolation of 2D scalar data
 */
#include <stdlib.h>
#include <stdio.h> /* Needed for printf debugging purposes */
#include "../ascot5.h"
#include "interp2Dexpl.h"
#include "interp2D.h"
#include "spline1D.h"

/**
 * @brief Calculate bicubic spline interpolation coefficients for scalar 2D data
 *
 * This function calculates the bicubic spline interpolation coefficients for
 * the given data and stores them in the data struct. The explicit cofficients
 * are calculated and stored.
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
void interp2Dexpl_init(interp2D_data* str, real* f, int n_r, int n_z,
		       real r_min, real r_max, real r_grid,
		       real z_min, real z_max, real z_grid) {

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;//(r_max-r_min)/(n_r-1);
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;//(z_max-z_min)/(n_z-1);
    str->c = malloc(n_z*n_r*16*sizeof(real));

    /* Declare and allocate the needed variables */
    int i_r;
    int i_z;
    int i_c;
    real* f_r = malloc(n_r*sizeof(real));
    real* f_z = malloc(n_z*sizeof(real));
    real* c_r = malloc((n_r-1)*4*sizeof(real));
    real* c_z = malloc((n_z-1)*4*sizeof(real));
    int i_s;
    int i_ct;

    /* Bicubic spline surface over rz-grid */
    /* Cubic spline along r for each z */
    for(i_z=0; i_z<n_z; i_z++) {
	for(i_r=0; i_r<n_r; i_r++) {
	    f_r[i_r] = f[i_z*n_r+i_r];
	}
	spline1D(f_r,n_r,0,c_r);
	for(i_r=0; i_r<n_r-1; i_r++) {
	    for(i_c=0; i_c<4; i_c++) {
		i_ct = i_c;
		str->c[i_z*n_r*16+i_r*16+i_c] = c_r[i_r*4+i_ct];
	    }
	}
    }
    
    /* Four cubic splines along z for each r using four different data sets */
    for(i_r=0; i_r<n_r-1; i_r++) {
	/* s0, s1, s2, s3 */
	for(i_s=0; i_s<4; i_s++) {
	    for(i_z=0; i_z<n_z; i_z++) {
		f_z[i_z] = str->c[i_z*n_r*16+i_r*16+i_s];
	    }
	    spline1D(f_z,n_z,0,c_z);
	    for(i_z=0; i_z<n_z-1; i_z++) {
		i_ct = 0;
		for(i_c=i_s; i_c<16; i_c=i_c+4) {
		    str->c[i_z*n_r*16+i_r*16+i_c] = c_z[i_z*4+i_ct];
		    i_ct++;
		}
	    }
	}
    }

    /* Free allocated memory */
    free(f_r);
    free(f_z);
    free(c_r);
    free(c_z);
}

/**
 * @brief Evaluate interpolated value of 2D scalar field
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bicubic spline interpolation coefficients of the explicit form.
 * 
 * @todo Check if discrepency to ascot4. Maybe different bc or because 2D from psi.
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param z z-coordinate
 */
void interp2Dexpl_eval_B(real* B, interp2D_data* str, real r, real z) {
    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dr2 = dr*dr;
    real dr3 = dr2*dr;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dz2 = dz*dz;
    real dz3 = dz2*dz;
    int n = i_z*str->n_r*16+i_r*16;

    *B =      str->c[n+ 0]+dr*str->c[n+ 1]+dr2*str->c[n+ 2]+dr3*str->c[n+ 3]
	 +dz*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7])
	+dz2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
	+dz3*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]);
}

/**
 * @brief Evaluate interpolated value of 2D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 2D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the explicit form.
 * 
 * @todo Check discrepency to ascot4. Maybe different bc or because 2D from psi.
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param z z-coordinate
 */
void interp2Dexpl_eval_dB(real* B_dB, interp2D_data* str, real r, real z) {
    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dr2 = dr*dr;
    real dr3 = dr2*dr;
    real rgi = 1.0/str->r_grid;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dz2 = dz*dz;
    real dz3 = dz2*dz;
    real zgi = 1.0/str->r_grid;
    int n = i_z*str->n_r*16+i_r*16;

    /* f */
    B_dB[0] = str->c[n+ 0]+dr*str->c[n+ 1]+dr2*str->c[n+ 2]+dr3*str->c[n+ 3]
	 +dz*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7])
	+dz2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
	+dz3*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]);

    /* df/dr */
    B_dB[1] = rgi*(str->c[n+ 1]+2*dr*str->c[n+ 2]+3*dr2*str->c[n+ 3]
		   +dz*(str->c[n+ 5]+2*dr*str->c[n+ 6]+3*dr2*str->c[n+ 7])
		   +dz2*(str->c[n+ 9]+2*dr*str->c[n+10]+3*dr2*str->c[n+11])
		   +dz3*(str->c[n+13]+2*dr*str->c[n+14]+3*dr2*str->c[n+15]));

    /* df/dz */
    B_dB[2] = zgi*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7]
		   +2*dz*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
		   +3*dz2*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]));

    /* d2f/dr^2 */
    B_dB[3] = rgi*rgi*(2*str->c[n+ 2]+6*dr*str->c[n+ 3]
		       +dz*(2*str->c[n+ 6]+6*dr*str->c[n+ 7])
		       +dz2*(2*str->c[n+10]+6*dr*str->c[n+11])
		       +dz3*(2*str->c[n+14]+6*dr*str->c[n+15]));

    /* d2f/dz^2 */
    B_dB[4] = zgi*zgi*(2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
		       +6*dz*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]));

    /* d2f/dzdr */
    B_dB[5] = rgi*zgi*(str->c[n+ 5]+2*dr*str->c[n+ 6]+3*dr2*str->c[n+ 7]
		       +2*dz*(str->c[n+ 9]+2*dr*str->c[n+10]+3*dr2*str->c[n+11])
		       +3*dz2*(str->c[n+13]+2*dr*str->c[n+14]+3*dr2*str->c[n+15]));
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
void interp2Dexpl_free(interp2D_data* str) {
    free(str->c);
}
