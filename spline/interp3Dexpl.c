/**
 * @file interp3Dexpl.c
 * @brief Tricubic spline interpolation, i.e. cubic spline interpolation of 3D scalar data
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "interp3Dexpl.h"
#include "interp3D.h"
#include "spline1D.h"

/**
 * @brief Calculate tricubic spline interpolation coefficients for scalar 3D data
 *
 * This function calculates the tricubic spline interpolation coefficients for
 * the given data and stores them in the data struct. The explicit cofficients
 * are calculated and stored.
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
void interp3Dexpl_init(interp3D_data* str, real* f, int n_r, int n_phi, int n_z,
		       real r_min, real r_max, real r_grid,
		       real phi_min, real phi_max, real phi_grid,
		       real z_min, real z_max, real z_grid) {

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_phi = n_phi;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;//(r_max-r_min)/(n_r-1);
    str->phi_min = phi_min;
    str->phi_max = phi_max;
    str->phi_grid = phi_grid;//(phi_max-phi_min)/n_phi;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;//(z_max-z_min)/(n_z-1);
    str->c = malloc(n_phi*n_z*n_r*64*sizeof(real));

    /* Declare and allocate the needed variables */
    int i_r;
    int i_phi;
    int i_z;
    int i_c;
    real* f_r = malloc(n_r*sizeof(real));
    real* f_phi = malloc(n_phi*sizeof(real));
    real* f_z = malloc(n_z*sizeof(real));
    real* c_r = malloc((n_r-1)*4*sizeof(real));
    real* c_phi = malloc(n_phi*4*sizeof(real));
    real* c_z = malloc((n_z-1)*4*sizeof(real));
    int i_s;
    int i_ss;
    int i_ct;

    /* Bicubic spline surface over rz-grid for each phi */
    for(i_phi=0; i_phi<n_phi; i_phi++) {
	/* Cubic spline along r for each z */
	for(i_z=0; i_z<n_z; i_z++) {
	    for(i_r=0; i_r<n_r; i_r++) {
		f_r[i_r] = f[i_phi*n_z*n_r+i_z*n_r+i_r];
	    }
	    spline1D(f_r,n_r,0,c_r);
	    for(i_r=0; i_r<n_r-1; i_r++) {
		for(i_c=0; i_c<4; i_c++) {
		    i_ct = i_c;
		    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] = c_r[i_r*4+i_ct];
		}
	    }
	}
	
	/* Four cubic splines along z for each r using four different data sets */
	for(i_r=0; i_r<n_r-1; i_r++) {
	    /* s0, s1, s2, s3 */
	    for(i_s=0; i_s<4; i_s++) {
		for(i_z=0; i_z<n_z; i_z++) {
		    f_z[i_z] = str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_s];
		}
		spline1D(f_z,n_z,0,c_z);
		for(i_z=0; i_z<n_z-1; i_z++) {
		    i_ct = 0;
		    for(i_c=i_s; i_c<16; i_c=i_c+4) {
			str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] = c_z[i_z*4+i_ct];
			i_ct++;
		    }
		}
	    }
	}
    }

    /* Cubic spline along phi for each rz-pair to find the coefficients
       of the tricubic spline volume */
    for(i_z=0; i_z<n_z-1; i_z++) {
	for(i_r=0; i_r<n_r-1; i_r++) {
	    for(i_ss=0; i_ss<4; i_ss++) {
		for(i_s=0; i_s<4; i_s++) {
		    for(i_phi=0; i_phi<n_phi; i_phi++) {
			f_phi[i_phi] = str->c[i_phi*n_z*n_r*64+i_z*n_r*64
					    +i_r*64+(i_ss*4+i_s)];
		    }
		    spline1D(f_phi,n_phi,1,c_phi);
		    for(i_phi=0; i_phi<n_phi; i_phi++) {
			i_ct = 0;
			for(i_c=4*i_ss+i_s; i_c<64; i_c=i_c+16) {
			    str->c[i_phi*n_z*n_r*64+i_z*n_r*64+i_r*64+i_c] 
				= c_phi[i_phi*4+i_ct];
			    i_ct++;
			}
		    }
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
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * tricubic spline interpolation coefficients of the explicit form.
 * 
 * @todo Check if discrepency to ascot4. Maybe different bc or because 2D from psi.
 * @todo Error checking
 *
 * @param B variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
void interp3Dexpl_eval_B(real* B, interp3D_data* str, real r, real phi, real z) {

    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI - phi;}

    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dr2 = dr*dr;
    real dr3 = dr2*dr;
    int i_phi = (phi-str->phi_min)/str->phi_grid;
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid;
    real dphi2 = dphi*dphi;
    real dphi3 = dphi2*dphi;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dz2 = dz*dz;
    real dz3 = dz2*dz;
    int n = i_phi*str->n_z*str->n_r*64+i_z*str->n_r*64+i_r*64;

    *B = (
	             str->c[n+ 0]+dr*str->c[n+ 1]+dr2*str->c[n+ 2]+dr3*str->c[n+ 3]
	        +dz*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7])
	       +dz2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
	       +dz3*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]))
	 +dphi*(
	             str->c[n+16]+dr*str->c[n+17]+dr2*str->c[n+18]+dr3*str->c[n+19]
	        +dz*(str->c[n+20]+dr*str->c[n+21]+dr2*str->c[n+22]+dr3*str->c[n+23])
	       +dz2*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27])
	       +dz3*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31]))
    	+dphi2*(
		      str->c[n+32]+dr*str->c[n+33]+dr2*str->c[n+34]+dr3*str->c[n+35]
		 +dz*(str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39])
    	        +dz2*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
    	        +dz3*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
    	+dphi3*(
		      str->c[n+48]+dr*str->c[n+49]+dr2*str->c[n+50]+dr3*str->c[n+51]
		 +dz*(str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55])
    	        +dz2*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
    	        +dz3*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63]));
}

/**
 * @brief Evaluate interpolated value of 3D scalar field and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 3D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the explicit form.
 * 
 * @todo Check if discrepency to ascot4. Maybe different bc or because 2D from psi.
 * @todo Error checking
 *
 * @param B_dB array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
void interp3Dexpl_eval_dB(real* B_dB, interp3D_data* str, real r, real phi, real z) {
    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI - phi;}

    int i_r = (r-str->r_min)/str->r_grid;
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid;
    real dr2 = dr*dr;
    real dr3 = dr2*dr;
    real rgi = 1.0/str->r_grid;
    int i_phi = (phi-str->phi_min)/str->phi_grid;
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid;
    real dphi2 = dphi*dphi;
    real dphi3 = dphi2*dphi;
    real phigi = 1.0/str->phi_grid;
    int i_z = (z-str->z_min)/str->z_grid;
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid;
    real dz2 = dz*dz;
    real dz3 = dz2*dz;
    real zgi = 1.0/str->z_grid;
    int n = i_phi*str->n_z*str->n_r*64+i_z*str->n_r*64+i_r*64;

    /* f */
    B_dB[0] = (
	              str->c[n+ 0]+dr*str->c[n+ 1]+dr2*str->c[n+ 2]+dr3*str->c[n+ 3]
	         +dz*(str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7])
	        +dz2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
	        +dz3*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]))
	 +dphi*(
	              str->c[n+16]+dr*str->c[n+17]+dr2*str->c[n+18]+dr3*str->c[n+19]
	         +dz*(str->c[n+20]+dr*str->c[n+21]+dr2*str->c[n+22]+dr3*str->c[n+23])
	        +dz2*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27])
	        +dz3*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31]))
    	+dphi2*(
		      str->c[n+32]+dr*str->c[n+33]+dr2*str->c[n+34]+dr3*str->c[n+35]
		 +dz*(str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39])
    	        +dz2*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
    	        +dz3*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
    	+dphi3*(
		      str->c[n+48]+dr*str->c[n+49]+dr2*str->c[n+50]+dr3*str->c[n+51]
		 +dz*(str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55])
    	        +dz2*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
    	        +dz3*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63]));

    /* df/dr */
    B_dB[1] = rgi*( 
	       (
                      str->c[n+ 1]+2*dr*str->c[n+ 2]+3*dr2*str->c[n+ 3]
                 +dz*(str->c[n+ 5]+2*dr*str->c[n+ 6]+3*dr2*str->c[n+ 7])
                +dz2*(str->c[n+ 9]+2*dr*str->c[n+10]+3*dr2*str->c[n+11])
                +dz3*(str->c[n+13]+2*dr*str->c[n+14]+3*dr2*str->c[n+15]))
         +dphi*(
                      str->c[n+17]+2*dr*str->c[n+18]+3*dr2*str->c[n+19]
                 +dz*(str->c[n+21]+2*dr*str->c[n+22]+3*dr2*str->c[n+23])
                +dz2*(str->c[n+25]+2*dr*str->c[n+26]+3*dr2*str->c[n+27])
                +dz3*(str->c[n+29]+2*dr*str->c[n+30]+3*dr2*str->c[n+31]))
        +dphi2*(
                      str->c[n+33]+2*dr*str->c[n+34]+3*dr2*str->c[n+35]
                 +dz*(str->c[n+37]+2*dr*str->c[n+38]+3*dr2*str->c[n+39])
                +dz2*(str->c[n+41]+2*dr*str->c[n+42]+3*dr2*str->c[n+43])
                +dz3*(str->c[n+45]+2*dr*str->c[n+46]+3*dr2*str->c[n+47]))
        +dphi3*(
                      str->c[n+49]+2*dr*str->c[n+50]+3*dr2*str->c[n+51]
                 +dz*(str->c[n+53]+2*dr*str->c[n+54]+3*dr2*str->c[n+55])
                +dz2*(str->c[n+57]+2*dr*str->c[n+58]+3*dr2*str->c[n+59])
		+dz3*(str->c[n+61]+2*dr*str->c[n+62]+3*dr2*str->c[n+63])));

    /* df/dphi */
    B_dB[2] = phigi*(
	         (
                        str->c[n+16]+dr*str->c[n+17]+dr2*str->c[n+18]+dr3*str->c[n+19]
                   +dz*(str->c[n+20]+dr*str->c[n+21]+dr2*str->c[n+22]+dr3*str->c[n+23])
                  +dz2*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27])
                  +dz3*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31]))
         +2*dphi*(
                        str->c[n+32]+dr*str->c[n+33]+dr2*str->c[n+34]+dr3*str->c[n+35]
                   +dz*(str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39])
                  +dz2*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
                  +dz3*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
        +3*dphi2*(
                        str->c[n+48]+dr*str->c[n+49]+dr2*str->c[n+50]+dr3*str->c[n+51]
                   +dz*(str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55])
                  +dz2*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
		  +dz3*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63])));

    /* df/dz */
    B_dB[3] = zgi*(
	      (
                        str->c[n+ 4]+dr*str->c[n+ 5]+dr2*str->c[n+ 6]+dr3*str->c[n+ 7]
                 +2*dz*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
                +3*dz2*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]))
        +dphi*(
                        str->c[n+20]+dr*str->c[n+21]+dr2*str->c[n+22]+dr3*str->c[n+23]
                 +2*dz*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27])
                +3*dz2*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31]))
        +dphi2*(
                        str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39]
                 +2*dz*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
                +3*dz2*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
        +dphi3*(
                        str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55]
                 +2*dz*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
		+3*dz2*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63])));

    /* d2f/dr^2 */
    B_dB[4] =  rgi*rgi*(
	      (
                      2*str->c[n+ 2]+6*dr*str->c[n+ 3]
		 +dz*(2*str->c[n+ 6]+6*dr*str->c[n+ 7])
		+dz2*(2*str->c[n+10]+6*dr*str->c[n+11])
		+dz3*(2*str->c[n+14]+6*dr*str->c[n+15]))
	 +dphi*(
                      2*str->c[n+18]+6*dr*str->c[n+19]
		 +dz*(2*str->c[n+22]+6*dr*str->c[n+23])
		+dz2*(2*str->c[n+26]+6*dr*str->c[n+27])
		+dz3*(2*str->c[n+30]+6*dr*str->c[n+31]))
        +dphi2*(
                      2*str->c[n+34]+6*dr*str->c[n+35]
		 +dz*(2*str->c[n+38]+6*dr*str->c[n+39])
		+dz2*(2*str->c[n+42]+6*dr*str->c[n+43])
		+dz3*(2*str->c[n+46]+6*dr*str->c[n+47]))
        +dphi3*(
                      2*str->c[n+50]+6*dr*str->c[n+51]
		 +dz*(2*str->c[n+54]+6*dr*str->c[n+55])
		+dz2*(2*str->c[n+58]+6*dr*str->c[n+59])
		      +dz3*(2*str->c[n+62]+6*dr*str->c[n+63])));

    /* d2f/dphi^2 */
    B_dB[5] = phigi*phigi*(
	      2*(
                       str->c[n+32]+dr*str->c[n+33]+dr2*str->c[n+34]+dr3*str->c[n+35]
		  +dz*(str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39])
		 +dz2*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
		 +dz3*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
        +6*dphi*(
                       str->c[n+48]+dr*str->c[n+49]+dr2*str->c[n+50]+dr3*str->c[n+51]
       		  +dz*(str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55])
		 +dz2*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
		 +dz3*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63])));

    /* d2f/dz^2 */
    B_dB[6] = zgi*zgi*(
	       (
                    2*(str->c[n+ 8]+dr*str->c[n+ 9]+dr2*str->c[n+10]+dr3*str->c[n+11])
	        +6*dz*(str->c[n+12]+dr*str->c[n+13]+dr2*str->c[n+14]+dr3*str->c[n+15]))
         +dphi*(
                    2*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27])
	        +6*dz*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31]))
        +dphi2*(
                    2*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
		+6*dz*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
        +dphi3*(
                    2*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
	        +6*dz*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63])));

    /* d2f/dphidr */
    B_dB[7] = rgi*phigi*(
		 (
                        str->c[n+17]+2*dr*str->c[n+18]+3*dr2*str->c[n+19]
		   +dz*(str->c[n+21]+2*dr*str->c[n+22]+3*dr2*str->c[n+23])
		  +dz2*(str->c[n+25]+2*dr*str->c[n+26]+3*dr2*str->c[n+27])
		  +dz3*(str->c[n+29]+2*dr*str->c[n+30]+3*dr2*str->c[n+31]))
         +2*dphi*(
                        str->c[n+33]+2*dr*str->c[n+34]+3*dr2*str->c[n+35]
		   +dz*(str->c[n+37]+2*dr*str->c[n+38]+3*dr2*str->c[n+39])
		  +dz2*(str->c[n+41]+2*dr*str->c[n+42]+3*dr2*str->c[n+43])
		  +dz3*(str->c[n+45]+2*dr*str->c[n+46]+3*dr2*str->c[n+47]))
        +3*dphi2*(
                        str->c[n+49]+2*dr*str->c[n+50]+3*dr2*str->c[n+51]
		   +dz*(str->c[n+53]+2*dr*str->c[n+54]+3*dr2*str->c[n+55])
		  +dz2*(str->c[n+57]+2*dr*str->c[n+58]+3*dr2*str->c[n+59])
		  +dz3*(str->c[n+61]+2*dr*str->c[n+62]+3*dr2*str->c[n+63])));

    /* d2f/dzdr */
    B_dB[8] =  rgi*zgi*(
	       (
                        str->c[n+ 5]+2*dr*str->c[n+ 6]+3*dr2*str->c[n+ 7]
		 +2*dz*(str->c[n+ 9]+2*dr*str->c[n+10]+3*dr2*str->c[n+11])
		+3*dz2*(str->c[n+13]+2*dr*str->c[n+14]+3*dr2*str->c[n+15]))
	 +dphi*(
                        str->c[n+21]+2*dr*str->c[n+22]+3*dr2*str->c[n+23]
		 +2*dz*(str->c[n+25]+2*dr*str->c[n+26]+3*dr2*str->c[n+27])
	        +3*dz2*(str->c[n+29]+2*dr*str->c[n+30]+3*dr2*str->c[n+31]))
        +dphi2*(
                        str->c[n+37]+2*dr*str->c[n+38]+3*dr2*str->c[n+39]
		 +2*dz*(str->c[n+41]+2*dr*str->c[n+42]+3*dr2*str->c[n+43])
		+3*dz2*(str->c[n+45]+2*dr*str->c[n+46]+3*dr2*str->c[n+47]))
        +dphi3*(
                        str->c[n+53]+2*dr*str->c[n+54]+3*dr2*str->c[n+55]
		 +2*dz*(str->c[n+57]+2*dr*str->c[n+58]+3*dr2*str->c[n+59])
			+3*dz2*(str->c[n+61]+2*dr*str->c[n+62]+3*dr2*str->c[n+63])));

    /* d2f/dzdphi */
    B_dB[9] = phigi*zgi*(
		 (
                          str->c[n+20]+dr*str->c[n+21]+dr2*str->c[n+22]+dr3*str->c[n+23]
		   +2*dz*(str->c[n+24]+dr*str->c[n+25]+dr2*str->c[n+26]+dr3*str->c[n+27])
		  +3*dz2*(str->c[n+28]+dr*str->c[n+29]+dr2*str->c[n+30]+dr3*str->c[n+31]))
	 +2*dphi*(
                          str->c[n+36]+dr*str->c[n+37]+dr2*str->c[n+38]+dr3*str->c[n+39]
		   +2*dz*(str->c[n+40]+dr*str->c[n+41]+dr2*str->c[n+42]+dr3*str->c[n+43])
		  +3*dz2*(str->c[n+44]+dr*str->c[n+45]+dr2*str->c[n+46]+dr3*str->c[n+47]))
        +3*dphi2*(
                         str->c[n+52]+dr*str->c[n+53]+dr2*str->c[n+54]+dr3*str->c[n+55]
		  +2*dz*(str->c[n+56]+dr*str->c[n+57]+dr2*str->c[n+58]+dr3*str->c[n+59])
		 +3*dz2*(str->c[n+60]+dr*str->c[n+61]+dr2*str->c[n+62]+dr3*str->c[n+63])));
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
void interp3Dexpl_free(interp3D_data* str) {
    free(str->c);
}
