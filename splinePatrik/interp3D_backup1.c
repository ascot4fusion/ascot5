/**
 * @file interp3D.c
 * @file Cubic spline interpolation of 3D magnetic field, tricubic
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "B_field.h"
#include "spline1D.c"

/* This function calculates the interpolation coefficients for a
   tricubic spline interpolation of a 3D magnetic field.*/
real* interp3D(real* f, int nr, int np, int nz, real* c) {
    /* We cast and initialize the needed variables */
    int ir;
    int ip;
    int iz;
    int ic;
    int nn = nz*nr;
    int nnn = nz*np*nr;
    real* ft = malloc(nr*sizeof(real)); // Best method? Maybe 3 custom length arrays.(vv)
    real* ct = malloc(4*nr*sizeof(real));
    int is;
    int iss;
    int ict;
    
    /* Bicubic spline surface over z and r for each phi */
    for(ip=0; ip<np; ip++) {
	/* Cubic spline along r for each z */
	for(iz=0; iz<nz; iz++) {
	    ft = realloc(ft, nr*sizeof(real));
	    for(ir=0; ir<nr; ir++) {
		ft[ir] = f[ip*nz*nr+iz*nr+ir];
	    }
	    ct = realloc(ct, 4*nr*sizeof(real)); // How many slots needed?
	    spline1D(ft,nr,0,ct);
	    for(ic=0; ic<4; ic++) {
		ict = ic;
		for(ir=0; ir<nr-1; ir++) { // Number of indies?
		    c[ic*nnn+ip*nn+iz*nr+ir] = ct[ict*nr+ir];
		}
	    }
	}
	
	/* Four cubic splines along z for each r using four different data sets */
	for(ir=0; ir<nr-1; ir++) {
	    /* s0, s1, s2, s3 */
	    for(is=0; is<4; is++) {
		ict = 0;
		ft = realloc(ft, nz*sizeof(real));
		for(iz=0; iz<nz; iz++) {
		    ft[iz] = c[is*nnn+ip*nn+iz*nr+ir];
		}
		ct = realloc(ct, 4*nz*sizeof(real)); // How many slots needed?
		spline1D(ft,nz,0,ct);
		for(ic=is; ic<16; ic=ic+4) {
		    for(iz=0; iz<nz-1; iz++) { // Number of indies?
			c[ic*nnn+ip*nn+iz*nr+ir] = ct[ict*nz+iz];
		    }
		    ict++;
		}
	    }
	}
    }

    /* Cubic spline along phi for each pair of z and r to find the coefficients
       of the tricubic spline volume */
    for(iz=0; iz<nz-1; iz++) {
	for(ir=0; ir<nr-1; ir++) {
	    for(iss=0; iss<4; iss++) {
		for(is=0; is<4; is++) {
		    ict = 0;
		    ft = realloc(ft, np*sizeof(real));
		    for(ip=0; ip<np; ip++) {
			ft[ip] = c[(4*iss+is)*nnn+ip*nn+iz*nr+ir];
		    }
		    ct = realloc(ct, 4*np*sizeof(real));
		    spline1D(ft,np,1,ct);
		    for(ic=4*iss+is; ic<64; ic=ic+16) {
			for(ip=0; ip<np; ip++) {
			    c[ic*nnn+ip*nn+iz*nr+ir] = ct[ict*np+ip];
			}
			ict++;
		    }
		}
	    }
	}
    }

    /* We free allocated memory */
    free(ft);
    free(ct);
}
