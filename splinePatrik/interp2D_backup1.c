/**
 * @file interp2D.c
 * @file Cubic spline interpolation of 2D magnetic field component, tricubic
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "B_field.h"
#include "spline1D.c"

/* This function calculates the interpolation coefficients for a
   bicubic spline interpolation of a 2D magnetic field component.*/

real* interp2D(real* f, int nr, int nz, real* c) {
    /* We cast and initialize the needed variables */
    int ir;
    int iz;
    int ic;
    int nn = nz*nr;
    /* real* c = malloc(16*nn*sizeof(real)); // As parameter? Unused elements at edges? */
    real* ft;
    real* ct;
    int is;
    int ict;

    /* Cubic spline along r for each z */
    for(iz=0; iz<nz; iz++) {
	ft = malloc(nr*sizeof(real)); // Change this!
	for(ir=0; ir<nr; ir++) {
	    ft[ir] = f[iz*nr+ir];
	}
	ct = malloc(4*nr*sizeof(real));
	spline1D(ft,nr,0,ct);
	for(ic=0; ic<4; ic++) {
	    ict = ic;
	    for(ir=0; ir<nr; ir++) {
		c[ic*nn+iz*nr+ir] = ct[ict*nr+ir];
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
		ft[iz] = c[is*nn+iz*nr+ir];
	    }
	    ct = realloc(ct, 4*nz*sizeof(real));
	    spline1D(ft,nz,0,ct);
	    for(ic=is; ic<16; ic=ic+4) {
		for(iz=0; iz<nz; iz++) {
		    c[ic*nn+iz*nr+ir] = ct[ict*nz+iz];
		}
		ict++;
	    }
	}
    }

    /* We free allocated memory */
    free(ft);
    free(ct);

    /* /\* We return the coefficient *\/ */
    /* return c; // Probably c needs to be given as parameter to avoid malloc problem? */
}
