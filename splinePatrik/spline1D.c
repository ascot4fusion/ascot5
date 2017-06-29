/**
 * @file spline1D.c
 * @file Cubic spline interpolation of a 1D data set
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"

/* This function calculates the interpolation coefficients
   for a 1D data set using one of two possible boundary conditions.
   Function returns pointer to coefficient array. */

void spline1D(real* f, int n, int bc, real* c) {
    /* We cast and initialize the needed variables */
    int i;
    real* Y = malloc(n*sizeof(real));
    real* p = malloc((n-1)*sizeof(real));
    real* D = malloc(n*sizeof(real));

    if(bc==0) {
	/* NATURAL (d2=0) */
	/* Initialize right side of eq */
	Y[0] = 3*(f[1]-f[0]);
	for(i=1; i<n-1; i++) {
	    Y[i] = 3*(f[i+1]-f[i-1]);
	}
	Y[n-1] = 3*(f[n-1]-f[n-2]);
	/* Thomas algorithm */
	/* Forward sweep */
	p[0] = 1.0/2;
	Y[0] = Y[0]/2;
	for(i=1; i<n-1; i++) {
	    p[i] = 1/(4-p[i-1]);
	    Y[i] = (Y[i]-Y[i-1])/(4-p[i-1]);
	}
	Y[n-1] = (Y[n-1]-Y[n-2])/(2-p[n-2]);
	/* Back substitution */
        D[n-1] = Y[n-1];
	for(i=n-2; i>-1; i--) {
	    D[i] = Y[i]-p[i]*D[i+1];
	}
    }
    else if(bc==1) {
	/* PERIODIC */
	/* Initialize some additional necessary variables */
	real l = 1.0;
	real dlast = 4.0;
	real* r = malloc((n-2)*sizeof(real));
	real blast;
	/* Initialize right side of eq */
	Y[0] = 3*(f[1] - f[n-1]);
	for(i=1; i<n-1; i++) {
	    Y[i] = 3*(f[i+1]-f[i-1]);
	}
	Y[n-1] = 3*(f[0]-f[n-2]);
	/* Forward sweep */
	p[0] = 1.0/4;
	r[0] = 1.0/4;
	Y[0] = Y[0]/4;
	for(i=1; i<n-2; i++) {
	    dlast = dlast-l*r[i-1];
	    Y[n-1] = Y[n-1]-l*Y[i-1];
	    l = -l*p[i-1];
	    p[i] = 1/(4-p[i-1]);
	    r[i] = (-r[i-1])/(4-p[i-1]);
	    Y[i] = (Y[i]-Y[i-1])/(4-p[i-1]);
	}
	blast = 1.0-l*p[n-3];
	dlast = dlast-l*r[n-3];
	Y[n-1] = Y[n-1]-l*Y[n-3];
	p[n-2] = (1-r[n-3])/(4-p[n-3]);
	Y[n-2] = (Y[n-2]-Y[n-3])/(4-p[n-3]);
	Y[n-1] = (Y[n-1]-blast*Y[n-2])/(dlast-blast*p[n-2]);
	/* Back substitution */
	D[n-1] = Y[n-1];
	D[n-2] = Y[n-2]-p[n-2]*D[n-1];
	for(i=n-3; i>-1; i--) {
	    D[i] = Y[i]-p[i]*D[i+1]-r[i]*D[n-1];
	}
	/* Free malloc */
	free(r);
	/* Period closing spline coefs */
	c[(n-1)*4] = f[n-1];
	c[(n-1)*4+1] = D[n-1];
	c[(n-1)*4+2] = 3*(f[0]-f[n-1])-2*D[n-1]-D[0];
	c[(n-1)*4+3] = 2*(f[n-1]-f[0])+D[n-1]+D[0];
    }
    /* else { */
    /* 	return -1; // How should this be done? */
    /* } */

    for(i=0; i<n-1; i++) {
	c[i*4] = f[i];
	c[i*4+1] = D[i];
	c[i*4+2] = 3*(f[i+1]-f[i])-2*D[i]-D[i+1];
        c[i*4+3] = 2*(f[i]-f[i+1])+D[i]+D[i+1];
    }

    /* Free mallocs */
    free(Y);
    free(p);
    free(D);
}
