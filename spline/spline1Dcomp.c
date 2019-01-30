/**
 * @file spline1Dcomp.c
 * @brief Cubic spline interpolation of a 1D data set, compact form
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "spline1Dcomp.h"

/**
 * @brief Calculate compact tricubic spline interpolation coefficients in 1D
 *
 * This function calculates the compact cubic interpolation coefficients for a
 * 1D data set using one of two possible boundary conditions. Function returns
 * a pointer to the coefficient array.
 *
 * @param f 1D data to be interpolated
 * @param n number of data points
 * @param bc boundary condition flag (0) natural, (1) periodic
 * @param c array for coefficient storage with length 2*n
 */
void spline1Dcomp(real* f, int n, int bc, real* c) {

    /* Array for RHS of matrix equation         */
    real* Y = malloc(n*sizeof(real));
    /* Array superdiagonal values               */
    real* p = malloc((n-1)*sizeof(real));
    /* Array for 2nd derivative vector to solve */
    real* D = malloc(n*sizeof(real));

    if(bc==0) {
        /** NATURAL (second derivative is zero) **/

        /* Initialize RHS of equation */
        Y[0] = 0.0;
        for(int i=1; i<n-1; i++) {
            Y[i] = 6 * (f[i+1] - 2 * f[i] + f[i-1]);
        }
        Y[n-1] = 0.0;

        /* Thomas algorithm is used*/

        /* Forward sweep */
        p[0] = 0.0;
        Y[0] = Y[0]/2;
        for(int i=1; i<n-1; i++) {
            p[i] =               1 / (4 - p[i-1]);
            Y[i] = (Y[i] - Y[i-1]) / (4 - p[i-1]);
        }
        Y[n-1] = 0.0;

        /* Back substitution */
        D[n-1] = Y[n-1];
        for(int i=n-2; i>-1; i--) {
            D[i] = Y[i] - p[i] * D[i+1];
        }
    }
    else if(bc==1) {
        /** PERIODIC **/

        /* Value that starts from lower left corner and moves right */
        real l     = 1.0;
        /* Last diagonal value        */
        real dlast = 4.0;
        /* Matrix right column values */
        real* r    = malloc((n-2)*sizeof(real));
        /* Last subdiagonal value     */
        real blast;

        /* Initialize RHS of equation */
        Y[0] = 6 * (f[1] - 2 * f[0] + f[n-1]);
        for(int i=1; i<n-1; i++) {
            Y[i] = 6 * (f[i+1] - 2 * f[i] + f[i-1]);
        }
        Y[n-1] = 6 * (f[0] - 2 * f[n-1] + f[n-2]);

        /* Simplified Gauss elimination is used (own algorithm) */

        /* Forward sweep */
        p[0] =  1.0/4;
        r[0] =  1.0/4;
        Y[0] = Y[0]/4;
        for(int i=1; i<n-2; i++) {
            dlast  = dlast  - l * r[i-1];
            Y[n-1] = Y[n-1] - l * Y[i-1];
            l      = -l * p[i-1];
            p[i]   =               1 / (4 - p[i-1]);
            r[i]   =         -r[i-1] / (4 - p[i-1]);
            Y[i]   = (Y[i] - Y[i-1]) / (4 - p[i-1]);
        }
        blast  = 1.0    - l * p[n-3];
        dlast  = dlast  - l * r[n-3];
        Y[n-1] = Y[n-1] - l * Y[n-3];

        p[n-2] =              (1 - r[n-3]) / (4 - p[n-3]);
        Y[n-2] =         (Y[n-2] - Y[n-3]) / (4 - p[n-3]);
        Y[n-1] = (Y[n-1] - blast * Y[n-2]) / (dlast - blast * p[n-2]);

        /* Back substitution */
        D[n-1] = Y[n-1];
        D[n-2] = Y[n-2] - p[n-2] * D[n-1];
        for(int i=n-3; i>-1; i--) {
            D[i] = Y[i] - p[i] * D[i+1] - r[i] * D[n-1];
        }

        free(r);
    }

    /* Store solved spline coefficients. With compact we store all the way to n,
       because evaluation requires it. This also means that the periodic bc
       clause (above) doesn't need a separate code segment for storage of
       period-closing spline coefficients, like with the explicit algorithm. */
    for(int i=0; i<n; i++) {
        c[i*2]   = f[i];
        c[i*2+1] = D[i];
    }

    /* Free allocated memory */
    free(Y);
    free(p);
    free(D);
}
