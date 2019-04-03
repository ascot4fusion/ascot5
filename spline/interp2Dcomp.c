/**
 * @file interp2Dcomp.c
 * @brief Bicubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "interp.h"
#include "spline.h"

/**
 * @brief Calculate bicubic spline interpolation coefficients for scalar 2D data
 *
 * This function calculates the bicubic spline interpolation coefficients and
 * stores them in a pre-allocated array. Compact cofficients are calculated.
 *
 * For each data point four coefficients are stored for spline-interpolation.
 *
 * @param c allocated array of length n_y*n_x*4 to store the coefficients
 * @param f 2D data to be interpolated
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis (0) natural (1) periodic
 * @param bc_y boundary condition for y axis (0) natural (1) periodic
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 *
 * @return zero if initialization succeeded
 */
int interp2Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
			    real y_min, real y_max) {

    /* Check boundary conditions and evaluate grid interval */
    real x_grid, y_grid;
    if(bc_x == PERIODICBC || bc_x == NATURALBC) {
        x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    }
    else {
        return 1;
    }

    if(bc_y == PERIODICBC || bc_y == NATURALBC) {
        y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );
    }
    else {
        return 1;
    }

    /* Allocate helper quantities */
    real* f_x = malloc(n_x*sizeof(real));
    real* f_y = malloc(n_y*sizeof(real));
    real* c_x = malloc(n_x*2*sizeof(real));
    real* c_y = malloc(n_y*2*sizeof(real));

    if(f_x == NULL || f_y == NULL || c_x == NULL || c_y == NULL) {
        return 1;
    }

    /* Calculate bicubic spline surface coefficients, i.e second derivatives.
       For each grid cell (i_x, i_y), there are four coefficients:
       [f, fxx, fyy, fxxyy]. Note how we account for normalized grid. */

    /* Cubic spline along x for each y, using f values to get fxx */
    for(int i_y=0; i_y<n_y; i_y++) {
	/* fxx */
        for(int i_x=0; i_x<n_x; i_x++) {
            f_x[i_x] = f[i_y*n_x+i_x];
        }
        splinecomp(f_x, n_x, bc_x, c_x);
        for(int i_x=0; i_x<n_x; i_x++) {
            c[i_y*n_x*4 + i_x*4    ] = c_x[i_x*2];
            c[i_y*n_x*4 + i_x*4 + 1] = c_x[i_x*2+1] / (x_grid*x_grid);
        }
    }

    /* Two cubic splines along y for each x, using f and fxx to get fyy and
       fxxyy */
    for(int i_x=0; i_x<n_x; i_x++) {

        /* fyy */
        for(int i_y=0; i_y<n_y; i_y++) {
            f_y[i_y] =  f[i_y*n_x + i_x];
        }
        splinecomp(f_y, n_y, bc_y, c_y);
        for(int i_y=0; i_y<n_y; i_y++) {
            c[i_y*n_x*4+i_x*4+2] = c_y[i_y*2+1]/(y_grid*y_grid);
        }

        /* fxxyy */
        for(int i_y=0; i_y<n_y; i_y++) {
            f_y[i_y] =  c[i_y*n_x*4 + i_x*4 + 1];
        }
        splinecomp(f_y, n_y, bc_y, c_y);
        for(int i_y=0; i_y<n_y; i_y++) {
            c[i_y*n_x*4 + i_x*4 + 3] = c_y[i_y*2 + 1] / (y_grid*y_grid);
        }
    }

    /* Free allocated memory */
    free(f_x);
    free(f_y);
    free(c_x);
    free(c_y);

    return 0;
}

/**
 * @brief Initialize a bicubic spline
 *
 * @param str pointer to spline to be initialized
 * @param c array where coefficients are stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 */
void interp2Dcomp_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
			      real y_min, real y_max) {

    /* Calculate grid spacings */
    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );

    /* Initialize the interp2D_data struct */
    str->n_x    = n_x;
    str->n_y    = n_y;
    str->bc_x   = bc_x;
    str->bc_y   = bc_y;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of a 2D field
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bicubic spline interpolation coefficients of the compact form.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 *
 * @return zero on success and one if (x,y) point is outside the domain.
 */
int interp2Dcomp_eval_f(real* f, interp2D_data* str, real x, real y) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable */
    int i_x   = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx   = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    /* Helper variables */
    real dx3  =  dx * (dx*dx - 1.0);
    real dxi  = 1.0 - dx;
    real dxi3 = dxi * (dxi*dxi - 1.0);
    real xg2  = str->x_grid*str->x_grid;

    /* Index for y variable */
    int i_y   = (y - str->y_min) / str->y_grid;
    /* Normalized y coordinate in current cell */
    real dy   = ( y - (str->y_min + i_y*str->y_grid) ) / str->y_grid;
    /* Helper variables */
    real dy3  =  dy * (dy*dy - 1.0);
    real dyi  = 1.0 - dy;
    real dyi3 = dyi * (dyi*dyi - 1.0);
    real yg2  = str->y_grid*str->y_grid;

    int n  = i_y*str->n_x*4+i_x*4; /* Index jump to cell       */
    int x1 = 4;                    /* Index jump one x forward */
    int y1 = str->n_x*4;           /* Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the domain. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }

    if(!err) {
        *f = (
            dxi*(dyi*str->c[n]+dy*str->c[n+y1])+
            dx*(dyi*str->c[n+x1]+dy*str->c[n+y1+x1]))
            +(xg2/6)*(
                dxi3*(dyi*str->c[n+1] + dy*str->c[n+y1+1])+
                dx3*(dyi*str->c[n+x1+1] + dy*str->c[n+y1+x1+1]))
            +(yg2/6)*(
                dxi*(dyi3*str->c[n+2]+dy3*str->c[n+y1+2])+
                dx*(dyi3*str->c[n+x1+2]+dy3*str->c[n+y1+x1+2]))
            +(xg2*yg2/36)*(
                dxi3*(dyi3*str->c[n+3]+dy3*str->c[n+y1+3])+
                dx3*(dyi3*str->c[n+x1+3]+dy3*str->c[n+y1+x1+3]));
    }

    return err;
}

/**
 * @brief Evaluate interpolated value and 1st and 2nd derivatives of 2D field
 *
 * This function evaluates the interpolated value of a 2D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 *
 * The evaluated  values are returned in an array with following elements:
 * - f_df[0] = f
 * - f_df[1] = f_x
 * - f_df[2] = f_y
 * - f_df[3] = f_xx
 * - f_df[4] = f_yy
 * - f_df[5] = f_xy
 *
 * @param f_df array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 *
 * @return zero on success and one if (x,y) point is outside the grid.
 */
int interp2Dcomp_eval_df(real* f_df, interp2D_data* str, real x, real y) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }
    if(str->bc_y == PERIODICBC) {
        y = fmod(y - str->y_min, str->y_max - str->y_min) + str->y_min;
        y = y + (y < str->y_min) * (str->y_max - str->y_min);
    }

    /* Index for x variable */
    int i_x   = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx     = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    /* Helper variables */
    real dx3    =  dx * (dx*dx - 1.0);
    real dx3dx  = 3*dx*dx - 1;
    real dxi    = 1.0 - dx;
    real dxi3   = dxi * (dxi*dxi - 1.0);
    real dxi3dx = -3*dxi*dxi + 1;
    real xg     = str->x_grid;
    real xg2    = xg*xg;
    real xgi    = 1.0/xg;

    /* Index for y variable */
    int i_y   = (y - str->y_min) / str->y_grid;
    /* Normalized y coordinate in current cell */
    real dy     = ( y - (str->y_min + i_y*str->y_grid) ) / str->y_grid;
    /* Helper variables */
    real dy3    =  dy * (dy*dy - 1.0);
    real dy3dy  = 3*dy*dy - 1;
    real dyi    = 1.0 - dy;
    real dyi3   = dyi * (dyi*dyi - 1.0);
    real dyi3dy = -3*dyi*dyi + 1;
    real yg     = str->y_grid;
    real yg2    = yg*yg;
    real ygi    = 1.0/yg;

    int n  = i_y*str->n_x*4+i_x*4; /**< Index jump to cell       */
    int x1 = 4;                    /**< Index jump one x forward */
    int y1 = str->n_x*4;           /**< Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the domain. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && (y < str->y_min || y > str->y_max) ) {
        err = 1;
    }

    if(!err) {
        /* f */
        f_df[0] = (
            dxi*(dyi*str->c[n]+dy*str->c[n+y1])+
            dx*(dyi*str->c[n+x1]+dy*str->c[n+y1+x1]))
            +(xg2/6)*(
                dxi3*(dyi*str->c[n+1] + dy*str->c[n+y1+1])+
                dx3*(dyi*str->c[n+x1+1] + dy*str->c[n+y1+x1+1]))
            +(yg2/6)*(
                dxi*(dyi3*str->c[n+2]+dy3*str->c[n+y1+2])+
                dx*(dyi3*str->c[n+x1+2]+dy3*str->c[n+y1+x1+2]))
            +(xg2*yg2/36)*(
                dxi3*(dyi3*str->c[n+3]+dy3*str->c[n+y1+3])+
                dx3*(dyi3*str->c[n+x1+3]+dy3*str->c[n+y1+x1+3]));

        /* df/dx */
        f_df[1] = xgi*(
            -(dyi*str->c[n]  +dy*str->c[n+y1])
            +(dyi*str->c[n+x1]+dy*str->c[n+y1+x1]))
            +(xg/6)*(
                dxi3dx*(dyi*str->c[n+1]  +dy*str->c[n+y1+1])+
                dx3dx*(dyi*str->c[n+x1+1]+dy*str->c[n+y1+x1+1]))
            +(xgi*yg2/6)*(
                -(dyi3*str->c[n+2]  +dy3*str->c[n+y1+2])
                +(dyi3*str->c[n+x1+2]+dy3*str->c[n+y1+x1+2]))
            +(xg*yg2/36)*(
                dxi3dx*(dyi3*str->c[n+3]  +dy3*str->c[n+y1+3])+
                dx3dx*(dyi3*str->c[n+x1+3]+dy3*str->c[n+y1+x1+3]));

        /* df/dy */
        f_df[2] = ygi*(
            dxi*(-str->c[n]  +str->c[n+y1])+
            dx*(-str->c[n+x1]+str->c[n+y1+x1]))
            +(xg2*ygi/6)*(
                dxi3*(-str->c[n+1]  +str->c[n+y1+1])+
                dx3*(-str->c[n+x1+1]+str->c[n+y1+x1+1]))
            +(yg/6)*(
                dxi*(dyi3dy*str->c[n+2]  +dy3dy*str->c[n+y1+2])+
                dx*(dyi3dy*str->c[n+x1+2]+dy3dy*str->c[n+y1+x1+2]))
            +(xg2*yg/36)*(
                dxi3*(dyi3dy*str->c[n+3]  +dy3dy*str->c[n+y1+3])+
                dx3*(dyi3dy*str->c[n+x1+3]+dy3dy*str->c[n+y1+x1+3]));

        /* d2f/dx2 */
        f_df[3] = (
            dxi*(dyi*str->c[n+1]  +dy*str->c[n+y1+1])+
            dx*(dyi*str->c[n+x1+1]+dy*str->c[n+y1+x1+1]))
            +(yg2/6)*(
                dxi*(dyi3*str->c[n+3]  +dy3*str->c[n+y1+3])+
                dx*(dyi3*str->c[n+x1+3]+dy3*str->c[n+y1+x1+3]));

        /* d2f/dy2 */
        f_df[4] = (
              dxi*(dyi*str->c[n+2]  +dy*str->c[n+y1+2])+
              dx*(dyi*str->c[n+x1+2]+dy*str->c[n+y1+x1+2]))
        +xg2/6*(
            dxi3*(dyi*str->c[n+3]  +dy*str->c[n+y1+3])+
            dx3*(dyi*str->c[n+x1+3]+dy*str->c[n+y1+x1+3]));

        /* d2f/dydx */
        f_df[5] = xgi*ygi*(
            str->c[n]  -str->c[n+y1]
            -str->c[n+x1]+str->c[n+y1+x1])
            +(xg/6*ygi)*(
                dxi3dx*(-str->c[n+1]  +str->c[n+y1+1])+
                dx3dx*(-str->c[n+x1+1]+str->c[n+y1+x1+1]))
            +(xgi/6*yg)*(
                -(dyi3dy*str->c[n+2]  +dy3dy*str->c[n+y1+2])
                +(dyi3dy*str->c[n+x1+2]+dy3dy*str->c[n+y1+x1+2]))
            +(xg*yg/36)*(
                dxi3dx*(dyi3dy*str->c[n+3]  +dy3dy*str->c[n+y1+3])+
                dx3dx*(dyi3dy*str->c[n+x1+3]+dy3dy*str->c[n+y1+x1+3]));
    }
    
    return err;
}
