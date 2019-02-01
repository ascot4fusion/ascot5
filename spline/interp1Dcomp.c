/**
 * @file interp1Dcomp.c
 * @brief Cubic spline interpolation in compact form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "interp.h"
#include "spline1Dcomp.h"

/**
 * @brief Calculate cubic spline interpolation coefficients for scalar 1D data
 *
 * This function calculates the cubic spline interpolation coefficients for
 * the given data and stores them in the data struct. Compact  cofficients are
 * calculated directly.
 *
 * @param c allocated array of length n_x*2 to store the coefficients
 * @param f 1D data to be interpolated
 * @param n_x number of data points in the x axis
 * @param bc_x boundary condition for the x axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 */
int interp1Dcomp_init_coeff(real* c, real* f, int n_x, int bc_x,
                            real x_min, real x_max) {

    /* Check boundary conditions and evaluate grid interval */
    real x_grid;
    if(bc_x == PERIODICBC || bc_x == NATURALBC) {
        x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    }
    else {
        return 1;
    }

    if(c == NULL) {
        return 1;
    }
    spline1Dcomp(f, n_x, bc_x, c);
    for(int i_x=0; i_x<n_x; i_x++) {
        c[i_x*2]     = c[i_x*2];
        c[i_x*2 + 1] = c[i_x*2+1] / (x_grid*x_grid);
    }

    return 0;
}

/**
 * @brief Initialize a cubic spline
 *
 * @param str pointer to spline to be initialized
 * @param c array where coefficients are stored
 * @param n_x number of data points in the x direction
 * @param bc_x boundary condition for x axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 */
void interp1Dcomp_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x, real x_min, real x_max) {

    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );

    str->n_x    = n_x;
    str->bc_x   = bc_x;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of 1D scalar field
 *
 * This function evaluates the interpolated value of a 1D scalar field using
 * bicubic spline interpolation coefficients of the compact form.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 *
 * @return zero on success and one if x point is outside the grid.
 */
int interp1Dcomp_eval_f(real* f, interp1D_data* str, real x) {

    /* Make sure periodic coordinates are within [max, min] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* index for x variable */
    int i_x   = (x-str->x_min) / str->x_grid;
    /**< Normalized x coordinate in current cell */
    real dx   = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    real dxi  = 1.0 - dx;
    real dxi3 = dxi * (dxi*dxi - 1.0);
    real xg2  = str->x_grid*str->x_grid;

    int n  = i_x*2; /* Index jump to cell       */
    int x1 = 2;     /* Index jump one x forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }

    if(!err) {
        *f = (
              dxi*str->c[n] +
              dx*str->c[n+x1] +
              (xg2/6)*(dxi3*str->c[n+1] + dxi3*str->c[n+x1+1])
              );
    }
    return err;
}

/**
 * @brief Evaluate interpolated value of 1D and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 1D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the compact form.
 *
 * The evaluated  values are returned in an array with following elements:
 * - f_df[0] = f
 * - f_df[1] = f_x
 * - f_df[2] = f_xx
 *
 * @param f_df array in which to place the evaluated values
 * @param str data struct for data interpolation
 * @param x x-coordinate
 *
 * @return zero on success and one if (x,y) point is outside the grid.
 */
int interp1Dcomp_eval_df(real* f_df, interp1D_data* str, real x) {

    /* Make sure periodic coordinates are within [max, min] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /**< index for x variable */
    int i_x     = (x - str->x_min) / str->x_grid;
    /**< Normalized x coordinate in current cell */
    real dx     = ( x - (str->x_min + i_x*str->x_grid) ) / str->x_grid;
    real dx3dx  = 3*dx*dx - 1;
    real dxi    = 1.0 - dx;
    real dxi3   = dxi * (dxi*dxi - 1);
    real dxi3dx = -3*dxi*dxi + 1;
    real xg     = str->x_grid;
    real xg2    = xg*xg;
    real xgi    = 1.0 / xg;

    int n  = i_x*2; /* Index jump to cell */
    int x1 = 2;     /* Index jump one x forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }

    if(!err) {
        /* f */
        f_df[0] = (
                   dxi*str->c[n] +
                   dx*str->c[n+x1] +
                   (xg2/6)*(dxi3*str->c[n+1] + dxi3*str->c[n+x1+1])
                   );

        /* df/dr */
        f_df[1] = (xgi*(str->c[n+x1] - str->c[n])  +
                   (xg/6)*(dx3dx*str->c[n+x1+1] +dxi3dx*str->c[n+1])
                   );

        /* d2f/dr2 */
        f_df[2] = (
                   dxi*str->c[n+1] + dx*str->c[n+x1+1]
                   );
    }
    return err;
}
