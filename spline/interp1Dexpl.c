/**
 * @file interp1Dexpl.c
 * @brief Cubic spline interpolation in explicit form
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "interp.h"
#include "spline.h"

/**
 * @brief Calculate cubic spline interpolation coefficients for scalar 1D data
 *
 * This function calculates the cubic spline interpolation coefficients and
 * stores them in a pre-allocated array. Explicit cofficients are calculated.
 *
 * @param c allocated array of length n_x*4 to store the coefficients
 * @param f 1D data to be interpolated
 * @param n_x number of data points in the x axis
 * @param bc_x boundary condition for the x axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 */
int interp1Dexpl_init_coeff(real* c, real* f, int n_x, int bc_x,
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

    /* Calculate cubic spline coefficients. For each grid cell i_x, there are
       four coefficients. Note how we account for normalized grid. */

    /* Cubic spline along x, using f values to get a total of four
       coefficients */
    splineexpl(f, n_x, bc_x, c);

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
void interp1Dexpl_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x, real x_min, real x_max) {

    /* Calculate grid spacing */
    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );

    /* Initialize the interp1D_data struct */
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
 * bicubic spline interpolation coefficients of the explicit form.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 *
 * @return zero on success and one if x point is outside the domain.
 */
int interp1Dexpl_eval_f(real* f, interp1D_data* str, real x) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* Index for x variable */
    int i_x = (x-str->x_min)/str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid;
    /* Helper variables */
    real dx2 = dx*dx;
    real dx3 = dx2*dx;

    int n = i_x*4; /* Index jump to cell */

    int err = 0;

    /* Check that the coordinate is within the grid. */
    if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }

    if(!err) {
        *f = str->c[n+0]+dx*str->c[n+1]+dx2*str->c[n+2]+dx3*str->c[n+3];
    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 1D and its 1st and 2nd derivatives
 *
 * This function evaluates the interpolated value of a 1D scalar field and
 * its 1st and 2nd derivatives using bicubic spline interpolation coefficients
 * of the explicit form.
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
int interp1Dexpl_eval_df(real* f_df, interp1D_data* str, real x) {

    /* Make sure periodic coordinates are within [min, max] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* Index for x variable */
    int i_x = (x-str->x_min)/str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx = (x-(str->x_min+i_x*str->x_grid))/str->x_grid;
    /* Helper variables */
    real dx2 = dx*dx;
    real dx3 = dx2*dx;
    real xgi = 1.0/str->x_grid;
    
    int n = i_x*4; /* Index jump to cell */

    int err = 0;

    /* Check that the coordinate is within the grid. */
    if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }

    if(!err) {
        /* f */
        f_df[0] = str->c[n+0]+dx*str->c[n+1]+dx2*str->c[n+2]+dx3*str->c[n+3];

        /* df/dx */
        f_df[1] = xgi*(str->c[n+1]+2*dx*str->c[n+2]+3*dx2*str->c[n+3]);

        /* d2f/dx2 */
        f_df[2] = xgi*xgi*(2*str->c[n+2]+6*dx*str->c[n+3]);
    }

    return err;
}
