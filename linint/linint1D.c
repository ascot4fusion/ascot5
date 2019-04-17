/**
 * @file linint1D.c
 * @brief Linear interpolation
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../spline/interp.h"
#include "linint.h"

/**
 * @brief Initialize linear interpolation struct for scalar 1D data
 *
 * @param str pointer to struct to be initialized
 * @param c array where data is stored
 * @param n_x number of data points in the x direction
 * @param bc_x boundary condition for x axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 */
void linint1D_init(linint1D_data* str, real* c,
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
 * linear interpolation.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 *
 * @return zero on success and one if x point is outside the grid.
 */
int linint1D_eval_f(real* f, linint1D_data* str, real x) {

    /* Make sure periodic coordinates are within [max, min] region. */
    if(str->bc_x == PERIODICBC) {
        x = fmod(x - str->x_min, str->x_max - str->x_min) + str->x_min;
        x = x + (x < str->x_min) * (str->x_max - str->x_min);
    }

    /* index for x variable */
    int i_x   = (x-str->x_min) / str->x_grid;
    /**< Normalized x coordinate in current cell */
    real dx = ( x - (str->x_min + i_x*str->x_grid)) / str->x_grid;

    int x1 = 1;   /* Index jump one r forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && (x < str->x_min || x > str->x_max) ) {
        err = 1;
    }

    if(!err) {
        *f = str->c[i_x]*(1 - dx) + str->c[i_x+x1]*dx;
    }
    return err;
}
