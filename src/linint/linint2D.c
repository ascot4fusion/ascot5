/**
 * @file linint2D.c
 * @brief Bilinear interpolation
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../math.h"
#include "linint.h"

/**
 * @brief Initialize linear interpolation struct for scalar 2D data
 *
 * @param str pointer to struct to be initialized
 * @param c array where data is stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 */
void linint2D_init(linint2D_data* str, real* c,
                   int n_x, int n_y,
                   int bc_x, int bc_y,
                   real x_min, real x_max,
                   real y_min, real y_max) {

    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );

    str->n_x    = n_x;
    str->n_y    = n_y;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of 2D scalar field
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bilinear interpolation.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 *
 * @return zero on success and one if (x,y) point is outside the grid.
 */
int linint2D_eval_f(real* f, linint2D_data* str, real x, real y) {
    real c00, c01, c10, c11;
    real c0, c1;

    /* Make sure periodic coordinates are within [max, min] region. */
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

    /* Index for y variable */
    int i_y   = (y - str->y_min) / str->y_grid;
    /* Normalized y coordinate in current cell */
    real dy   = ( y - (str->y_min + i_y*str->y_grid) ) / str->y_grid;

    int n  = i_y*str->n_x + i_x; /* Index jump to cell       */
    int x1 = 1;                  /* Index jump one x forward */
    int y1 = str->n_x;           /* Index jump one y forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
    if( str->bc_x == PERIODICBC && i_x == str->n_x-1 ) {
        x1 = -(str->n_x-1)*x1;
    }
    else if( str->bc_x == NATURALBC && !(x >= str->x_min && x <= str->x_max) ) {
        err = 1;
    }
    if( str->bc_y == PERIODICBC && i_y == str->n_y-1 ) {
        y1 = -(str->n_y-1)*y1;
    }
    else if( str->bc_y == NATURALBC && !(y >= str->y_min && y <= str->y_max) ) {
        err = 1;
    }

    if(!err) {
        /* Values at grid cell corners */
        c00 = str->c[n];
        c10 = str->c[n + x1];
        c01 = str->c[n + y1];
        c11 = str->c[n + y1 + x1];
        /* Interpolate along x */
        c0 = c00*(1 - dx) + c10*dx;
        c1 = c01*(1 - dx) + c11*dx;
        /* Finally we interpolate these values along y */
        *f = c0*(1 - dy) + c1*dy;
    }

    return err;
}
