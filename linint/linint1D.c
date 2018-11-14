/**
 * @file linint1D.c
 * @brief Linear interpolation
 */
#include <stdlib.h>
#include <string.h>         /* For memcpy */
#include <stdio.h> /* Needed for printf debugging purposes */
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "linint1D.h"

/**
 * @brief Initialize linear interpolation struct for scalar 1D data
 *
 * @param str data struct for data interpolation
 * @param f data to be interpolated
 * @param n_r number of data points in the r direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 * @param r_grid grid size of the r axis
 */
int linint1D_init(linint1D_data* str, real* f, int n_r,
                  real r_min, real r_max, real r_grid,
                  int periodic) {

    int err = 0;

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    if(periodic) {
        str->r_max += str->r_grid;    /**< We can evaluate outside the last point as well */
    }
    str->f = f;

    return err;
}

/**
 * @brief Evaluate interpolated value of 1D scalar field
 *
 * This function evaluates the interpolated value of a 1D scalar field using
 * linear interpolation.
 *
 * @param val variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 */
integer linint1D_eval(real* val, linint1D_data* str, real r) {
    real c0, c1;
    int i_r = (r - str->r_min)/str->r_grid;     /**< index for r variable */
    int r1 = 1;                                 /**< Index jump one r forward */
    if(i_r == str->n_r-1) {
      r1 = -(str->n_r-1)*r1;                    /**< If last cell, index jump to 1st r */
    }
    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
                                                               current cell */
    int err = 0;
    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max) {
        err = 1;
    }
    else {
        c0 = str->f[i_r];
        c1 = str->f[i_r+r1];
        val[0] = c0*(1 - dr) + c1*dr;
    }

    return err;
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
void linint1D_free(linint1D_data* str) {
    free(str->f);
}
