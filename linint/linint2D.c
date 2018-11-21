/**
 * @file linint2D.c
 * @brief Bilinear interpolation
 */
#include <stdlib.h>
#include <string.h>         /* For memcpy */
#include <stdio.h> /* Needed for printf debugging purposes */
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "linint2D.h"

/**
 * @brief Initialize linear interpolation struct for scalar 2D data
 *
 * @param str data struct for data interpolation
 * @param f 2D data to be interpolated
 * @param n_r number of data points in the r direction
 * @param n_z number of data points in the z direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 * @param r_grid grid size of the r axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 * @param z_grid grid size of the z axis
 */
int linint2D_init(linint2D_data* str, real* f, int n_r, int n_z,
                  real r_min, real r_max, real r_grid,
                  real z_min, real z_max, real z_grid) {

    int err = 0;

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;
    str->c = malloc(n_r*n_z*sizeof(real));
    for(int i = 0; i < n_r*n_z; i++) {
        str->c[i] = f[i];
    }

    return err;
}

/**
 * @brief Evaluate interpolated value of 2D scalar field
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bilinear interpolation.
 *
 * @param val variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param z z-coordinate
 */
integer linint2D_eval(real* val, linint2D_data* str, real r, real z) {
    real c00, c01, c10, c11;
    real c0, c1;

    int i_r = (r - str->r_min)/str->r_grid;     /**< index for r variable */
    int i_z = (z - str->z_min)/str->z_grid;     /**< index for z variable */

    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
                                                               current cell */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
                                                               current cell */

    int z1 = str->n_r;                      /**< Index jump one z forward */

    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
       || z < str->z_min || z > str->z_max) {
        err = 1;
    }
    else {
        /* Values at grid cell corners */
        c00 = str->c[i_z*z1 + i_r];
        c10 = str->c[i_z*z1 + (i_r + 1)];
        c01 = str->c[(i_z + 1)*z1 + i_r];
        c11 = str->c[(i_z + 1)*z1 + (i_r + 1)];
        /* Interpolate along r */
        c0 = c00*(1 - dr) + c10*dr;
        c1 = c01*(1 - dr) + c11*dr;
        /* Finally we interpolate these values along z */
        val[0] = c0*(1 - dz) + c1*dz;
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
void linint2D_free(linint2D_data* str) {
    free(str->c);
}
