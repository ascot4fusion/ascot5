/**
 * @file linint3D.c
 * @brief Trilinear interpolation
 */
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "linint.h"

/**
 * @brief Initialize linear interpolation struct for scalar 3D data
 *
 * @param str pointer to struct to be initialized
 * @param c array where data is stored
 * @param n_x number of data points in the x direction
 * @param n_y number of data points in the y direction
 * @param n_z number of data points in the z direction
 * @param bc_x boundary condition for x axis
 * @param bc_y boundary condition for y axis
 * @param bc_y boundary condition for z axis
 * @param x_min minimum value of the x axis
 * @param x_max maximum value of the x axis
 * @param y_min minimum value of the y axis
 * @param y_max maximum value of the y axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 */
void linint3D_init(linint3D_data* str, real* c,
                   int n_x, int n_y, int n_z,
                   int bc_x, int bc_y, int bc_z,
                   real x_min, real x_max,
                   real y_min, real y_max,
                   real z_min, real z_max) {

    real x_grid = (x_max - x_min) / ( n_x - 1 * (bc_x == NATURALBC) );
    real y_grid = (y_max - y_min) / ( n_y - 1 * (bc_y == NATURALBC) );
    real z_grid = (z_max - z_min) / ( n_z - 1 * (bc_z == NATURALBC) );

    str->n_x    = n_x;
    str->n_y    = n_y;
    str->n_z    = n_z;
    str->bc_x   = bc_x;
    str->bc_y   = bc_y;
    str->bc_z   = bc_z;
    str->x_min  = x_min;
    str->x_max  = x_max;
    str->x_grid = x_grid;
    str->y_min  = y_min;
    str->y_max  = y_max;
    str->y_grid = y_grid;
    str->z_min  = z_min;
    str->z_max  = z_max;
    str->z_grid = z_grid;
    str->c      = c;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * trilinear interpolation.
 *
 * @param f variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 */
int linint3D_eval_f(real* f, linint3D_data* str, real x, real y, real z) {
    real c000, c100, c001, c101, c010, c110, c011, c111;
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
    if(str->bc_z == PERIODICBC) {
        z = fmod(z - str->z_min, str->z_max - str->z_min) + str->z_min;
        z = z + (z < str->z_min) * (str->z_max - str->z_min);
    }

    /* index for x variable */
    int i_x   = (x - str->x_min) / str->x_grid;
    /* Normalized x coordinate in current cell */
    real dx   = (x - (str->x_min + i_x*str->x_grid)) / str->x_grid;

    /* index for y variable */
    int i_y   = (y - str->y_min) / str->y_grid;
    /* Normalized y coordinate in current cell */
    real dy   = (y - (str->y_min + i_y*str->y_grid)) / str->y_grid;

    /* index for z variable */
    int i_z   = (z - str->z_min) / str->z_grid;
    /* Normalized z coordinate in current cell */
    real dz   = (z - (str->z_min + i_z*str->z_grid)) / str->z_grid;

    int x1 = 1;                 /* Index jump one x forward */
    int y1 = str->n_z*str->n_x; /* Index jump one y forward */
    int z1 = str->n_x;          /* Index jump one z forward */

    int err = 0;

    /* Enforce periodic BC or check that the coordinate is within the grid. */
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
    if( str->bc_z == PERIODICBC && i_z == str->n_z-1 ) {
        z1 = -(str->n_z-1)*z1;
    }
    else if( str->bc_z == NATURALBC && (z < str->z_min || z > str->z_max) ) {
        err = 1;
    }

    if(!err) {
        /* Values at grid cell corners */
        c000 = str->c[i_y*y1 + i_z*z1 + i_x];
        c100 = str->c[i_y*y1 + i_z*z1 + (i_x + 1)];
        c001 = str->c[(i_y + 1)*y1 + i_z*z1 + i_x];
        c101 = str->c[(i_y + 1)*y1 + i_z*z1 + (i_x + 1)];
        c010 = str->c[i_y*y1 + (i_z + 1)*z1 + i_x];
        c110 = str->c[i_y*y1 + (i_z + 1)*z1 + (i_x + 1)];
        c011 = str->c[(i_y + 1)*y1 + (i_z + 1)*z1 + i_x];
        c111 = str->c[(i_y + 1)*y1 + (i_z + 1)*z1 + (i_x + 1)];
        /* Interpolate along x */
        c00 = c000*(1 - dx) + c100*dx;
        c01 = c001*(1 - dx) + c101*dx;
        c10 = c010*(1 - dx) + c110*dx;
        c11 = c011*(1 - dx) + c111*dx;
        /* Interpolate these values along z */
        c0 = c00*(1 - dz) + c10*dz;
        c1 = c01*(1 - dz) + c11*dz;
        /* Finally we interpolate these values along y */
        *f = c0*(1 - dy) + c1*dy;
    }

    return err;
}
