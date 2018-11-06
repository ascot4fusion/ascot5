/**
 * @file linint3D.c
 * @brief Trilinear interpolation
 */
#include <stdlib.h>
#include <string.h>         /* For memcpy */
#include <stdio.h> /* Needed for printf debugging purposes */
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "linint3D.h"

/**
 * @brief Initialize linear interpolation struct for scalar 3D data
 *
 * @param str data struct for data interpolation
 * @param f 3D data to be interpolated
 * @param n_r number of data points in the r direction
 * @param n_phi number of data points in the phi direction
 * @param n_z number of data points in the z direction
 * @param r_min minimum value of the r axis
 * @param r_max maximum value of the r axis
 * @param r_grid grid size of the r axis
 * @param phi_min minimum value of the phi axis
 * @param phi_max maximum value of the phi axis
 * @param phi_grid grid size of the phi axis
 * @param z_min minimum value of the z axis
 * @param z_max maximum value of the z axis
 * @param z_grid grid size of the z axis
 */
int linint3D_init(linint3D_data* str, real* f, int n_r, int n_phi, int n_z,
		       real r_min, real r_max, real r_grid,
		       real phi_min, real phi_max, real phi_grid,
		       real z_min, real z_max, real z_grid) {

    int err = 0;

    /* Initialize and fill the data struct */
    str->n_r = n_r;
    str->n_phi = n_phi;
    str->n_z = n_z;
    str->r_min = r_min;
    str->r_max = r_max;
    str->r_grid = r_grid;
    str->phi_min = phi_min;
    str->phi_max = phi_max;
    str->phi_grid = phi_grid;
    str->z_min = z_min;
    str->z_max = z_max;
    str->z_grid = z_grid;
    str->f = f;

    return err;
}

/**
 * @brief Evaluate interpolated value of 3D scalar field
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * trilinear interpolation.
 * 
 * @param val variable in which to place the evaluated value
 * @param str data struct for data interpolation
 * @param r r-coordinate
 * @param phi phi-coordinate
 * @param z z-coordinate
 */
integer linint3D_eval(real* val, linint3D_data* str, real r, real phi, real z) {
    real c000, c100, c001, c101, c010, c110, c011, c111;
    real c00, c01, c10, c11;
    real c0, c1;
    /** Make sure phi is in interval [0,2pi) */
    phi = fmod(phi,CONST_2PI);
    if(phi < 0){phi = CONST_2PI + phi;}

    int i_r = (r - str->r_min)/str->r_grid;     /**< index for r variable */
    int i_phi = (phi - str->phi_min)/str->phi_grid; /**< index for phi variable */
    int i_z = (z - str->z_min)/str->z_grid;     /**< index for z variable */

    real dr = (r-(str->r_min+i_r*str->r_grid))/str->r_grid; /**< Normalized r coordinate in
							       current cell */
    real dphi = (phi-(str->phi_min+i_phi*str->phi_grid))/str->phi_grid; /**< Normalized phi
									   coordinate in
									   current cell */
    real dz = (z-(str->z_min+i_z*str->z_grid))/str->z_grid; /**< Normalized z coordinate in
							       current cell */

    int phi1 = str->n_z*str->n_r;           /**< Index jump one phi forward */
    if(i_phi==str->n_phi-1) {
	phi1 = -(str->n_phi-1)*phi1;          /**< If last cell, index jump to 1st phi */
    }
    int z1 = str->n_r;                      /**< Index jump one z forward */
	   
    int err = 0;

    /* Check that the point is not outside the evaluation regime */
    if(r < str->r_min || r > str->r_max
	|| z < str->z_min || z > str->z_max) {
	err = 1;
    }
    else {
        /* Values at grid cell corners */
        c000 = str->f[i_phi*phi1 + i_z*z1 + i_r];
        c100 = str->f[i_phi*phi1 + i_z*z1 + (i_r + 1)];
        c001 = str->f[(i_phi + 1)*phi1 + i_z*z1 + i_r];
        c101 = str->f[(i_phi + 1)*phi1 + i_z*z1 + (i_r + 1)];
        c010 = str->f[i_phi*phi1 + (i_z + 1)*z1 + i_r];
        c110 = str->f[i_phi*phi1 + (i_z + 1)*z1 + (i_r + 1)];
        c011 = str->f[(i_phi + 1)*phi1 + (i_z + 1)*z1 + i_r];
        c111 = str->f[(i_phi + 1)*phi1 + (i_z + 1)*z1 + (i_r + 1)];
        /* Interpolate along r */
        c00 = c000*(1 - dr) + c100*dr;
        c01 = c001*(1 - dr) + c101*dr;
        c10 = c010*(1 - dr) + c110*dr;
        c11 = c011*(1 - dr) + c111*dr;
        /* Interpolate these values along z */
        c0 = c00*(1 - dz) + c10*dz;
        c1 = c01*(1 - dz) + c11*dz;
        /* Finally we interpolate these values along phi */
        val[0] = c0*(1 - dphi) + c1*dphi;
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
void linint3D_free(linint3D_data* str) {
    free(str->f);
}
