/**
 * @file interp_data.h
 * Interpolation data structures.
 */
#ifndef INTERP_DATA_H
#define INTERP_DATA_H

#include "defines.h"
#include <stddef.h>

/**
 * Boundary conditions for the interpolation.
 */
enum boundaryCondition
{
    NATURALBC = 0, /**< Second derivative is zero at the boundary.            */
    PERIODICBC = 1 /**< Function has same value and derivatives at both ends. */
};

/**
 * Number of coefficients stored for each data point in splines.
 */
enum SplineSize
{
    NSIZE_COMP1D = 2,
    NSIZE_COMP2D = 4,
    NSIZE_COMP3D = 8,
};

/**
 * Cubic interpolation struct.
 */
typedef struct
{
    size_t n_x;  /**< Number of x grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real *c;     /**< Pointer to array with spline coefficients.              */
} Spline1D;

/**
 * Bicubic interpolation struct.
 */
typedef struct
{
    size_t n_x;  /**< Number of x grid points.                                */
    size_t n_y;  /**< Number of y grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    int bc_y;    /**< Boundary condition for y coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real y_min;  /**< Minimum y coordinate in the grid.                       */
    real y_max;  /**< Maximum y coordinate in the grid.                       */
    real y_grid; /**< Interval between two adjacent points in y grid.         */
    real *c;     /**< Pointer to array with spline coefficients.              */
} Spline2D;

/**
 * Tricubic interpolation struct.
 */
typedef struct
{
    size_t n_x;  /**< Number of x grid points.                                */
    size_t n_y;  /**< Number of y grid points.                                */
    size_t n_z;  /**< Number of z grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    int bc_y;    /**< Boundary condition for y coordinate.                    */
    int bc_z;    /**< Boundary condition for z coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real y_min;  /**< Minimum y coordinate in the grid.                       */
    real y_max;  /**< Maximum y coordinate in the grid.                       */
    real y_grid; /**< Interval between two adjacent points in y grid.         */
    real z_min;  /**< Minimum z coordinate in the grid.                       */
    real z_max;  /**< Maximum z coordinate in the grid.                       */
    real z_grid; /**< Interval between two adjacent points in z grid.         */
    real *c;     /**< Pointer to array with spline coefficients.              */
} Spline3D;

/**
 * 1D linear interpolation struct.
 */
typedef struct
{
    size_t n_x;  /**< Number of x grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real *c;     /**< Pointer to array with interpolant values.               */
} Linear1D;

/**
 * @brief 2D linear interpolation struct.
 */
typedef struct
{
    size_t n_x;  /**< Number of x grid points.                                */
    size_t n_y;  /**< Number of y grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    int bc_y;    /**< Boundary condition for y coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real y_min;  /**< Minimum y coordinate in the grid.                       */
    real y_max;  /**< Maximum y coordinate in the grid.                       */
    real y_grid; /**< Interval between two adjacent points in y grid.         */
    real *c;     /**< Pointer to array with interpolant values.               */
} Linear2D;

/**
 * @brief 3D linear interpolation struct.
 */
typedef struct
{
    size_t n_x;  /**< Number of x grid points.                                */
    size_t n_y;  /**< Number of y grid points.                                */
    size_t n_z;  /**< Number of z grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    int bc_y;    /**< Boundary condition for y coordinate.                    */
    int bc_z;    /**< Boundary condition for z coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real y_min;  /**< Minimum y coordinate in the grid.                       */
    real y_max;  /**< Maximum y coordinate in the grid.                       */
    real y_grid; /**< Interval between two adjacent points in y grid.         */
    real z_min;  /**< Minimum z coordinate in the grid.                       */
    real z_max;  /**< Maximum z coordinate in the grid.                       */
    real z_grid; /**< Interval between two adjacent points in z grid.         */
    real *c;     /**< Pointer to array with interpolant values.               */
} Linear3D;

#endif
