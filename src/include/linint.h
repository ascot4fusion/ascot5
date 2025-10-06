/**
 * Linear interpolation library.
 *
 * Linear interpolation interpolates data given on a uniform grid.
 * Each axis may have either natural or periodic boundary condition.
 *
 * The interpolant does not need any initialization on the host. After
 * offloading is done, call linintXD_init() which assigns the values and
 * a pointer to the data to a linint struct.
 */
#ifndef LININT_H
#define LININT_H
#include "defines.h"
#include "parallel.h"
#include "interp.h"

/**
 * 1D linear interpolation struct.
 */
typedef struct
{
    int n_x;     /**< Number of x grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real *c;     /**< Pointer to array with interpolant values.               */
} linint1D_data;

/**
 * @brief 2D linear interpolation struct.
 */
typedef struct
{
    int n_x;     /**< Number of x grid points.                                */
    int n_y;     /**< Number of y grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    int bc_y;    /**< Boundary condition for y coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real y_min;  /**< Minimum y coordinate in the grid.                       */
    real y_max;  /**< Maximum y coordinate in the grid.                       */
    real y_grid; /**< Interval between two adjacent points in y grid.         */
    real *c;     /**< Pointer to array with interpolant values.               */
} linint2D_data;

/**
 * @brief 3D linear interpolation struct.
 */
typedef struct
{
    int n_x;     /**< Number of x grid points.                                */
    int n_y;     /**< Number of y grid points.                                */
    int n_z;     /**< Number of z grid points.                                */
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
} linint3D_data;

/**
 * Initialize linear interpolation struct for scalar 1D data.
 *
 * @param str Pointer to struct to be initialized.
 * @param c Array where data is stored.
 * @param n_x Number of data points in the x direction.
 * @param bc_x Boundary condition for x axis.
 * @param x_min Minimum value of the x axis.
 * @param x_max Maximum value of the x axis.
 */
void linint1D_init(
    linint1D_data *str, real *c, int n_x, int bc_x, real x_min, real x_max);

/**
 * Initialize linear interpolation struct for scalar 2D data.
 *
 * @param str Pointer to struct to be initialized.
 * @param c Array where data is stored.
 * @param n_x Number of data points in the x direction.
 * @param n_y Number of data points in the y direction.
 * @param bc_x Boundary condition for x axis.
 * @param bc_y Boundary condition for y axis.
 * @param x_min Minimum value of the x axis.
 * @param x_max Maximum value of the x axis.
 * @param y_min Minimum value of the y axis.
 * @param y_max Maximum value of the y axis.
 */
void linint2D_init(
    linint2D_data *str, real *c, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

/**
 * Initialize linear interpolation struct for scalar 3D data.
 *
 * @param str Pointer to struct to be initialized.
 * @param c Array where data is stored.
 * @param n_x Number of data points in the x direction.
 * @param n_y Number of data points in the y direction.
 * @param n_z Number of data points in the z direction.
 * @param bc_x Boundary condition for x axis.
 * @param bc_y Boundary condition for y axis.
 * @param bc_z Boundary condition for z axis.
 * @param x_min Minimum value of the x axis.
 * @param x_max Maximum value of the x axis.
 * @param y_min Minimum value of the y axis.
 * @param y_max Maximum value of the y axis.
 * @param z_min Minimum value of the z axis.
 * @param z_max Maximum value of the z axis.
 */
void linint3D_init(
    linint3D_data *str, real *c, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

/**
 * Evaluate interpolated value of 1D scalar field.
 *
 * This function evaluates the interpolated value of a 1D scalar field using
 * linear interpolation.
 *
 * @param f Variable in which to place the evaluated value.
 * @param str Data struct for data interpolation.
 * @param x The query point x coordinate.
 *
 * @return Zero on success and one if the point is outside the grid.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
int linint1D_eval_f(real *f, linint1D_data *str, real x);
DECLARE_TARGET_END

/**
 * Evaluate interpolated value of 2D scalar field.
 *
 * This function evaluates the interpolated value of a 2D scalar field using
 * bilinear interpolation.
 *
 * @param f variable in which to place the evaluated value.
 * @param str data struct for data interpolation.
 * @param x The query point x coordinate.
 * @param y The query point y coordinate.
 *
 * @return Zero on success and one if the point is outside the grid.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
int linint2D_eval_f(real *f, linint2D_data *str, real x, real y);
DECLARE_TARGET_END

/**
 * @brief Evaluate interpolated value of 3D scalar field.
 *
 * This function evaluates the interpolated value of a 3D scalar field using
 * trilinear interpolation.
 *
 * @param f Variable in which to place the evaluated value.
 * @param str Data struct for data interpolation.
 * @param x The query point x coordinate.
 * @param y The query point y coordinate.
 * @param z The query point z coordinate.
 */
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
int linint3D_eval_f(real *f, linint3D_data *str, real x, real y, real z);
DECLARE_TARGET_END

#endif
