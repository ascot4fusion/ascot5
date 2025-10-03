/**
 * @file interp.h
 * Spline interpolation library
 *
 * Spline interpolation fits cubic splines on data given on a uniform grid. Each
 * axis may have either natural or periodic boundary condition.
 *
 * There exists two representations for the splines: compact and explicit. Both
 * give identical results but the difference is that compact requires fewer
 * coefficients to be stored (and fetched from the memory at each evaluation)
 * than the explicit. On the other hand, explicit requires less computations
 * per evaluations, but compact is usually faster, and conserves memory,
 * especially in 3D. Therefore compact splines are preferred.
 *
 * To initialize splines, first call corresponding init_coeff() function which
 * evaluates coefficients to a pre-allocated array (i.e. to the offload array).
 * Then (after offloading is done) call init_spline() which assigns the
 * coefficients to a spline struct.
 *
 * In order to allocate the array for storing the coefficients, one needs to
 * know how many coefficients are stored per data grid point:
 *
 * - 1D compact  2, explicit 4
 * - 2D compact  4, explicit 16
 * - 3D compact  8, explicit 64
 */
#ifndef INTERP_H
#define INTERP_H
#include "ascot5.h"
#include "error.h"
#include "offload.h"

/**
 * Boundary conditions for the spline interpolation.
 */
enum boundaryCondition
{
    NATURALBC = 0, /**< Second derivative is zero at both ends.               */
    PERIODICBC = 1 /**< Function has same value and derivatives at both ends. */
};

/**
 * Number of coefficients stored for each data point.
 */
enum splinesize
{
    NSIZE_COMP1D = 2,
    NSIZE_COMP2D = 4,
    NSIZE_COMP3D = 8,
    NSIZE_EXPL1D = 4,
    NSIZE_EXPL2D = 16,
    NSIZE_EXPL3D = 64
};

/**
 * Cubic interpolation struct.
 */
typedef struct
{
    int n_x;     /**< Number of x grid points.                                */
    int bc_x;    /**< Boundary condition for x coordinate.                    */
    real x_min;  /**< Minimum x coordinate in the grid.                       */
    real x_max;  /**< Maximum x coordinate in the grid.                       */
    real x_grid; /**< Interval between two adjacent points in x grid.         */
    real *c;     /**< Pointer to array with spline coefficients.              */
} interp1D_data;

/**
 * Bicubic interpolation struct.
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
    real *c;     /**< Pointer to array with spline coefficients.              */
} interp2D_data;

/**
 * Tricubic interpolation struct.
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
    real *c;     /**< Pointer to array with spline coefficients.              */
} interp3D_data;

/**
 * Calculate cubic spline interpolation coefficients for scalar 1D data.
 *
 * This function calculates the cubic spline interpolation coefficients and
 * stores them in a pre-allocated array. Compact cofficients are calculated.
 *
 * @param c Allocated array of length n_x*2 to store the coefficients.
 * @param f Data to be interpolated.
 * @param n_x Number of data points in the x axis.
 * @param bc_x Boundary condition for the x axis.
 * @param x_min Minimum value of the x axis.
 * @param x_max Maximum value of the x axis.
 */
int interp1Dcomp_init_coeff(
    real *c, real *f, int n_x, int bc_x, real x_min, real x_max);

int interp2Dcomp_init_coeff(
    real *c, real *f, int n_x, int n_y, int bc_x, int bc_y, real x_min,
    real x_max, real y_min, real y_max);

int interp3Dcomp_init_coeff(
    real *c, real *f, int n_x, int n_y, int n_z, int bc_x, int bc_y, int bc_z,
    real x_min, real x_max, real y_min, real y_max, real z_min, real z_max);

int interp1Dexpl_init_coeff(
    real *c, real *f, int n_x, int bc_x, real x_min, real x_max);

int interp2Dexpl_init_coeff(
    real *c, real *f, int n_x, int n_y, int bc_x, int bc_y, real x_min,
    real x_max, real y_min, real y_max);

int interp3Dexpl_init_coeff(
    real *c, real *f, int n_x, int n_y, int n_z, int bc_x, int bc_y, int bc_z,
    real x_min, real x_max, real y_min, real y_max, real z_min, real z_max);

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
void interp1Dcomp_init_spline(
    interp1D_data *str, real *c, int n_x, int bc_x, real x_min, real x_max);

void interp2Dcomp_init_spline(
    interp2D_data *str, real *c, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

void interp3Dcomp_init_spline(
    interp3D_data *str, real *c, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

void interp1Dexpl_init_spline(
    interp1D_data *str, real *c, int n_x, int bc_x, real x_min, real x_max);

void interp2Dexpl_init_spline(
    interp2D_data *str, real *c, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

void interp3Dexpl_init_spline(
    interp3D_data *str, real *c, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

int interp1Dcomp_setup(
    interp1D_data *str, real *f, int n_x, int bc_x, real x_min, real x_max);

int interp2Dcomp_setup(
    interp2D_data *str, real *f, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

int interp3Dcomp_setup(
    interp3D_data *str, real *f, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dcomp_eval_f(real *f, interp1D_data *str, real x);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dcomp_eval_f(real *f, interp2D_data *str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dcomp_eval_f(real *f, interp3D_data *str, real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dexpl_eval_f(real *f, interp1D_data *str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dexpl_eval_f(real *f, interp2D_data *str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dexpl_eval_f(real *f, interp3D_data *str, real x, real y, real z);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dcomp_eval_df(real *f_df, interp1D_data *str, real x);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dcomp_eval_df(real *f_df, interp2D_data *str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dcomp_eval_df(
    real *f_df, interp3D_data *str, real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp1Dexpl_eval_df(real *f_df, interp1D_data *str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp2Dexpl_eval_df(real *f_df, interp2D_data *str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
a5err interp3Dexpl_eval_df(
    real *f_df, interp3D_data *str, real x, real y, real z);
#endif
