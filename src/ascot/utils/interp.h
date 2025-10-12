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
#include "defines.h"
#include "parallel.h"
#include "interp_data.h"
#include <stddef.h>

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
    real *c, real *f, int n_x, int bc_x);

int interp2Dexpl_init_coeff(
    real *c, real *f, int n_x, int n_y, int bc_x, int bc_y);

int interp3Dexpl_init_coeff(
    real *c, real *f, int n_x, int n_y, int n_z, int bc_x, int bc_y, int bc_z);

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
    Spline1D *str, real *c, int n_x, int bc_x, real x_min, real x_max);

void interp2Dcomp_init_spline(
    Spline2D *str, real *c, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

void interp3Dcomp_init_spline(
    Spline3D *str, real *c, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

void interp1Dexpl_init_spline(
    Spline1D *str, real *c, int n_x, int bc_x, real x_min, real x_max);

void interp2Dexpl_init_spline(
    Spline2D *str, real *c, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

void interp3Dexpl_init_spline(
    Spline3D *str, real *c, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

int interp1Dcomp_setup(
    Spline1D *str, real *f, int n_x, int bc_x, real x_min, real x_max);

int interp2Dcomp_setup(
    Spline2D *str, real *f, int n_x, int n_y, int bc_x, int bc_y,
    real x_min, real x_max, real y_min, real y_max);

int interp3Dcomp_setup(
    Spline3D *str, real *f, int n_x, int n_y, int n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp1Dcomp_eval_f(real *f, Spline1D *str, real x);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp2Dcomp_eval_f(real *f, Spline2D *str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp3Dcomp_eval_f(real *f, Spline3D *str, real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp1Dexpl_eval_f(real *f, Spline1D *str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp2Dexpl_eval_f(real *f, Spline2D *str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp3Dexpl_eval_f(real *f, Spline3D *str, real x, real y, real z);

GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp1Dcomp_eval_df(real *f_df, Spline1D *str, real x);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp2Dcomp_eval_df(real *f_df, Spline2D *str, real x, real y);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp3Dcomp_eval_df(
    real *f_df, Spline3D *str, real x, real y, real z);
DECLARE_TARGET_END

DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp1Dexpl_eval_df(real *f_df, Spline1D *str, real x);
DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp2Dexpl_eval_df(real *f_df, Spline2D *str, real x, real y);
DECLARE_TARGET_SIMD_UNIFORM(str)
err_t interp3Dexpl_eval_df(
    real *f_df, Spline3D *str, real x, real y, real z);

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
int linint1D_eval_f(real *f, Linear1D *str, real x);
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
int linint2D_eval_f(real *f, Linear2D *str, real x, real y);
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
int linint3D_eval_f(real *f, Linear3D *str, real x, real y, real z);
DECLARE_TARGET_END

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
    Linear1D *str, real *c, size_t n_x, int bc_x, real x_min, real x_max);

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
    Linear2D *str, real *c, size_t n_x, size_t n_y, int bc_x, int bc_y,
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
    Linear3D *str, real *c, size_t n_x, size_t n_y, size_t n_z, int bc_x, int bc_y,
    int bc_z, real x_min, real x_max, real y_min, real y_max, real z_min,
    real z_max);


/**
 * Calculate compact cubic spline interpolation coefficients in 1D.
 *
 * This function calculates the compact cubic interpolation coefficients for a
 * 1D data set using one of the boundary conditions.
 *
 * @param f Data to be interpolated.
 * @param n Number of data points.
 * @param bc Boundary condition.
 * @param c Calculated coefficients.
 */
void splinecomp(real* f, size_t n, int bc, real c[2*n]);

#endif
