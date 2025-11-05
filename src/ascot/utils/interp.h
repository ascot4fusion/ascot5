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
#include "interp_data.h"
#include "parallel.h"
#include <stddef.h>

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
void Linear1D_init(
    Linear1D *linear, size_t nx, int xbc, real xlim[2], real c[nx]);

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
void Linear2D_init(
    Linear2D *linear, size_t nx, size_t ny, int xbc, int ybc, real xlim[2],
    real ylim[2], real c[nx * ny]);

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
void Linear3D_init(
    Linear3D *linear, size_t nx, size_t ny, size_t nz, int xbc, int ybc,
    int zbc, real xlim[2], real ylim[2], real zlim[2], real c[nx * ny * nz]);

GPU_DECLARE_TARGET_SIMD_UNIFORM(linear)
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
int Linear1D_eval_f(real f[1], Linear1D *linear, real x);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(linear)
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
int Linear2D_eval_f(real f[1], Linear2D *linear, real x, real y);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(linear)
/**
 * Evaluate interpolated value of 3D scalar field.
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
int Linear3D_eval_f(real f[1], Linear3D *linear, real x, real y, real z);
DECLARE_TARGET_END

/**
 * Initialize 1D spline interpolant from scalar data.
 *
 * @param spline Spline to initialize.
 * @param nx Number of grid points in x.
 * @param xbc Boundary condition for x.
 * @param xlim Grid limits in x.
 * @param f 1D data to be interpolated.
 *
 * @return Zero if the initialization succeeded.
 */
int Spline1D_init(
    Spline1D *spline, size_t nx, BoundaryCondition xbc, real xlim[2],
    real f[nx]);

/**
 * Initialize 2D spline interpolant from scalar data.
 *
 * @param spline Spline to initialize.
 * @param nx Number of grid points in x.
 * @param ny Number of grid points in y.
 * @param xbc Boundary condition for x.
 * @param ybc Boundary condition for y.
 * @param xlim Grid limits in x.
 * @param ylim Grid limits in y.
 * @param f 2D data to be interpolated.
 *        Layout: [i_x*ny + i_y] (C order).
 *
 * @return Zero if the initialization succeeded.
 */
int Spline2D_init(
    Spline2D *spline, size_t nx, size_t ny, BoundaryCondition xbc,
    BoundaryCondition ybc, real xlim[2], real ylim[2], real f[nx * ny]);

/**
 * Initialize 3D spline interpolant from scalar data.
 *
 * @param spline Spline to initialize.
 * @param nx Number of grid points in x.
 * @param ny Number of grid points in y.
 * @param nz Number of grid points in z.
 * @param xbc Boundary condition for x.
 * @param ybc Boundary condition for y.
 * @param zbc Boundary condition for z.
 * @param xlim Grid limits in x.
 * @param ylim Grid limits in y.
 * @param zlim Grid limits in z.
 * @param f 3D data to be interpolated.
 *        Layout: [i_x*ny*nz + i_y*nz + i_z] (C order)
 *
 * @return Zero if the initialization succeeded.
 */
int Spline3D_init(
    Spline3D *spline, size_t nx, size_t ny, size_t nz, BoundaryCondition xbc,
    BoundaryCondition ybc, BoundaryCondition zbc, real xlim[2], real ylim[2],
    real zlim[2], real f[nx * ny * nz]);

GPU_DECLARE_TARGET_SIMD_UNIFORM(spline)
/**
 * Interpolate 1D data using cubic splines.
 *
 * @param f Interpolated value.
 * @param spline The spline interpolant.
 * @param x Query point x coordinate.
 *
 * @return Zero on success and one if the point is outside the domain.
 */
int Spline1D_eval_f(real f[1], Spline1D *spline, real x);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(spline)
/**
 * Interpolate 2D data using cubic splines.
 *
 * @param f Interpolated value.
 * @param spline The spline interpolant.
 * @param x Query point x coordinate.
 * @param y Query point y coordinate.
 *
 * @return Zero on success and one if the point is outside the domain.
 */
int Spline2D_eval_f(real f[1], Spline2D *spline, real x, real y);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(spline)
/**
 * Interpolate 3D data using cubic splines.
 *
 * @param f Interpolated value.
 * @param spline The spline interpolant.
 * @param x Query point x coordinate.
 * @param y Query point y coordinate.
 * @param y Query point z coordinate.
 *
 * @return Zero on success and one if the point is outside the domain.
 */
int Spline3D_eval_f(real f[1], Spline3D *spline, real x, real y, real z);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(spline)
/**
 * Interpolate 1D data and the first and second order derivatives using cubic
 * splines.
 *
 * @param f_df Interpolated value and derivatives.
 *        Layout: [f, f_x, f_xx].
 * @param spline The spline interpolant.
 * @param x Query point x coordinate.
 *
 * @return Zero on success and one if the point is outside the domain.
 */
int Spline1D_eval_f_df(real f_df[3], Spline1D *spline, real x);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(spline)
/**
 * Interpolate 2D data and the first and second order derivatives using cubic
 * splines.
 *
 * @param f_df Interpolated value and derivatives.
 *        Layout: [f, f_x, f_y, f_xx, f_yy, f_xy].
 * @param spline The spline interpolant.
 * @param x Query point x coordinate.
 * @param y Query point y coordinate.
 *
 * @return Zero on success and one if the point is outside the domain.
 */
int Spline2D_eval_f_df(real f_df[6], Spline2D *spline, real x, real y);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(spline)
/**
 * Interpolate 3D data and the first and second order derivatives using cubic
 * splines.
 *
 * @param f_df Interpolated value and derivatives.
 *        Layout: [f, f_x, f_y, f_z, f_xx, f_yy, f_zz, f_xy, f_xz, f_yz].
 * @param spline The spline interpolant.
 * @param x Query point x coordinate.
 * @param y Query point y coordinate.
 * @param z Query point z coordinate.
 *
 * @return Zero on success and one if the point is outside the domain.
 */
int Spline3D_eval_f_df(real f_df[10], Spline3D *spline, real x, real y, real z);
DECLARE_TARGET_END

/**
 * Calculate compact cubic spline interpolation coefficients in 1D for uniformly
 * spaced data.
 *
 * This function calculates the compact cubic interpolation coefficients for a
 * 1D dataset using one of the boundary conditions.
 *
 * For natural boundary condition, the algorithm enforces
 *
 * f''_0 = f''_{n-1} = 0.
 *
 * The internal spline coefficients satisfy the tridiagonal system:
 *
 *   f''_{i-1} + 4f''_i + f''_{i+1}
 *   = 6\,(f_{i+1} - 2f_i + f_{i-1}), \qquad i = 1, \dots, n-2.
 *
 * This system is solved using a simplified **Thomas algorithm**.
 *
 * For periodic boundary condition, the spline enforces continuity of the
 * function and its first and second derivatives at the endpoints:
 *
 *   f_0 = f_{n},
 *   f'_0 = f'_{n},
 *   f''_0 = f''_{n}
 *
 *
 * The resulting **cyclic tridiagonal system** is:
 *
 *  4f''_0 + f''_1 + f''_{n-2} &= 6(f_1 - 2f_0 + f_{n-2}) \\
 *  f''_{i-1} + 4f''_i + f''_{i+1} &= 6(f_{i+1} - 2f_i + f_{i-1}), \quad i = 1,
 * \dots, n-3 \\ f''_{n-3} + 4f''_{n-2} + f''_0 &= 6(f_0 - 2f_{n-2} + f_{n-3})
 *
 * This system is solved by a custom elimination algorithm derived from Gaussian
 * elimination adapted for periodic constraints.
 *
 * @param n Number of data points.
 * @param c Preallocated array for storing the calculated coefficients.
 *        Layout: c[2*i] = f_i, c[2*i+1] = f''_i.
 * @param f Data to be interpolated.
 * @param bc Boundary condition.
 */
int solve_compact_cubic_spline(size_t n, real c[2 * n], real f[n], int bc);

#endif
