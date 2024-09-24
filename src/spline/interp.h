/**
 * @file interp.h
 * @brief Spline interpolation library
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
#include "../offload_acc_omp.h"
#include "../ascot5.h"
#include "../error.h"

/**
 * @brief Boundary conditions for the spline interpolation.
 */
enum boundaryCondition {
    NATURALBC  = 0, /**< Second derivative is zero at both ends               */
    PERIODICBC = 1  /**< Function has same value and derivatives at both ends */
};

/**
 * @brief Number of coefficients stored for each data point.
 */
enum splinesize {
    NSIZE_COMP1D =  2,
    NSIZE_COMP2D =  4,
    NSIZE_COMP3D =  8,
    NSIZE_EXPL1D =  4,
    NSIZE_EXPL2D = 16,
    NSIZE_EXPL3D = 64
};

/**
 * @brief Cubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp1D_data;

/**
 * @brief Bicubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int n_y;     /**< number of y grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real y_min;  /**< minimum y coordinate in the grid               */
    real y_max;  /**< maximum y coordinate in the grid               */
    real y_grid; /**< interval between two adjacent points in y grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp2D_data;

/**
 * @brief Tricubic interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int n_y;     /**< number of y grid points                        */
    int n_z;     /**< number of z grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    int bc_y;    /**< boundary condition for y coordinate            */
    int bc_z;    /**< boundary condition for z coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real y_min;  /**< minimum y coordinate in the grid               */
    real y_max;  /**< maximum y coordinate in the grid               */
    real y_grid; /**< interval between two adjacent points in y grid */
    real z_min;  /**< minimum z coordinate in the grid               */
    real z_max;  /**< maximum z coordinate in the grid               */
    real z_grid; /**< interval between two adjacent points in z grid */
    real* c;     /**< pointer to array with spline coefficients      */
} interp3D_data;

int interp1Dcomp_init_coeff(real* c, real* f,
                            int n_x, int bc_x,
                            real x_min, real x_max);

int interp2Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

int interp3Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max);

int interp1Dexpl_init_coeff(real* c, real* f,
                            int n_x, int bc_x,
                            real x_min, real x_max);

int interp2Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

int interp3Dexpl_init_coeff(real* c, real* f,
                            int n_x, int n_y, int n_z,
                            int bc_x, int bc_y, int bc_z,
                            real x_min, real x_max,
                            real y_min, real y_max,
                            real z_min, real z_max);

void interp1Dcomp_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x,
                              real x_min, real x_max);

void interp2Dcomp_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);

void interp3Dcomp_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max);

void interp1Dexpl_init_spline(interp1D_data* str, real* c,
                              int n_x, int bc_x,
                              real x_min, real x_max);

void interp2Dexpl_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);

void interp3Dexpl_init_spline(interp3D_data* str, real* c,
                              int n_x, int n_y, int n_z,
                              int bc_x, int bc_y, int bc_z,
                              real x_min, real x_max,
                              real y_min, real y_max,
                              real z_min, real z_max);

#ifndef GPU
#pragma omp declare simd uniform(str)
#else
DECLARE_TARGET
#endif
a5err interp1Dcomp_eval_f(real* f, interp1D_data* str, real x);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(str)
#else
DECLARE_TARGET
#endif
a5err interp2Dcomp_eval_f(real* f, interp2D_data* str, real x, real y);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(str)
#else
DECLARE_TARGET
#endif
a5err interp3Dcomp_eval_f(real* f, interp3D_data* str,
                         real x, real y, real z);
DECLARE_TARGET_END

#pragma omp declare simd uniform(str)
a5err interp1Dexpl_eval_f(real* f, interp1D_data* str, real x);
#pragma omp declare simd uniform(str)
a5err interp2Dexpl_eval_f(real* f, interp2D_data* str, real x, real y);
#pragma omp declare simd uniform(str)
a5err interp3Dexpl_eval_f(real* f, interp3D_data* str,
                          real x, real y, real z);

#ifndef GPU
#pragma omp declare simd uniform(str)
#else
DECLARE_TARGET
#endif
a5err interp1Dcomp_eval_df(real* f_df, interp1D_data* str, real x);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(str)
#else
DECLARE_TARGET
#endif
a5err interp2Dcomp_eval_df(real* f_df, interp2D_data* str, real x, real y);
DECLARE_TARGET_END
#ifndef GPU
#pragma omp declare simd uniform(str)
#else
DECLARE_TARGET
#endif
a5err interp3Dcomp_eval_df(real* f_df, interp3D_data* str,
                           real x, real y, real z);
DECLARE_TARGET_END

#pragma omp declare simd uniform(str)
a5err interp1Dexpl_eval_df(real* f_df, interp1D_data* str, real x);
#pragma omp declare simd uniform(str)
a5err interp2Dexpl_eval_df(real* f_df, interp2D_data* str, real x, real y);
#pragma omp declare simd uniform(str)
a5err interp3Dexpl_eval_df(real* f_df, interp3D_data* str,
                           real x, real y, real z);
#endif
