/**
 * @file linint.h
 * @brief Linear interpolation library
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
#include "../offload_acc_omp.h"
#include "../ascot5.h"
#include "../spline/interp.h"

/**
 * @brief 1D interpolation struct.
 */
typedef struct {
    int n_x;     /**< number of x grid points                        */
    int bc_x;    /**< boundary condition for x coordinate            */
    real x_min;  /**< minimum x coordinate in the grid               */
    real x_max;  /**< maximum x coordinate in the grid               */
    real x_grid; /**< interval between two adjacent points in x grid */
    real* c;     /**< pointer to array with interpolant values       */
} linint1D_data;

/**
 * @brief 2D interpolation struct.
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
    real* c;     /**< pointer to array with interpolant values       */
} linint2D_data;

/**
 * @brief 3D interpolation struct.
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
    real* c;     /**< pointer to array with interpolant values       */
} linint3D_data;

#pragma omp declare target
void linint1D_init(linint1D_data* str, real* c,
                   int n_x, int bc_x,
                   real x_min, real x_max);

void linint2D_init(linint2D_data* str, real* c,
                   int n_x, int n_y, int bc_x, int bc_y,
                   real x_min, real x_max,
                   real y_min, real y_max);

void linint3D_init(linint3D_data* str, real* c,
                   int n_x, int n_y, int n_z,
                   int bc_x, int bc_y, int bc_z,
                   real x_min, real x_max,
                   real y_min, real y_max,
                   real z_min, real z_max);

#pragma omp declare simd uniform(str)
DECLARE_TARGET
int linint1D_eval_f(real* f, linint1D_data* str, real x);
DECLARE_TARGET_END
#pragma omp declare simd uniform(str)
DECLARE_TARGET
int linint2D_eval_f(real* f, linint2D_data* str, real x, real y);
DECLARE_TARGET_END
#pragma omp declare simd uniform(str)
DECLARE_TARGET
int linint3D_eval_f(real* f, linint3D_data* str,
                    real x, real y, real z);
DECLARE_TARGET_END
#pragma omp end declare target
#endif
