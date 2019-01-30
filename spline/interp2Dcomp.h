/**
 * @file interp2Dcomp.h
 * @brief Header file for interp2Dcomp.c
 */
#ifndef INTERP2DCOMP_H
#define INTERP2DCOMP_H
#include "../ascot5.h"

enum boundaryCondition {
    NATURALBC  = 0,
    PERIODICBC = 1
};

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

int interp2Dcomp_init_coeff(real* c, real* f,
                            int n_x, int n_y, int bc_x, int bc_y,
                            real x_min, real x_max,
                            real y_min, real y_max);

#pragma omp declare target
void interp2Dcomp_init_spline(interp2D_data* str, real* c,
                              int n_x, int n_y, int bc_x, int bc_y,
                              real x_min, real x_max,
                              real y_min, real y_max);
#pragma omp declare simd uniform(str)
int interp2Dcomp_eval_f(real* f, interp2D_data* str, real x, real y);
#pragma omp declare simd uniform(str)
int interp2Dcomp_eval_df(real* f_df, interp2D_data* str, real x, real y);
#pragma omp end declare target
#endif
