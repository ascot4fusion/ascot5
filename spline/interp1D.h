/**
 * @file interp1D.h
 * @brief Header file for 1D interpolation
 */
#ifndef INTERP1D_H
#define INTERP1D_H

/**
 * @brief Interpolation struct
 */
typedef struct {
    int n_r;                  /**< number of r grid points */
    real r_min;               /**< minimum r coordinate in the grid */
    real r_max;               /**< r grid interval (r_max-r_min)/(n_r-1) */
    real r_grid;              /**< r grid interval (r_max-r_min)/(n_r-1) */
    real* c;                  /**< pointer to array with spline coefficients */
} interp1D_data;

#endif
