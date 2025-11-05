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
typedef enum
{
    NATURALBC,  /**< Second derivative is zero at the boundary.               */
    PERIODICBC, /**< Function has same value and derivatives at both ends.    */
} BoundaryCondition;

/**
 * Number of spline coefficients needed for each point in data.
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
    size_t nx;    /**< Number of grid points in x.                            */
    int xbc;      /**< Boundary condition in x.                               */
    real dx;      /**< Grid interval in x.                                    */
    real xlim[2]; /**< Grid limits in x.                                      */

    /**
     * Spline coefficients.
     *
     * Layout: [ix * 2 + c], where c = 0 for f and c = 1 for fxx.
     */
    real *c;
} Spline1D;

/**
 * Bicubic interpolation struct.
 */
typedef struct
{
    size_t nx;    /**< Number of grid points in x.                            */
    size_t ny;    /**< Number of grid points in y.                            */
    int xbc;      /**< Boundary condition in x.                               */
    int ybc;      /**< Boundary condition in y.                               */
    real dx;      /**< Grid interval in x.                                    */
    real dy;      /**< Grid interval in y.                                    */
    real xlim[2]; /**< Grid limits in x.                                      */
    real ylim[2]; /**< Grid limits in y.                                      */

    /**
     * Spline coefficients.
     *
     * Layout: [(ix * ny + iy) * 4 + c] (C order), where c = 0 for f, c = 1
     * for fxx, c = 2 for fyy, and c = 3 for fxy.
     */
    real *c;
} Spline2D;

/**
 * Tricubic interpolation struct.
 */
typedef struct
{
    size_t nx;    /**< Number of grid points in x.                            */
    size_t ny;    /**< Number of grid points in y.                            */
    size_t nz;    /**< Number of grid points in z.                            */
    int xbc;      /**< Boundary condition in x.                               */
    int ybc;      /**< Boundary condition in y.                               */
    int zbc;      /**< Boundary condition in z.                               */
    real dx;      /**< Grid interval in x.                                    */
    real dy;      /**< Grid interval in y.                                    */
    real dz;      /**< Grid interval in z.                                    */
    real xlim[2]; /**< Grid limits in x.                                      */
    real ylim[2]; /**< Grid limits in y.                                      */
    real zlim[2]; /**< Grid limits in z.                                      */

    /**
     * Spline coefficients.
     *
     * Layout: [(ix * ny * nz + iy * nz + iz) * 8 + c] (C order), where c = 0
     * for f, c = 1 for fxx, c = 2 for fyy, c = 3 for fzz, c = 4 for fxxyy,
     * c = 5 for fxxzz, c = 6 for fyyzz, and c = 7 for fxxyyzz.
     */
    real *c;
} Spline3D;

/**
 * 1D linear interpolation struct.
 */
typedef struct
{
    size_t nx;    /**< Number of grid points in x.                            */
    int xbc;      /**< Boundary condition in x.                               */
    real dx;      /**< Grid interval in x.                                    */
    real xlim[2]; /**< Grid limits in x.                                      */
    real *c;      /**< Interpolant values at the grid points (nx,).           */
} Linear1D;

/**
 * 2D linear interpolation struct.
 */
typedef struct
{
    size_t nx;    /**< Number of grid points in x.                            */
    size_t ny;    /**< Number of grid points in y.                            */
    int xbc;      /**< Boundary condition in x.                               */
    int ybc;      /**< Boundary condition in y.                               */
    real dx;      /**< Grid interval in x.                                    */
    real dy;      /**< Grid interval in y.                                    */
    real xlim[2]; /**< Grid limits in x.                                      */
    real ylim[2]; /**< Grid limits in y.                                      */
    real *c;      /**< Interpolant values at the grid points (nx,).           */
} Linear2D;

/**
 * 3D linear interpolation struct.
 */
typedef struct
{
    size_t nx;    /**< Number of grid points in x.                            */
    size_t ny;    /**< Number of grid points in y.                            */
    size_t nz;    /**< Number of grid points in z.                            */
    int xbc;      /**< Boundary condition in x.                               */
    int ybc;      /**< Boundary condition in y.                               */
    int zbc;      /**< Boundary condition in z.                               */
    real dx;      /**< Grid interval in x.                                    */
    real dy;      /**< Grid interval in y.                                    */
    real dz;      /**< Grid interval in z.                                    */
    real xlim[2]; /**< Grid limits in x.                                      */
    real ylim[2]; /**< Grid limits in y.                                      */
    real zlim[2]; /**< Grid limits in z.                                      */
    real *c;      /**< Interpolant values at the grid points (nx,).           */
} Linear3D;

#endif
