/**
 * @file mhd_data.h
 * MHD input data types.
 */
#ifndef MHD_DATA_H
#define MHD_DATA_H

#include "interp_data.h"

/**
 * MHD input data types.
 */
typedef enum Mhd_type
{
    MHD_STATIONARY = 1, /**< Corresponds to MhdStationary.                    */
    MHD_DYNAMIC         /**< Corresponds to MhdDynamic.                       */
} Mhd_type;

/**
 * MHD data with stationary eigenmodes.
 *
 * The eigenmodes are a function of rho.
 */
typedef struct
{
    size_t n;        /**< Number of modes.                                    */
    int *nmode;      /**< Toroidal mode numbers.                              */
    int *mmode;      /**< Poloidal mode numbers.                              */
    real *amplitude; /**< Mode amplitudes.                                    */
    real *omega;     /**< Mode frequency [rad/s].                             */
    real *phase;     /**< Mode phase [rad].                                   */
    Spline1D *alpha; /**< Eigenmodes for magnetic perturbation [m].           */
    Spline1D *phi;   /**< Eigenmodes for electric perturbation [V].           */
} MhdStationary;

/**
 * MHD parameters
 */
typedef struct
{
    size_t n;        /**< Number of modes.                                    */
    int *nmode;      /**< Toroidal mode numbers.                              */
    int *mmode;      /**< Poloidal mode numbers.                              */
    real *amplitude; /**< Mode amplitudes.                                    */
    real *omega;     /**< Mode frequency [rad/s].                             */
    real *phase;     /**< Mode phase [rad].                                   */
    Spline2D *alpha; /**< Eigenmodes for magnetic perturbation [m].           */
    Spline2D *phi;   /**< Eigenmodes for electric perturbation [V].           */
} MhdDynamic;

/**
 * MHD simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the ``type`` field.
 */
typedef struct
{
    MhdStationary *stationary; /**< Stationary eigenmodes.                    */
    MhdDynamic *dynamic;       /**< Dynamic eigenmodes.                       */
    Mhd_type type;             /**< Current Mhd input type.                   */
} Mhd;

#endif
