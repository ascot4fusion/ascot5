/**
 * @file hist.h
 * @brief Header file for hist.c
 */
#ifndef HIST_H
#define HIST_H

#include <stdlib.h>
#include "../ascot5.h"
#include "../particle.h"

#define HIST_ALLDIM 16

/**
 * @brief Quantities that can be used as histogram axis coordinates.
 */
typedef enum {
    R,      /**< The R coordinate in cylindrical basis [m].                   */
    PHI,    /**< The phi coordinate in cylindrical basis [rad].               */
    Z,      /**< The z coordinate in cylindrical basis [m].                   */
    RHO,    /**< Square root of normalized poloidal flux [1].                 */
    THETA,  /**< Poloidal angle [rad].                                        */
    PPAR,   /**< Momentum component parallel to the magnetic field [kg*m/s].  */
    PPERP,  /**< Momentum component orthogonal to the magnetic field [kg*m/s].*/
    PR,     /**< Momentum R-component [kg*m/s].                               */
    PPHI,   /**< Momentum phi-component [kg*m/s].                             */
    PZ,     /**< Momentum z-component [kg*m/s].                               */
    EKIN,   /**< Kinetic energy [J].                                          */
    XI,     /**< Pitch [1].                                                   */
    MU,     /**< Magnetic moment [J/T].                                       */
    PTOR,   /**< Canonical toroidal angular momentum [kg*m/s].                */
    TIME,   /**< Time instant (laboratory time) [s].                          */
    CHARGE, /**< Charge state [e].                                            */
} hist_coordinate;

/**
 * @brief Coordinate axis for the histogram.
 */
typedef struct {
    hist_coordinate name; /**< Coordinate mapped to this axis                 */
    real min;             /**< Lower limit of the coordinate interval         */
    real max;             /**< Upper limit of the coordinate interval         */
    size_t n;             /**< Number of bins in this axis                    */
} hist_axis;

/**
 * @brief Histogram parameters.
 *
 * The bins are stored as a flattened array. The element (i0,i1,...,in)
 * corresponds to bins[i0*strides[0] + i1*strides[1] + ... + in].
 */
typedef struct {
    hist_axis axes[HIST_ALLDIM];   /**< The coordinate axes.                  */
    size_t strides[HIST_ALLDIM-1]; /**< Row length for each dimension.        */
    size_t nbin;                   /**< Number of bins.                       */
    real* bins;                    /**< The bin array.                        */
} histogram;

int hist_init(histogram* data, int dimensions, hist_coordinate* coordinates,
              real* binmin, real* binmax, size_t* nbin);
void hist_free(histogram* data);
void hist_offload(histogram* data);
void hist_update_fo(histogram* hist, particle_simd_fo* p_f,
                    particle_simd_fo* p_i);
void hist_update_gc(histogram* hist, particle_simd_gc* p_f,
                    particle_simd_gc* p_i);
#endif
