/**
 * @file diag_data.h
 * Diagnostics data structures.
 */
#ifndef DIAG_DATA_H
#define DIAG_DATA_H

#include "defines.h"
#include <stddef.h>

/**
 * Orbit diagnostics data.
 *
 */
typedef struct
{

    /** Marker R coordinate (layout: [imrk*npoint + ipnt]) [m].               */
    real *r;

    /** Marker z coordinate (layout: [imrk*npoint + ipnt]) [m].               */
    real *z;

    /** Marker phi coordinate (layout: [imrk*npoint + ipnt]) [rad].           */
    real *phi;

    /**
     * First momentum coordinate (layout: [imrk*npoint + ipnt]).
     *
     * - pr for gyro-orbit.
     * - ppara for guiding-center.
     */
    real *p1;

    /**
     * Second momentum coordinate (layout: [imrk*npoint + ipnt]).
     *
     * - pphi for gyro-orbit.
     * - mu for guiding-center.
     */
    real *p2;

    /**
     * Third momentum coordinate (layout: [imrk*npoint + ipnt]).
     *
     * - pz for gyro-orbit.
     * - gyro-angle for guiding-center.
     */
    real *p3;

    /** Marker mileage (layout: [imrk*npoint + ipnt]) [s].                    */
    real *mileage;

    /** Time marker was last updated (layout: [imrk]) [s].                    */
    real *stamp;

    /** Marker ID (layout: [imrk*npoint + ipnt]) [s].                         */
    size_t *id;

    /** Index (ipoint) of the last recorded point (layout: [imrk]) [s].       */
    size_t *idx;

    /** Marker charge state (layout: [imrk*npoint + ipnt]) [e].               */
    int *charge;

    /**
     * Poincaré flag (layout: [imrk*npoint + ipnt]).
     *
     * This flag identifies which Poincaré plane this data point corresponds to.
     * The index starts with zero.
     *
     * - Indices 0 <= i < ntoroidal correspond to toroidal Poincaré planes.
     * - Indices ntoroidal <= i < ntoroidal + npoloidal correspond to poloidal
     *   Poincaré planes.
     * - Indices ntoroidal + npoloidal <= i correspond to radial Poincaré
     *   surfaces.
     *
     * In addition this flag is assigned sign that indicates whether the
     * particle was travelling along the positive or negative direction
     * poloidally.
     */
    int *poincare;

    /**
     * Simulation mode when the record was made (layout: [imrk*npoint + ipnt]).
     *
     * This flag is required to store orbits in hybrid mode consistently.
     */
    int *simmode;

    /** Number of points to record.                                           */
    size_t npoint;

    /** Number of toroidal Poincaré planes.                                   */
    size_t ntoroidal;

    /** Number of poloidal Poincaré planes.                                   */
    size_t npoloidal;

    /** Number of radial Poincaré surfaces.                                   */
    size_t nradial;

    /**
     * Interval for recording orbit [s].
     *
     * Negative to record Poincaré plane crossings instead.
     */
    real interval;

    /** Poloidal angles of the toroidal Poincaré planes (ntoroidal) [rad].    */
    real *toroidal;

    /** Toroidal angles of the poloidal Poincaré planes (npoloidal) [rad].    */
    real *poloidal;

    /** Poloidal angles of the toroidal Poincaré planes (ntoroidal) [rad].    */
    real *radial;
} DiagOrbit;

/**
 * Available histogram axis coordinates.
 */
typedef enum
{
    HIST_R,      /**< The R coordinate in cylindrical basis [m].              */
    HIST_PHI,    /**< The phi coordinate in cylindrical basis [rad].          */
    HIST_Z,      /**< The z coordinate in cylindrical basis [m].              */
    HIST_RHO,    /**< Square root of normalized poloidal flux [1].            */
    HIST_THETA,  /**< Poloidal angle [rad].                                   */
    HIST_PPAR,   /**< Momentum parallel to the magnetic field [kg*m/s].       */
    HIST_PPERP,  /**< Momentum orthogonal to the magnetic field [kg*m/s].     */
    HIST_PR,     /**< Momentum R-component [kg*m/s].                          */
    HIST_PPHI,   /**< Momentum phi-component [kg*m/s].                        */
    HIST_PZ,     /**< Momentum z-component [kg*m/s].                          */
    HIST_EKIN,   /**< Kinetic energy [J].                                     */
    HIST_XI,     /**< Pitch [1].                                              */
    HIST_MU,     /**< Magnetic moment [J/T].                                  */
    HIST_PTOR,   /**< Canonical toroidal angular momentum [kg*m/s].           */
    HIST_TIME,   /**< Time instant (laboratory time) [s].                     */
    HIST_CHARGE, /**< Charge state [e].                                       */
    HIST_NDIM,   /**< Number of coordinates.                                  */
} HistCoordinate;

/**
 * Coordinate axis for the histogram.
 */
typedef struct
{
    real min;  /**< Lower limit of the coordinate interval.                   */
    real max;  /**< Upper limit of the coordinate interval.                   */
    size_t n;  /**< Number of bins in this axis.                              */
    int coord; /**< Coordinate mapped to this axis.                           */
} HistAxis;

/**
 * Histogram data.
 *
 * The bins are stored as a flattened array. The element (i0,i1,...,in)
 * corresponds to bins[i0*strides[0] + i1*strides[1] + ... + in].
 */
typedef struct
{
    /**
     * Total number of bins in the histogram.
     */
    size_t nbin;

    /**
     * Row length for each dimension.
     *
     * The strides are used to access the bins in the flattened array. For
     * example: ``bins[i0*strides[0] + i1*strides[1] + ... + in]``, where the
     * fastest running index corresponds to ``HIST_R``, and otherwise the
     * strides are in the same order as in the enum ``HistCoordinate``.
     */
    size_t strides[HIST_NDIM - 1];

    /**
     * The coordinate axes.
     *
     * All possible axes are included, but those that are not used can be set as
     * dummies by setting their number of bins to zero. The axes are listed here
     * in the same order as in the enum ``HistCoordinate``.
     */
    HistAxis axes[HIST_NDIM];

    /**
     * The histogram bins.
     *
     * Each bin stores the weighted marker count. The layout is:
     * f(x0, x1, ..., xn) = bins[i0*strides[0] + i1*strides[1] + ... + in]
     * (C order).
     */
    real *bins;
} DiagHist;

/**
 * Diagnostics data struct.
 */
typedef struct
{
    DiagOrbit *orbit; /**< Orbit diagnostics.                                 */
    DiagHist *hist;   /**< List of histogram diagnostics to collect.          */
    size_t nhist;     /**< Number of histogram diagnostics.                   */
} Diagnostics;

#endif