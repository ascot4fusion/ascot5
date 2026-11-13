/**
 * Data structures for atomic data.
 */
#ifndef ATOMIC_DATA_H
#define ATOMIC_DATA_H

#include "interp_data.h"

/**
 * Interpolator for atomic reaction rate coefficients.
 *
 * Union is used because electron impact ionization rate coefficient is given
 * in terms of the temperature, ion impact ionization and charge-exchange rate
 * coefficients are in terms of temperature and beam energy, and beam stopping
 * rate coefficient is in terms of temperature, beam energy and density. Also
 * any of the abscissae can be logarithmic. The value of the function is
 * logarithmic always.
 */
typedef struct
{
    int logscale[3];       /**< Logarithmic scale flags for x, y and z.       */
    union {
        Linear1D linear1d; /**< 1D (T) linear interpolator.                   */
        Linear2D linear2d; /**< 2D (T, E) linear interpolator.                */
        Linear3D linear3d; /**< 3D (T, E, n) linear interpolator.             */
    };
} AtomicInterpolator;

/**
 * Atomic reaction rate coefficients.
 */
typedef struct
{
    AtomicInterpolator *EI;  /**< Electron-impact ionization interpolator.    */
    AtomicInterpolator *II;  /**< Ion-impact ionization interpolator.         */
    AtomicInterpolator *CX;  /**< Charge-exchange interpolator.               */
    AtomicInterpolator *BMS; /**< Effective beam stopping interpolator.       */
} Atomic;

#endif
