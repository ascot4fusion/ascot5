/**
 * @file boozer_data.h
 * Boozer data struct.
 */
#ifndef BOOZER_DATA_H
#define BOOZER_DATA_H

#include "interp_data.h"
#include <stddef.h>

/**
 * Data for mapping between the cylindrical and Boozer coordinates.
 *
 * Separatrix coordinates are needed in order not to evaluate Boozer coordinates
 * in private plasma region.
 */
typedef struct
{
    size_t nrz;     /**< Number of separatrix points.                         */
    real psilim[2]; /**< Minimum psi in other fields.                         */
    real *rlim;     /**< Separatrix R coordinates [m].                        */
    real *zlim;     /**< Separatrix z coordinates [m].                        */

    /**
     * Boozer poloidal angle as a function of psi and geometrical poloidal
     * angle [rad].
     */
    Spline2D theta;

    /**
     * Difference between the cylindrical and Boozer toroidal angle as a
     * function of psi and Boozer poloidal angle [rad].
     */
    Spline2D nu;
} Boozer;

#endif
