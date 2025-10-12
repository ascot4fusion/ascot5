/**
 * @file neutral_data.h
 * Neutral data types.
 */
#ifndef NEUTRAL_DATA_H
#define NEUTRAL_DATA_H

#include "interp_data.h"

/**
 * Neutral data types.
 */
typedef enum Neutral_type
{
    NEUTRAL_RADIAL = 1, /**< Corresponds to NeutralRadial.                    */
    NEUTRAL_ARBITRARY,  /**< Corresponds to NeutralArbitrary.                 */
} Neutral_type;

/**
 * Parameters for radial neutral distribution.
 */
typedef struct
{
    size_t n;     /**< Number of neutral species.                             */
    int *anum;    /**< Neutral species mass number.                           */
    int *znum;    /**< Neutral species charge number.                         */
    Linear1D *n0; /**< Density interpolation struct for each species.         */
    Linear1D *T0; /**< Temperature intepolation struct for each species.      */
} NeutralRadial;

/**
 * Parameters for arbitrary neutral distribution.
 */
typedef struct
{
    size_t n;     /**< Number of neutral species.                             */
    int *anum;    /**< Neutral species mass number.                           */
    int *znum;    /**< Neutral species charge number.                         */
    Linear3D *n0; /**< Density interpolation struct for each species.         */
    Linear3D *T0; /**< Temperature intepolation struct for each species.      */
} NeutralArbitrary;

/**
 * Neutral simulation data.
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the ``type`` field.
 */
typedef struct
{
    NeutralRadial *radial;       /**< Radial neutral distribution.            */
    NeutralArbitrary *arbitrary; /**< Arbitrary neutral distribution.         */
    Neutral_type type;           /**< Current neutral data type.              */
} Neutral;

#endif
