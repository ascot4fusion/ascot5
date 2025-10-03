/**
 * @file neutral.h
 * Neutral interface
 *
 * This is an interface through which neutral data is initialized and accessed.
 * Reading e.g. from disk is done elsewhere.
 *
 * To add a new neutral data instance, make sure these functions are implemented
 * and called from this interface, and that neutral.h contains enum type for the
 * new instance.
 *
 * The interface checks which instance given data corresponds to from the
 * "type"-field in neutral_offload_data or neutral_data that is given as an
 * argument, and calls the relevant function for that instance.
 */
#ifndef NEUTRAL_H
#define NEUTRAL_H

#include "ascot5.h"
#include "error.h"
#include "linint.h"
#include "offload.h"

/**
 * @brief Neutral data types
 */
typedef enum neutral_type
{
    neutral_type_1D, /**< Linearly-interpolated 1D neutral data          */
    neutral_type_3D, /**< Linearly-interpolated 3D neutral data          */
} neutral_type;

/**
 * @brief 1D neutral parameters on the target
 */
typedef struct
{
    int n_species;     /**< Number of neutral species                         */
    int *anum;         /**< Neutral species mass number                       */
    int *znum;         /**< Neutral species charge number                     */
    linint1D_data *n0; /**< Density interpolation struct for each species     */
    linint1D_data *t0; /**< Temperature intepolation struct for each species  */
} N0_1D_data;

/**
 * @brief 3D neutral parameters on the target
 */
typedef struct
{
    int n_species;     /**< Number of neutral species                         */
    int *anum;         /**< Neutral species mass number                       */
    int *znum;         /**< Neutral species charge number                     */
    linint3D_data *n0; /**< Density interpolation struct for each species     */
    linint3D_data *t0; /**< Temperature intepolation struct for each species  */
} N0_3D_data;

/**
 * @brief Neutral simulation data
 *
 * The intended usage is that only single type of data is used at a time. This
 * is declared using the `type` field.
 */
typedef struct
{
    N0_1D_data *N01D;  /**< 1D field or NULL if not active           */
    N0_3D_data *N03D;  /**< 3D field or NULL if not active           */
    neutral_type type; /**< Neutral data type wrapped by this struct */
} neutral_data;

/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void neutral_free(neutral_data *data);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void neutral_offload(neutral_data *data);
DECLARE_TARGET_SIMD_UNIFORM(ndata)

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density n0 at the given coordinates.
 *
 * This is a SIMD function.
 *
 * @param n0 pointer where neutral density is stored [m^-3]
 * @param rho normalized poloidal flux coordinate
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param ndata pointer to neutral density data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err neutral_eval_n0(
    real *n0, real rho, real r, real phi, real z, real t, neutral_data *ndata);

/**
 * @brief Evaluate neutral temperature
 *
 * This function evaluates the neutral temperature t0 at the given coordinates.
 *
 * This is a SIMD function.
 *
 * @param t0 pointer where neutral temperature is stored [eV]
 * @param rho normalized poloidal flux coordinate
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param ndata pointer to neutral temperature data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
DECLARE_TARGET_SIMD_UNIFORM(ndata)
a5err neutral_eval_t0(
    real *t0, real rho, real r, real phi, real z, real t, neutral_data *ndata);

/**
 * @brief Get the number of neutral species
 *
 * Retrieve the number of how many neutral species the data contains.
 *
 * This is a SIMD function.
 *
 * @param ndata pointer to neutral data struct
 *
 * @return The number of neutral species
 */
DECLARE_TARGET_SIMD_UNIFORM(ndata)
int neutral_get_n_species(neutral_data *ndata);

#endif
