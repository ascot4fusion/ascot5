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

#include "defines.h"
#include "neutral_data.h"
#include "parallel.h"

/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void Neutral_free(Neutral *neutral);

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void Neutral_offload(Neutral *neutral);

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
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
DECLARE_TARGET_SIMD_UNIFORM(neutral)
err_t Neutral_eval_n0(
    real *n0, real rho, real r, real phi, real z, real t, Neutral *neutral);

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
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
DECLARE_TARGET_SIMD_UNIFORM(neutral)
err_t Neutral_eval_T0(
    real *T0, real rho, real r, real phi, real z, real t, Neutral *neutral);

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
DECLARE_TARGET_SIMD_UNIFORM(neutral)
int neutral_get_n_species(Neutral *neutral);

#endif
