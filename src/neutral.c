/**
 * @file neutral.c
 * @brief Neutral interface
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
#include <stdio.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "neutral.h"
#include "neutral/N0_1D.h"
#include "neutral/N0_3D.h"

/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void neutral_free(neutral_data* data) {
    switch(data->type) {
        case neutral_type_1D:
            N0_1D_free(&data->N01D);
            break;
        case neutral_type_3D:
            N0_3D_free(&data->N03D);
            break;
    }
}

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
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param ndata pointer to neutral density data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err neutral_eval_n0(real* n0, real rho, real r, real phi, real z, real t,
                      neutral_data* ndata) {
    a5err err = 0;

    switch(ndata->type) {
        case neutral_type_1D:
            err = N0_1D_eval_n0(n0, rho, &(ndata->N01D));
            break;
        case neutral_type_3D:
            err = N0_3D_eval_n0(n0, r, phi, z, &(ndata->N03D));
            break;
        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_NEUTRAL);
            break;
    }

    if(err) {
        /* Return some reasonable values to avoid further errors */
        n0[0] = 0;
    }

    return err;
}

/**
 * @brief Evaluate neutral temperature
 *
 * This function evaluates the neutral temperature t0 at the given coordinates.
 *
 * This is a SIMD function.
 *
 * @param t0 pointer where neutral temperature is stored [J]
 * @param rho normalized poloidal flux coordinate
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param ndata pointer to neutral temperature data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err neutral_eval_t0(real* t0, real rho, real r, real phi, real z, real t,
                      neutral_data* ndata) {
    a5err err = 0;

    switch(ndata->type) {
        case neutral_type_1D:
            err = N0_1D_eval_t0(t0, rho, &(ndata->N01D));
            break;
        case neutral_type_3D:
            err = N0_3D_eval_t0(t0, r, phi, z, &(ndata->N03D));
            break;
        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_NEUTRAL);
            break;
    }

    if(err) {
        /* Return some reasonable values to avoid further errors */
        t0[0] = 1;
    }

    return err;
}

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
int neutral_get_n_species(neutral_data* ndata) {
    int n = 0;
    switch(ndata->type) {
        case neutral_type_1D:
            n = N0_1D_get_n_species(&(ndata->N01D));
            break;
        case neutral_type_3D:
            n = N0_3D_get_n_species(&(ndata->N03D));
            break;
    }

    return n;
}
