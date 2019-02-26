/**
 * @file neutral.c
 * @brief Neutral density interface
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
#include "neutral/N0_3D.h"

/**
 * @brief Load neutral data and prepare parameters
 *
 * This function fills the relevant neutral offload struct with parameters and
 * allocates and fills the offload array. Sets offload array length in the
 * offload struct.
 *
 * The offload data has to have a type when this function is called as it should
 * be set when the offload data is constructed from inputs.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return Zero if initialization succeeded
 */
int neutral_init_offload(neutral_offload_data* offload_data,
                         real** offload_array) {
    int err = 0;

    switch(offload_data->type) {
        case neutral_type_3D:
            err = N0_3D_init_offload(&(offload_data->N03D), offload_array);
            offload_data->offload_array_length =
                offload_data->N03D.offload_array_length;
            break;
        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized neutral data type.");
            err = 1;
            break;
    }
    if(!err) {
        print_out(VERBOSE_IO, "Estimated memory usage %.1f MB\n",
                  offload_data->offload_array_length
                  * sizeof(real) / (1024.0*1024.0) );
    }

    return err;
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void neutral_free_offload(neutral_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {
        case neutral_type_3D:
            N0_3D_free_offload(&(offload_data->N03D), offload_array);
            break;
    }
}

/**
 * @brief Initialize neutral data struct on target
 *
 * This function copies the neutral data parameters from the offload struct to
 * the struct on target and sets the neutral data pointers to correct offsets in
 * the offload array.
 *
 * @param ndata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @return Zero if initialization succeeded
 */
int neutral_init(neutral_data* ndata, neutral_offload_data* offload_data,
                  real* offload_array) {
    int err = 0;

    switch(offload_data->type) {
        case neutral_type_3D:
            N0_3D_init(&(ndata->N03D),
                       &(offload_data->N03D), offload_array);
            break;
        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized neutral data type.\n");
            err = 1;
            break;
    }
    ndata->type = offload_data->type;

    return err;
}

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density n0 at the given coordinates.
 *
 * This is a SIMD function.
 *
 * @param n0 pointer where neutral density is stored [m^-3]
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param ndata pointer to neutral density data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err neutral_eval_n0(real* n0, real r, real phi, real z, real t,
                      neutral_data* ndata) {
    a5err err = 0;

    switch(ndata->type) {
        case neutral_type_3D:
            err = N0_3D_eval_n0(n0, r, phi, z, 0, &(ndata->N03D));
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
