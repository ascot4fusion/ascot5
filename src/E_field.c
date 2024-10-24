/**
 * @file E_field.c
 * @brief Electric field interface
 *
 * This is an interface through which electric field data is initialized and
 * accessed. Reading e.g. from disk is done elsewhere.
 *
 * To add a new electric field instance, make sure these functions are
 * implemented and called from this interface, and that E_field.h contains enum
 * type for the new instance.
 *
 * The interface checks which instance given data corresponds to from
 * E_field_offload_data.type and E_field_data.type from the struct that is given
 * as an argument, and calls the relevant function for that instance.
 */
#include <stdio.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "E_field.h"
#include "B_field.h"
#include "Efield/E_TC.h"
#include "Efield/E_1DS.h"

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void E_field_free(E_field_data* data) {
    switch(data->type) {
        case E_field_type_1DS:
            E_1DS_free(&data->E1DS);
            break;
        case E_field_type_TC:
            E_TC_free(&data->ETC);
            break;
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void E_field_offload(E_field_data* data) {
    switch(data->type) {
        case E_field_type_1DS:
            E_1DS_offload(&data->E1DS);
            break;
        case E_field_type_TC:
            E_TC_offload(&data->ETC);
            break;
    }
}

/**
 * @brief Evaluate electric field
 *
 * This function evaluates the electric field at the given coordinates. Note
 * that magnetic field data is also required in case electric field is e.g.
 * a flux quantity.
 *
 * The values are stored in the given array as:
 *
 * - E[0] = ER
 * - E[1] = Ephi
 * - E[2] = Ez
 *
 * This is a SIMD function.
 *
 * @param E pointer to array where electric field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param Edata pointer to electric field data struct
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err E_field_eval_E(real E[3], real r, real phi, real z, real t,
                     E_field_data* Edata, B_field_data* Bdata) {
    a5err err = 0;

    switch(Edata->type) {

        case E_field_type_1DS:
            err = E_1DS_eval_E(E, r, phi, z, &(Edata->E1DS), Bdata);
            break;

        case E_field_type_TC:
            err = E_TC_eval_E(E, r, phi, z, &(Edata->ETC), Bdata);
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_E_FIELD );
            break;
    }

    return err;
}
