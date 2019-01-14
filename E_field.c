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
#include "Efield/E_1D.h"
#include "Efield/E_1DS.h"
#include "Efield/E_3D.h"


/**
 * @brief Load electric field data and prepare parameters
 *
 * This function fills the relevant electric field offload struct with
 * parameters and allocates and fills the offload array.
 *
 * The offload data has to have a type when this function is called as it should
 * be set when the offload data is constructed from inputs.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded
 */
int E_field_init_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
    int err = 0;

    switch(offload_data->type) {

        case E_field_type_1D:
            err = E_1D_init_offload(&(offload_data->E1D), offload_array);
            offload_data->offload_array_length =
                offload_data->E1D.offload_array_length;
            break;

        case E_field_type_1DS:
            err = E_1DS_init_offload(&(offload_data->E1DS), offload_array);
            offload_data->offload_array_length =
                offload_data->E1DS.offload_array_length;
            break;

        case E_field_type_TC:
            err = E_TC_init_offload(&(offload_data->ETC), offload_array);
            offload_data->offload_array_length =
                offload_data->ETC.offload_array_length;
            break;

       case E_field_type_3D:
	    err = E_3D_init_offload(&(offload_data->E3D), offload_array);
	    offload_data->offload_array_length = 
	        offload_data->E3D.offload_array_length;
	    break;

        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized electric field type.");
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
void E_field_free_offload(E_field_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {

        case E_field_type_1D:
            E_1D_free_offload(&(offload_data->E1D), offload_array);
            break;

        case E_field_type_1DS:
            E_1DS_free_offload(&(offload_data->E1DS), offload_array);
            break;

        case E_field_type_TC:
            E_TC_free_offload(&(offload_data->ETC), offload_array);
            break;

        case E_field_type_3D:
            E_3D_free_offload(&(offload_data->E3D), offload_array);
	    break;

    }
}

/**
 * @brief Initialize electric field data struct on target
 *
 * This function copies the electric field parameters from the offload struct
 * to the struct on target and sets the electric field data pointers to correct
 * offsets in the offload array.
 *
 * This function returns error if the offload data has not been initialized.
 * The instances themselves should not return an error since all they do is
 * assign pointers and values.
 *
 * @param Edata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array the offload array
 *
 * @return Non-zero integer if offload was not initialized beforehand
 */
int E_field_init(E_field_data* Edata, E_field_offload_data* offload_data,
                  real* offload_array) {
    int err = 0;

    switch(offload_data->type) {

        case E_field_type_1D:
            E_1D_init(&(Edata->E1D), &(offload_data->E1D), offload_array);
            break;

        case E_field_type_1DS:
            E_1DS_init(&(Edata->E1DS), &(offload_data->E1DS), offload_array);
            break;

        case E_field_type_TC:
            E_TC_init(&(Edata->ETC), &(offload_data->ETC), offload_array);
            break;

        case E_field_type_3D:
            E_3D_init(&(Edata->E3D), &(offload_data->E3D), offload_array);
	    break;

        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized electric field type.\n");
            err = 1;
            break;
    }
    Edata->type = offload_data->type;

    return err;
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
 * @param Edata pointer to electric field data struct
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err E_field_eval_E(real E[3], real r, real phi, real z, E_field_data* Edata,
                     B_field_data* Bdata) {
    a5err err = 0;

    switch(Edata->type) {

        case E_field_type_1D:
            err = E_1D_eval_E(E, r, phi, z, &(Edata->E1D), Bdata);
            break;

        case E_field_type_1DS:
            err = E_1DS_eval_E(E, r, phi, z, &(Edata->E1DS), Bdata);
            break;

        case E_field_type_TC:
            err = E_TC_eval_E(E, r, phi, z, &(Edata->ETC), Bdata);
            break;

        case E_field_type_3D:
	    err = E_3D_eval_E(E, r, phi, z, &(Edata->E3D), Bdata);
	    break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_E_FIELD );
            break;
    }

    return err;
}
