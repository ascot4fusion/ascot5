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
#include "efield.h"
#include "defines.h"
#include "bfield.h"
#include "efield_cartesian.h"
#include "efield_potential1d.h"
#include "errors.h"
#include <stdio.h>

/**
 * Free allocated resources
 *
 * @param data pointer to the data struct
 */
void E_field_free(Efield *efield)
{
    switch (efield->type)
    {
    case E_field_type_potential1d:
        EfieldPotential1D_free(efield->potential1d);
        break;
    case E_field_type_cartesian:
        EfieldCartesian_free(efield->cartesian);
        break;
    }
}

/**
 * Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void E_field_offload(Efield *efield)
{
    switch (efield->type)
    {
    case E_field_type_potential1d:
        EfieldPotential1D_offload(efield->potential1d);
        break;
    case E_field_type_cartesian:
        EfieldCartesian_offload(efield->cartesian);
        break;
    }
}

/**
 * Evaluate electric field
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
 * @param efield pointer to electric field data struct
 * @param bfield pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err E_field_eval_E(
    real e[3], real r, real phi, real z, real t, Efield *efield, Bfield *bfield)
{
    a5err err = 0;
    (void)t; // Unused until dynamic efield is implemented

    switch (efield->type)
    {

    case E_field_type_potential1d:
        err =
            EfieldPotential1D_eval_e(e, r, phi, z, efield->potential1d, bfield);
        break;

    case E_field_type_cartesian:
        err = EfieldCartesian_eval_e(e, phi, efield->cartesian);
        break;

    default:
        /* Unregonized input. Produce error. */
        err = error_raise(ERR_UNKNOWN_INPUT, __LINE__, EF_E_FIELD);
        break;
    }

    return err;
}
