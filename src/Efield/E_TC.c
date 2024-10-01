/**
 * @file E_TC.c
 * @brief Trivial Cartesian Electric field
 *
 * Electric field that has constant x, y, and z components. Note that this field
 * is defined in Cartesian coordinates and not cylindrical. This field is
 * intended for testing purposes and to act as a dummy input.
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "../math.h"
#include "../error.h"
#include "../print.h"
#include "../B_field.h"
#include "E_TC.h"

/**
 * @brief Initialize electric field data and check inputs
 *
 * @param data pointer to the data struct
 * @param exyz electric field vector on cartesian basis [V/m]
 *
 * @return Zero to indicate initialization succeeded
 */
int E_TC_init(E_TC_data* data, real exyz[3]) {

    data->Exyz[0] = exyz[0];
    data->Exyz[1] = exyz[1];
    data->Exyz[2] = exyz[2];

    print_out(VERBOSE_IO, "\nTrivial Cartesian electric field (E_TC)\n");
    print_out(VERBOSE_IO, "E_x = %le, E_y = %le, E_z = %le\n",
              data->Exyz[0], data->Exyz[1], data->Exyz[2]);
    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void E_TC_free(E_TC_data* data) {
    // No resources were allocated
}

/**
 * @brief Evaluate electric field
 *
 * Even though this module represents a Cartesian electric field, the returned
 * values are given in cylindrical coordinates.
 *
 * @param E pointer to array where electric field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Edata pointer to magnetic field data struct
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Zero to indicate success
 */
a5err E_TC_eval_E(real E[3], real r, real phi, real z, E_TC_data* Edata,
                  B_field_data* Bdata) {
    math_vec_xyz2rpz(Edata->Exyz, E, phi);

    return 0;
}
