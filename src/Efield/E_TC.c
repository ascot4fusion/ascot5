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
 * The offload data struct must have the following fields initialized:
 * - E_TC_offload_data.Exyz
 *
 * There is nothing left to initialize.
 *
 * Instead this function only prints values of electric field components as a
 * sanity check.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return Zero to indicate initialization succeeded
 */
int E_TC_init_offload(E_TC_offload_data* offload_data,
                       real** offload_array) {
    // Do no initialization
    offload_data->offload_array_length = 0;
    *offload_array = NULL;

    print_out(VERBOSE_IO, "\nTrivial Cartesian electric field (E_TC)\n");
    print_out(VERBOSE_IO, "E_x = %le, E_y = %le, E_z = %le\n",
              offload_data->Exyz[0], offload_data->Exyz[1],
              offload_data->Exyz[2]);

    return 0;
}

/**
 * @brief Free offload array and return null pointer
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void E_TC_free_offload(E_TC_offload_data* offload_data,
                       real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize electric field simulation data
 *
 * @param Edata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void E_TC_init(E_TC_data* Edata, E_TC_offload_data* offload_data,
               real* offload_array) {
    Edata->Exyz = offload_data->Exyz;
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
