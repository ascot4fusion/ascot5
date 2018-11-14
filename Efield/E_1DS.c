/**
 * @author Joona Kontula joona.kontula@aalto.fi
 * @file E_1DS.c
 * @brief 1D spline electric field evaluation functions
 *
 * Radial electric field evaluated from 1D profile. Exactly as E_1D but the
 * potential is evaluated with splines.
 */
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "../B_field.h"
#include "E_1DS.h"

void E_1DS_init_offload(E_1DS_offload_data* offload_data,
                        real** offload_array) {
    // redundant function
}

/**
 * @brief Load 1DS electric field data and prepare parameters
 *
 * This function reads the 1D electric field data from input.erad, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @todo move spline initialization here
 */
int E_1DS_init_offload(E_1DS_offload_data* offload_data, real** offload_array) {
    print_out(VERBOSE_IO, "\nRadial electric field (E_1DS)");
    print_out(VERBOSE, "(n_rho, rho_min, rho_max) = (%d, %le, %le)\n",
              offload_data->n_rho,offload_data->rho_min,
              offload_data->rho_max);
    return 0;
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void E_1DS_free_offload(E_1DS_offload_data* offload_data,
                       real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize 1D spline electric field data struct on target
 *
 * This function copies the 1D spline electric field parameters from the offload
 * struct to the struct on target and sets the 1D spline electric field data
 * pointers to correct offsets in the offload array.
 *
 * @todo Move spline initialization to offload
 *
 * @param Edata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void E_1DS_init(E_1DS_data* Edata, E_1DS_offload_data* offload_data,
               real* offload_array) {

    Edata->n_rho = offload_data->n_rho;
    Edata->rho_min = offload_data->rho_min;
    Edata->rho_max = offload_data->rho_max;
    Edata->rho_grid = (Edata->rho_max - Edata->rho_min)/(Edata->n_rho - 1);
    interp1Dcomp_init(&Edata->dV, &offload_array[Edata->n_rho],
                      Edata->n_rho, Edata->rho_min, Edata->rho_max,
                      Edata->rho_grid);
}

/**
 * @brief Evaluate 1D spline radial electric field
 *
 * This function evaluates the 1D spline potential gradient of the plasma at the
 * given radial coordinate using linear interpolation, and then calculates the
 * radial electric field by multiplying that with the rho-gradient.
 *
 * @param E array where the electric field will be stored (E_r -> E[1],
 *        E_phi -> E[1], E_z -> E[2])
 * @param rho_drho array where rho and components of gradrho a are stored
 *          (rho -> rho_drho[0], [gradrho]_r -> rho_drho[1],
 *          [gradrho]_phi -> rho_drho[2],
 *          [gradrho]_z -> rho_drho[3])
 * @param Edata pointer to electric field data
 */
a5err E_1DS_eval_E(real E[], real r, real phi, real z, E_1DS_data* Edata,
                   B_field_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real rho_drho[4];
    B_field_eval_rho_drho(rho_drho, r, phi, z, Bdata);
    /* Convert partial derivative to gradient */
    rho_drho[2] = rho_drho[2]/r;
    /* We set the field to zero if outside the profile. */
    if (rho_drho[0] < Edata->rho_min || rho_drho[0] > Edata->rho_max ) {
        E[0] = 0;
        E[1] = 0;
        E[2] = 0;
        return err;
    }
    real dV;
    interperr += interp1Dcomp_eval_B(&dV, &Edata->dV, rho_drho[0]);

    E[0] = dV * rho_drho[1];
    E[1] = dV * rho_drho[2];
    E[2] = dV * rho_drho[3];

    if(interperr) {err = error_raise( ERR_OUTSIDE_PLASMA, __LINE__ );}

    return err;
}
