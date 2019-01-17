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

/**
 * @brief Initialize 1DS electric field data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct and to construct
 * a cubic spline whose coefficients are stored in re-allocated offload array.
 *
 * The offload data struct must have the following fields initialized:
 * - E_1DS_offload_data.n_rho
 * - E_1DS_offload_data.rho_min
 * - E_1DS_offload_data.rho_max
 *
 * E_1DS_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[i] = dV_drho(rho_i)   [V]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int E_1DS_init_offload(E_1DS_offload_data* offload_data, real** offload_array) {
    /* Fill rest of the offload data struct */
    offload_data->rho_grid = (offload_data->rho_max - offload_data->rho_min)
        / (offload_data->n_rho - 1);

    /* Spline initialization. Use spline structs for temporary storage */
    int err = 0;
    int splinesize = 2 * offload_data->n_rho;
    interp1D_data dV_drho;
    err += interp1Dcomp_init(&dV_drho, *offload_array,
                             offload_data->n_rho, offload_data->rho_min,
                             offload_data->rho_max, offload_data->rho_grid);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Re-allocate the offload array and store the spline coefficients there */
    free(*offload_array);

    offload_data->offload_array_length = splinesize;
    *offload_array =
        (real*) malloc( offload_data->offload_array_length * sizeof(real) );

    for(int i = 0; i < splinesize; i++) {
        (*offload_array)[i] = dV_drho.c[i];
    }

    interp1Dcomp_free(&dV_drho);

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nRadial electric field (E_1DS)");
    print_out(VERBOSE_IO, "(n_rho, rho_min, rho_max) = (%d, %le, %le)\n",
              offload_data->n_rho,offload_data->rho_min,
              offload_data->rho_max);
    return err;
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
    Edata->dV.n_r    = offload_data->n_rho;
    Edata->dV.r_min  = offload_data->rho_min;
    Edata->dV.r_max  = offload_data->rho_max;
    Edata->dV.r_grid = offload_data->rho_grid;
    Edata->dV.c      = offload_array;
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
 *
 * @return zero if evaluation succeeded
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
    if (rho_drho[0] < Edata->dV.r_min || rho_drho[0] > Edata->dV.r_max ) {
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

    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_1DS );}

    return err;
}
