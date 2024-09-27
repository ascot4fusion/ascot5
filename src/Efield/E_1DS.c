/**
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
#include "../spline/interp.h"
#include "E_1DS.h"

/**
 * @brief Initialize 1DS electric field data
 *
 * @param data pointer to the data struct
 * @param nrho number of points in the rho grid
 * @param rhomin minimum rho value in the grid
 * @param rhomax maximum rho value in the grid
 * @param reff effective minor radius
 * @param dvdrho gradient of the potential
 *
 * @return zero if initialization succeeded
 */
int E_1DS_init(E_1DS_data* data, int nrho, real rhomin, real rhomax, real reff,
               real* dvdrho) {

    int err = 0;
    /* Scale derivatives by effective minor radius */
    real* temp = (real*) malloc( nrho * sizeof(real) );
    for(int i = 0; i < nrho; i++) {
        temp[i] = reff * dvdrho[i];
    }
    err = interp1Dcomp_setup(&data->dV, dvdrho, nrho, NATURALBC,
                             rhomin, rhomax);
    free(temp);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nRadial electric field (E_1DS)");
    print_out(VERBOSE_IO, "(n_rho, rho_min, rho_max) = (%d, %le, %le)\n",
              nrho, rhomin, rhomax);
    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void E_1DS_free_offload(E_1DS_data* data) {
    free(data->dV.c);
}

/**
 * @brief Evaluate 1D spline radial electric field
 *
 * This function evaluates the 1D spline potential gradient of the plasma at the
 * given radial coordinate using linear interpolation, and then calculates the
 * radial electric field by multiplying that with the rho-gradient. Gradient of
 * rho is obtained via magnetic field module.
 *
 * @param E array where the electric field will be stored (E_r -> E[1],
 *        E_phi -> E[1], E_z -> E[2])
 * @param r R-coordiante [m]
 * @param phi phi-coordinate [rad]
 * @param z z-coordiante [m]
 * @param Edata pointer to electric field data
 * @param Bdata pointer to magnetic field data
 *
 * @return zero if evaluation succeeded
 */
a5err E_1DS_eval_E(real E[3], real r, real phi, real z, E_1DS_data* Edata,
                   B_field_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real rho_drho[4];
    err = B_field_eval_rho_drho(rho_drho, r, phi, z, Bdata);
    if(err) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_1DS );
    }
    /* Convert partial derivative to gradient */
    rho_drho[2] = rho_drho[2]/r;
    /* We set the field to zero if outside the profile. */
    if (rho_drho[0] < Edata->dV.x_min || rho_drho[0] > Edata->dV.x_max ) {
        E[0] = 0;
        E[1] = 0;
        E[2] = 0;
        return err;
    }
    real dV;
    interperr += interp1Dcomp_eval_f(&dV, &Edata->dV, rho_drho[0]);

    E[0] = -dV * rho_drho[1];
    E[1] = -dV * rho_drho[2];
    E[2] = -dV * rho_drho[3];

    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_1DS );
    }

    return err;
}
