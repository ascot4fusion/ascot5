/**
 * @file N0_1D.c
 * @brief 1D neutral data with linear interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "N0_1D.h"
#include "../linint/linint.h"

/**
 * @brief Initialize offload data
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload data array
 *
 * @return zero if initialization succeeded
 */
int N0_1D_init_offload(N0_1D_offload_data* offload_data,
                       real** offload_array) {
    int N0_size = offload_data->n_rho;
    int T0_size = offload_data->n_rho;

    offload_data->offload_array_length = offload_data->n_species * N0_size
        + offload_data->n_species * T0_size;

    print_out(VERBOSE_IO, "\n1D neutral density and temperature (N0_1D)\n");
    print_out(VERBOSE_IO, "Grid:  nrho = %4.d   rhomin = %3.3f   rhomax = %3.3f\n",
              offload_data->n_rho,
              offload_data->rho_min, offload_data->rho_max);
    print_out(VERBOSE_IO,
              " Number of neutral species = %d\n",
              offload_data->n_species);
    print_out(VERBOSE_IO,
              "Species Z/A   (Maxwellian)\n");
    for(int i=0; i < offload_data->n_species; i++) {
        print_out(VERBOSE_IO,
                  "      %3d/%3d (%1d)    \n",
                  (int)(offload_data->znum[i]),
                  (int)(offload_data->anum[i]),
                  (int)(offload_data->maxwellian[i]));
    }

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
void N0_1D_free_offload(N0_1D_offload_data* offload_data,
                        real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize neutral data on target
 *
 * This function copies parameters from the offload struct to the struct on
 * target and sets the data pointers on target struct to correct offsets in
 * the offload array.
 *
 * Any initialization that requires any computations must have been done already
 * when the offload struct was initialized.
 *
 * @param ndata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void N0_1D_init(N0_1D_data* ndata, N0_1D_offload_data* offload_data,
                real* offload_array) {
    int N0_size = offload_data->n_rho;
    int T0_size = offload_data->n_rho;
    ndata->n_species  = offload_data->n_species;
    for(int i = 0; i < offload_data->n_species; i++) {
        ndata->anum[i]       = offload_data->anum[i];
        ndata->znum[i]       = offload_data->znum[i];
        ndata->maxwellian[i] = offload_data->maxwellian[i];

        linint1D_init(
            &ndata->n0[i], &offload_array[i * N0_size],
            offload_data->n_rho, NATURALBC,
            offload_data->rho_min, offload_data->rho_max);

        linint1D_init(
            &ndata->t0[i],
            &offload_array[i * T0_size + offload_data->n_species * N0_size],
            offload_data->n_rho, NATURALBC,
            offload_data->rho_min, offload_data->rho_max);
    }
}

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density at the given coordinates using
 * linear interpolation on the 1D neutral density data.
 *
 * @param n0 n0 value will be stored in n0[0]
 * @param rho normalized poloidal flux coordinate
 * @param species neutral species index
 * @param ndata pointer to neutral density data struct
 *
 * @return zero if evaluation succeeded
 */
a5err N0_1D_eval_n0(real* n0, real rho, int species,
                    N0_1D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint1D_eval_f(&n0[0], &ndata->n0[species], rho);

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}

/**
 * @brief Evaluate neutral temperature
 *
 * This function evaluates the neutral temperature at the given coordinates
 * using linear interpolation on the 1D neutral temperature data.
 *
 * @param t0 t0 value will be stored in t0[0]
 * @param rho normalized poloidal flux coordinate
 * @param species neutral species index
 * @param ndata pointer to neutral density data struct
 *
 * @return zero if evaluation succeeded
 */
a5err N0_1D_eval_t0(real* t0, real rho, int species,
                    N0_1D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint1D_eval_f(&t0[0], &ndata->t0[species], rho);

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}
