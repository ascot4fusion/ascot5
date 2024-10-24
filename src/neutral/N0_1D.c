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
 * @brief Initialize data
 *
 * @param data pointer to data struct
 * @param n_rho number of r grid points in the data
 * @param rho_min minimum r coordinate in the grid in the data [m]
 * @param rho_max maximum r coordinate in the grid in the data [m]
 * @param n_species number of neutral species
 * @param anum neutral species mass number
 * @param znum neutral species charge number
 * @param maxwellian is the species distribution Maxwellian or monoenergetic
 *
 * @return zero if initialization succeeded
 */
int N0_1D_init(N0_1D_data* data, int n_rho, real rho_min, real rho_max,
               int n_species, int* anum, int* znum, int* maxwellian,
               real* density, real* temperature) {

    data->n_species = n_species;
    data->anum = (int*) malloc(n_species * sizeof(int));
    data->znum = (int*) malloc(n_species * sizeof(int));
    data->maxwellian = (int*) malloc(n_species * sizeof(int));
    data->n0 = (linint1D_data*) malloc( n_species * sizeof(linint1D_data) );
    data->t0 = (linint1D_data*) malloc( n_species * sizeof(linint1D_data) );
    for(int i = 0; i < data->n_species; i++) {
        data->anum[i] = anum[i];
        data->znum[i] = znum[i];
        data->maxwellian[i] = maxwellian[i];

        real* c = (real*) malloc(n_rho * sizeof(real));
        for(int i = 0; i < n_rho; i++) {
            c[i] = density[i];
        }
        linint1D_init(&data->n0[i], c, n_rho, NATURALBC, rho_min, rho_max);
        c = (real*) malloc(n_rho * sizeof(real));
        for(int i = 0; i < n_rho; i++) {
            c[i] = temperature[i];
        }
        linint1D_init(&data->t0[i], c, n_rho, NATURALBC, rho_min, rho_max);
    }

    print_out(VERBOSE_IO, "\n1D neutral density and temperature (N0_1D)\n");
    print_out(VERBOSE_IO,
              "Grid:  nrho = %4.d   rhomin = %3.3f   rhomax = %3.3f\n",
              n_rho, rho_min, rho_max);
    print_out(VERBOSE_IO, " Number of neutral species = %d\n", data->n_species);
    print_out(VERBOSE_IO, "Species Z/A   (Maxwellian)\n");
    for(int i=0; i < data->n_species; i++) {
        print_out(VERBOSE_IO,
                  "      %3d/%3d (%1d)    \n",
                  (int)(data->znum[i]),
                  (int)(data->anum[i]),
                  (int)(data->maxwellian[i]));
    }

    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void N0_1D_free(N0_1D_data* data) {
    free(data->anum);
    free(data->znum);
    free(data->maxwellian);
    for(int i = 0; i < data->n_species; i++) {
        free(data->n0->c);
        free(data->t0->c);
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void N0_1D_offload(N0_1D_data* data) {
    //TODO: Implement
}

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density at the given coordinates using
 * linear interpolation on the 1D neutral density data.
 *
 * @param n0 n0 value will be stored in n0[0]
 * @param rho normalized poloidal flux coordinate
 * @param ndata pointer to neutral data struct
 *
 * @return zero if evaluation succeeded
 */
a5err N0_1D_eval_n0(real* n0, real rho, N0_1D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for(int i=0; i<ndata->n_species; i++) {
        interperr += linint1D_eval_f(&n0[i], &ndata->n0[i], rho);
    }

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
 * @param ndata pointer to neutral data struct
 *
 * @return zero if evaluation succeeded
 */
a5err N0_1D_eval_t0(real* t0, real rho, N0_1D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for(int i=0; i<ndata->n_species; i++) {
        interperr += linint1D_eval_f(&t0[i], &ndata->t0[i], rho);
    }

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}

/**
 * @brief Return number of neutral species
 *
 * @param ndata pointer to neutral data struct
 *
 * @return number of neutral species
 */
int N0_1D_get_n_species(N0_1D_data* ndata) {
    return ndata->n_species;
}
