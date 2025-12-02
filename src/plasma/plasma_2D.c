/**
 * @file plasma_2D.c
 * @brief 2D linearly interpolated plasma
 *
 * Plasma data which is defined in a (R,z) uniform grid from which the values
 * are interpolated linearly.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../error.h"
#include "../consts.h"
#include "../print.h"
#include "../linint/linint.h"
#include "plasma_2D.h"

/**
 * @brief Initialize 2D plasma data and check inputs
 *
 * @param data pointer to the data struct
 *
 * @return zero if initialization succes
 */
int plasma_2D_init(plasma_2D_data* data, int nr, int nz, int nion,
                   real r_min, real r_max, real z_min, real z_max,
                   int* anum, int* znum, real* mass, real* charge,
                   real* Te, real* Ti, real* ne, real* ni, real* vtor) {

    data->n_species = nion + 1;
    data->anum = (int*) malloc( nion*sizeof(int) );
    data->znum = (int*) malloc( nion*sizeof(int) );
    data->mass = (real*) malloc( (nion+1)*sizeof(real) );
    data->charge = (real*) malloc( (nion+1)*sizeof(real) );
    for(int i = 0; i < data->n_species; i++) {
        if(i < nion) {
            data->znum[i] = znum[i];
            data->anum[i] = anum[i];
        }
        data->mass[i] = mass[i];
        data->charge[i] = charge[i];
    }

    data->dens = (linint2D_data*) malloc( (nion + 1) * sizeof(linint2D_data) );
    data->temp = (linint2D_data*) malloc( 2 * sizeof(linint2D_data) );
    data->vtor = (linint2D_data*) malloc( sizeof(linint2D_data) );
    real* c = (real*) malloc(nr * nz * sizeof(real));
    for(int j = 0; j < nr * nz; j++) {
        c[j] = ne[j];
    }
    linint2D_init(
        &data->dens[0], c, nr, nz, NATURALBC, NATURALBC, r_min, r_max,
        z_min, z_max);

    c = (real*) malloc(nr * nz * sizeof(real));
    for(int j = 0; j < nr * nz; j++) {
        c[j] = Te[j];
    }
    linint2D_init(
        &data->temp[0], c, nr, nz, NATURALBC, NATURALBC, r_min, r_max,
        z_min, z_max);

    c = (real*) malloc(nr * nz * sizeof(real));
    for(int j = 0; j < nr * nz; j++) {
        c[j] = Ti[j];
    }
    linint2D_init(
        &data->temp[1], c, nr, nz, NATURALBC, NATURALBC, r_min, r_max,
        z_min, z_max);

    c = (real*) malloc(nr * nz * sizeof(real));
    for(int j = 0; j < nr * nz; j++) {
        c[j] = vtor[j];
    }
    linint2D_init(
        &data->vtor[0], c, nr, nz, NATURALBC, NATURALBC, r_min, r_max,
        z_min, z_max);

    for(int i = 0; i < nion; i++) {

        c = (real*) malloc(nr * nz * sizeof(real));
        for(int j = 0; j < nr * nz; j++) {
            c[j] = ni[j + i * nion];
        }
        linint2D_init(
            &data->dens[i+1], c, nr, nz, NATURALBC, NATURALBC, r_min, r_max,
            z_min, z_max);
    }

    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void plasma_2D_free(plasma_2D_data* data) {
    free(data->anum);
    free(data->znum);
    free(data->mass);
    free(data->charge);
    for(int i = 0; i < data->n_species; i++) {
        free(data->dens[i].c);
    }
    free(data->temp[0].c);
    free(data->temp[1].c);
    free(data->vtor[0].c);
    free(data->dens);
    free(data->temp);
    free(data->vtor);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void plasma_2D_offload(plasma_2D_data* data) {
    //TODO: Implement
}

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param temp pointer to where evaluated temperature [J] is stored
 * @param r R coordinate
 * @param z z coordinate
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_2D_eval_temp(real* temp, real r, real z, int species,
                          plasma_2D_data* pls_data) {
    a5err err = 0;
    int i = species > 0;
    int interperr = linint2D_eval_f(temp, &pls_data->temp[i], r, z);
    if(interperr) {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_2D);
    }
    return err;
}

/**
 * @brief Evaluate plasma density
 *
 * This function evaluates the density of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param dens pointer to where evaluated density [m^-3] is stored
 * @param r R coordinate
 * @param z z coordinate
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_2D_eval_dens(real* dens, real r, real z, int species,
                          plasma_2D_data* pls_data) {
    a5err err = 0;
    int interperr = linint2D_eval_f(dens, &pls_data->dens[species], r, z);
    if(interperr) {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_2D);
    }
    return err;
}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 * This function evaluates the density and temperature of all plasma species at
 * the given radial coordinate using linear interpolation.
 *
 * @param dens pointer to where interpolated densities [m^-3] are stored
 * @param temp pointer to where interpolated temperatures [J] are stored
 * @param r R coordinate
 * @param z z coordinate
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_2D_eval_densandtemp(real* dens, real* temp, real r, real z,
                                 plasma_2D_data* pls_data) {
    int interperr = 0;
    for(int i = 0; i < pls_data->n_species; i++) {
        int ision = i > 0;
        interperr += linint2D_eval_f(&temp[i], &pls_data->temp[ision], r, z);
        interperr += linint2D_eval_f(&dens[i], &pls_data->dens[i], r, z);
    }
    a5err err = 0;
    if(interperr) {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_2D);
    }
    return err;
}

/**
 * @brief Evalate plasma flow along the field lines
 *
 * @param vflow pointer where the flow value is stored [m/s]
 * @param r R coordinate
 * @param z z coordinate
 * @param pls_data pointer to plasma data
 */
a5err plasma_2D_eval_flow(real* vflow, real r, real z,
                          plasma_2D_data* pls_data) {
    a5err err = 0;
    int interperr = linint2D_eval_f(vflow, pls_data->vtor, r, z);
    if(interperr) {
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_2D);
    }
    *vflow *= r;
    return err;
}