/**
 * @file N0_3D.c
 * @brief 3D neutral data with trilinear interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "N0_3D.h"
#include "../linint/linint.h"

/**
 * @brief Initialize neutral data
 *
 * @param data pointer to data struct
 * @param n_r number of r grid points in the data
 * @param r_min minimum r coordinate in the grid in the data [m]
 * @param r_max maximum r coordinate in the grid in the data [m]
 * @param n_phi number of phi grid points in the data
 * @param phi_min minimum phi coordinate in the grid in the data [rad]
 * @param phi_max maximum phi coordinate in the grid in the data [rad]
 * @param n_z number of z grid points in the data
 * @param z_min minimum z coordinate in the grid in the data [m]
 * @param z_max maximum z coordinate in the grid in the data [m]
 * @param n_species number of neutral species
 * @param anum neutral species mass number
 * @param znum neutral species charge number
 * @param density neutral species-wise density [m^-3]
 * @param temperature neutral species-wise temperature [J]
 */
int N0_3D_init(N0_3D_data* data,
               int n_r, real r_min, real r_max,
               int n_phi, real phi_min, real phi_max,
               int n_z, real z_min, real z_max,
               int n_species, int* anum, int* znum,
               real* density, real* temperature) {
    data->n_species = n_species;
    data->anum = (int*) malloc(n_species * sizeof(int));
    data->znum = (int*) malloc(n_species * sizeof(int));
    data->n0 = (linint3D_data*) malloc( n_species * sizeof(linint3D_data) );
    data->t0 = (linint3D_data*) malloc( n_species * sizeof(linint3D_data) );
    int err = 0;
    for(int i = 0; i < data->n_species; i++) {
        data->anum[i]       = anum[i];
        data->znum[i]       = znum[i];

        real* c = (real*) malloc(n_r * n_phi * n_z * sizeof(real));
        err += c == NULL ? 1 : 0;
        for(int i = 0; i < n_r * n_phi * n_z; i++) {
            c[i] = density[i];
        }
        linint3D_init(&data->n0[i], c, n_r, n_phi, n_z,
                             NATURALBC, PERIODICBC, NATURALBC,
                             r_min, r_max, phi_min, phi_max, z_min, z_max);
        c = (real*) malloc(n_r * n_phi * n_z * sizeof(real));
        err += c == NULL ? 1 : 0;
        for(int i = 0; i < n_r * n_phi * n_z; i++) {
            c[i] = temperature[i];
        }
        linint3D_init(&data->t0[i], c, n_r, n_phi, n_z,
                      NATURALBC, PERIODICBC, NATURALBC,
                      r_min, r_max, phi_min, phi_max, z_min, z_max);
    }
    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void N0_3D_free(N0_3D_data* data) {
    free(data->anum);
    free(data->znum);
    for(int i = 0; i < data->n_species; i++) {
        free(data->n0[i].c);
        free(data->t0[i].c);
    }
    free(data->n0);
    free(data->t0);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void N0_3D_offload(N0_3D_data* data) {
    //TODO: Implement
}

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density at the given coordinates using
 * trilinear interpolation on the 3D neutral density data.
 *
 * @param n0 pointer where neutral density is stored [m^-3]
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param ndata pointer to neutral density data struct
 *
 * @return zero if evaluation succeeded
 */
a5err N0_3D_eval_n0(real* n0, real r, real phi, real z, N0_3D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for(int i=0; i<ndata->n_species; i++) {
        interperr += linint3D_eval_f(&n0[i], &ndata->n0[i], r, phi, z);
    }

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_3D);
    }

    return err;
}

/**
 * @brief Evaluate neutral temperature
 *
 * This function evaluates the neutral temperature at the given coordinates
 * using trilinear interpolation on the 3D neutral temperature data.
 *
 * @param t0 pointer where neutral temperature is stored [J]
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 *
 * @return zero if evaluation succeeded
 */
a5err N0_3D_eval_t0(real* t0, real r, real phi, real z, N0_3D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for(int i=0; i<ndata->n_species; i++) {
        interperr += linint3D_eval_f(&t0[i], &ndata->t0[i], r, phi, z);
    }

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_3D);
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
int N0_3D_get_n_species(N0_3D_data* ndata) {
    return ndata->n_species;
}
