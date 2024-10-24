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
 * @param ndata pointer to the data struct
 */
int N0_3D_init(N0_3D_data* data,
               int n_r, real r_min, real r_max,
               int n_phi, real phi_min, real phi_max,
               int n_z, real z_min, real z_max,
               int n_species, int* anum, int* znum, int* maxwellian,
               real* density, real* temperature) {
    data->n_species = n_species;
    data->anum = (int*) malloc(n_species * sizeof(int));
    data->znum = (int*) malloc(n_species * sizeof(int));
    data->maxwellian = (int*) malloc(n_species * sizeof(int));
    data->n0 = (linint3D_data*) malloc( n_species * sizeof(linint3D_data) );
    data->t0 = (linint3D_data*) malloc( n_species * sizeof(linint3D_data) );
    for(int i = 0; i < data->n_species; i++) {
        data->anum[i]       = anum[i];
        data->znum[i]       = znum[i];
        data->maxwellian[i] = maxwellian[i];

        real* c = (real*) malloc(n_r * n_phi * n_z * sizeof(real));
        for(int i = 0; i < n_r * n_phi * n_z; i++) {
            c[i] = density[i];
        }
        linint3D_init(&data->n0[i], c, n_r, n_phi, n_z,
                      NATURALBC, PERIODICBC, NATURALBC,
                      r_min, r_max, phi_min, phi_max, z_min, z_max);
        c = (real*) malloc(n_r * n_phi * n_z * sizeof(real));
        for(int i = 0; i < n_r * n_phi * n_z; i++) {
            c[i] = temperature[i];
        }
        linint3D_init(&data->t0[i], c, n_r, n_phi, n_z,
                      NATURALBC, PERIODICBC, NATURALBC,
                      r_min, r_max, phi_min, phi_max, z_min, z_max);
    }
    print_out(VERBOSE_IO, "\n3D neutral density and temperature (N0_3D)\n");
    print_out(VERBOSE_IO, "Grid:  nR = %4.d   Rmin = %3.3f   Rmax = %3.3f\n",
              n_r,
              r_min, r_max);
    print_out(VERBOSE_IO, "       nz = %4.d   zmin = %3.3f   zmax = %3.3f\n",
              n_z,
              z_min, z_max);
    print_out(VERBOSE_IO, "     nphi = %4.d phimin = %3.3f phimax = %3.3f\n",
              n_phi,
              phi_min, phi_max);
    print_out(VERBOSE_IO,
              " Number of neutral species = %d\n",
              data->n_species);
    print_out(VERBOSE_IO, "Species Z/A   (Maxwellian)\n");
    for(int i=0; i < data->n_species; i++) {
        print_out(VERBOSE_IO,
                  "      %3d/%3d (%1d)    \n",
                  data->znum[i], data->anum[i], data->maxwellian[i]);
    }

    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param offload_data pointer to the data struct
 */
void N0_3D_free(N0_3D_data* data) {
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
void N0_3D_offload(N0_3D_data* data) {
    //TODO: Implement
}

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density at the given coordinates using
 * trilinear interpolation on the 3D neutral density data.
 *
 * @param n0 n0 value will be stored in n0[0]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param ndata pointer to neutral data struct
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
 * @param t0 t0 value will be stored in t0[0]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param ndata pointer to neutral data struct
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
