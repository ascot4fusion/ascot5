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
 * @brief Initialize offload data
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload data array
 *
 * @return zero if initialization succeeded
 */
int N0_3D_init_offload(N0_3D_offload_data* offload_data,
                       real** offload_array) {
    int N0_size = offload_data->n_r * offload_data->n_phi * offload_data->n_z;
    int T0_size = offload_data->n_r * offload_data->n_phi * offload_data->n_z;

    offload_data->offload_array_length = offload_data->n_species * N0_size
        + offload_data->n_species * T0_size;

    print_out(VERBOSE_IO, "\n3D neutral density and temperature (N0_3D)\n");
    print_out(VERBOSE_IO, "Grid:  nR = %4.d   Rmin = %3.3f   Rmax = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max);
    print_out(VERBOSE_IO, "       nz = %4.d   zmin = %3.3f   zmax = %3.3f\n",
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max);
    print_out(VERBOSE_IO, "     nphi = %4.d phimin = %3.3f phimax = %3.3f\n",
              offload_data->n_phi,
              offload_data->phi_min, offload_data->phi_max);
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
void N0_3D_free_offload(N0_3D_offload_data* offload_data,
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
void N0_3D_init(N0_3D_data* ndata, N0_3D_offload_data* offload_data,
                real* offload_array) {
    int N0_size = offload_data->n_r * offload_data->n_phi * offload_data->n_z;
    int T0_size = offload_data->n_r * offload_data->n_phi * offload_data->n_z;
    ndata->n_species  = offload_data->n_species;
    for(int i = 0; i < offload_data->n_species; i++) {
        ndata->anum[i]       = offload_data->anum[i];
        ndata->znum[i]       = offload_data->znum[i];
        ndata->maxwellian[i] = offload_data->maxwellian[i];

        linint3D_init(
            &ndata->n0[i], &offload_array[i * N0_size],
            offload_data->n_r, offload_data->n_phi, offload_data->n_z,
            NATURALBC, PERIODICBC, NATURALBC,
            offload_data->r_min, offload_data->r_max,
            offload_data->phi_min, offload_data->phi_max,
            offload_data->z_min, offload_data->z_max);

        linint3D_init(
            &ndata->t0[i],
            &offload_array[i * T0_size + offload_data->n_species * N0_size],
            offload_data->n_r, offload_data->n_phi, offload_data->n_z,
            NATURALBC, PERIODICBC, NATURALBC,
            offload_data->r_min, offload_data->r_max,
            offload_data->phi_min, offload_data->phi_max,
            offload_data->z_min, offload_data->z_max);
    }
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
