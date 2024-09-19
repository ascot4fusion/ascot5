/**
 * @file E_3DS.c
 * @brief 3D electric field with tricubic interpolation
 *
 * This module represents an electric field where data is given in \f$R\phi z\f$
 * -grid from which it is interpolated with tricubic splines.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown. For \f$\phi\f$-grid, periodic boundary
 * conditions are used but it is user's responsibility to provide input
 * whose \f$\phi\f$-grid makes sense (in that it actually represents a periodic
 * field), i.e., \f$\phi_\mathrm{max}-\phi_\mathrm{min} = 2\pi/(N+1)\f$.
 * However, do note that in this module \f$\phi_\mathrm{max}\f$ is not the
 * "last" grid point but the second last, e.g. if \f$\phi_\mathrm{min}=0\f$
 * and \f$n_\phi = 360\f$, then \f$\phi_\mathrm{max}=359\f$ if periodicity is
 * \f$N=0\f$.
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "E_3DS.h"
#include "../spline/interp.h"

/**
 * @brief Initialize electric field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct.
 *
 * The offload data struct must have the following fields initialized:
 * - E_3DS_offload_data.n_r
 * - E_3DS_offload_data.n_z
 * - E_3DS_offload_data.r_min
 * - E_3DS_offload_data.r_max
 * - E_3DS_offload_data.z_min
 * - E_3DS_offload_data.z_max
 * - E_3DS_offload_data.n_phi
 * - E_3DS_offload_data.phi_min
 * - E_3DS_offload_data.phi_max
 *
 * E_3DS_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[                     z*En_r*En_z + j*En_r + i]
 *   = E_R(R_i, phi_z, z_j)   [V/m]
 * - offload_array[  En_r*En_z*En_phi + z*En_r*En_z + j*En_r + i]
 *   = E_phi(R_i, phi_z, z_j)   [V/m]
 * - offload_array[2*En_r*En_z*En_phi + z*En_r*En_z + j*En_r + i]
 *   = E_z(R_i, phi_z, z_j)   [V/m]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int E_3DS_init_offload(E_3DS_offload_data* offload_data, real** offload_array) {

    /* Spline initialization. */
    int err = 0;
    int E_size = offload_data->n_r * offload_data->n_z
      * offload_data->n_phi;

    /* Allocate enough space to store three 3D arrays and one 2D array */
    real* coeff_array = (real*) malloc( (3*NSIZE_COMP3D*E_size)*sizeof(real));
    real* E_r   = &(coeff_array[0*E_size*NSIZE_COMP3D]);
    real* E_phi = &(coeff_array[1*E_size*NSIZE_COMP3D]);
    real* E_z   = &(coeff_array[2*E_size*NSIZE_COMP3D]);

    err += interp3Dcomp_init_coeff(
        E_r, *offload_array + 0*E_size,
        offload_data->n_r, offload_data->n_phi,
        offload_data->n_z,
        NATURALBC, PERIODICBC, NATURALBC,
        offload_data->r_min,   offload_data->r_max,
        offload_data->phi_min, offload_data->phi_max,
        offload_data->z_min,   offload_data->z_max);

    err += interp3Dcomp_init_coeff(
        E_phi, *offload_array + 1*E_size,
        offload_data->n_r, offload_data->n_phi,
        offload_data->n_z,
        NATURALBC, PERIODICBC, NATURALBC,
        offload_data->r_min,   offload_data->r_max,
        offload_data->phi_min, offload_data->phi_max,
        offload_data->z_min,   offload_data->z_max);

    err += interp3Dcomp_init_coeff(
        E_z, *offload_array + 2*E_size,
        offload_data->n_r, offload_data->n_phi,
        offload_data->n_z,
        NATURALBC, PERIODICBC, NATURALBC,
        offload_data->r_min,   offload_data->r_max,
        offload_data->phi_min, offload_data->phi_max,
        offload_data->z_min,   offload_data->z_max);

    if(err) {
      print_err("Error: Failed to initialize splines.\n");
      return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);

    *offload_array = coeff_array;
        offload_data->offload_array_length = NSIZE_COMP3D*E_size*3;

    /* Print some sanity check on data */
    print_out(VERBOSE_IO,
              "\n3D electric field, tricubic interpolation (E_3DS)\n");
    print_out(VERBOSE_IO, "Grid: nR = %4.d Rmin = %3.3f Rmax = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f zmax = %3.3f\n",
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max);
    print_out(VERBOSE_IO, "     nphi = %4.d phimin = %3.3f phimax = %3.3f\n",
              offload_data->n_phi,
              offload_data->phi_min, offload_data->phi_max);

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
void E_3DS_free_offload(
    E_3DS_offload_data* offload_data, real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize electric field data struct on target
 *
 * This function copies the electric field parameters from the offload struct
 * to the struct on target and sets the electric field data pointers to
 * correct offsets in the offload array.
 *
 * @param EData pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void E_3DS_init(E_3DS_data* Edata, E_3DS_offload_data* offload_data,
               real* offload_array) {

    int E_size = NSIZE_COMP3D * offload_data->n_r
        * offload_data->n_z * offload_data->n_phi;

    /* Initialize spline structs from the coefficients */
    interp3Dcomp_init_spline(&Edata->E_r, &(offload_array[0*E_size]),
                             offload_data->n_r,
                             offload_data->n_phi,
                             offload_data->n_z,
                             NATURALBC, PERIODICBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->phi_min,
                             offload_data->phi_max,
                             offload_data->z_min,
                             offload_data->z_max);

    interp3Dcomp_init_spline(&Edata->E_phi, &(offload_array[1*E_size]),
                             offload_data->n_r,
                             offload_data->n_phi,
                             offload_data->n_z,
                             NATURALBC, PERIODICBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->phi_min,
                             offload_data->phi_max,
                             offload_data->z_min,
                             offload_data->z_max);

    interp3Dcomp_init_spline(&Edata->E_z, &(offload_array[2*E_size]),
                             offload_data->n_r,
                             offload_data->n_phi,
                             offload_data->n_z,
                             NATURALBC, PERIODICBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->phi_min,
                             offload_data->phi_max,
                             offload_data->z_min,
                             offload_data->z_max);

}

/**
 * @brief Evaluate electric field
 *
 * This function evaluates the electric field at the given coordinates using
 * trilinear interpolation on the 3D electric field data.
 *
 * @param E value will be stored in E[1] E[2] E[3]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Edata pointer to electric field data struct
 *
 */
a5err E_3DS_eval_E(real E[3], real r, real phi, real z,
                   E_3DS_data* Edata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += interp3Dcomp_eval_f(&E[0], &Edata->E_r, r, phi, z);
    interperr += interp3Dcomp_eval_f(&E[1], &Edata->E_phi, r, phi, z);
    interperr += interp3Dcomp_eval_f(&E[2], &Edata->E_z, r, phi, z);

    /* Test for E field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_3DS );
    }

    return err;
}
