/**
 * @file E_3DST.c @brief Time dependent 3D electric field with 4D spline interpolation
 *
 * This module represents an electric field where data is given in \f$R\phi z t\f$-
 * grid from which it is interpolated with tricubic splines.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown. For \f$\phi\f$-grid, periodic boundary
 * conditions are used but it is user's responsibility to provide input
 * whose \f$\phi\f$-grid makes sense (in that it actually represents a periodic
 * field), i.e., \f$\phi_\mathrm{max}-\phi_\mathrm{min} = 2\pi/(N+1)\f$.
 * However, do note that in this module \f$\phi_\mathrm{max}\f$ is not the
 * "last" grid point but the second last, e.g. if \f$\phi_\mathrm{min}=0\f$
 * and \f$n_\phi = 360\f$, then \f$\phi_\mathrm{max}=359\f$ if periodicity is
 * \f$N=0\f$. This module does not extrapolate in the time domain so an error
 * will be thrown or the marker will be stopped.
 *
 * @see E_field.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "E_3DST.h"
#include "../spline/interp.h"

/**
 * @brief Initialize electric field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct.
 *
 * The offload data struct must have the following fields initialized:
 * - E_3DST_offload_data.n_r
 * - E_3DST_offload_data.n_z
 * - E_3DST_offload_data.r_min
 * - E_3DST_offload_data.r_max
 * - E_3DST_offload_data.z_min
 * - E_3DST_offload_data.z_max
 * - E_3DST_offload_data.n_phi
 * - E_3DST_offload_data.phi_min
 * - E_3DST_offload_data.phi_max
 * - E_3DST_offload_data.n_t
 * - E_3DST_offload_data.t_min
 * - E_3DST_offload_data.t_max
 *
 * E_3DST_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[                          m*En_r*En_z*En_phi + k*En_r*En_z + j*En_r + i]
 *   = E_R(R_i, phi_k, z_j, t_m)   [V/m]
 * - offload_array[1*En_r*En_z*En_phi*En_t + m*En_r*En_z*En_phi + k*En_r*En_z + j*En_r + i]
 *   = E_phi(R_i, phi_k, z_j, t_m)   [V/m]
 * - offload_array[2*En_r*En_z*En_phi*En_t + m*En_r*En_z*En_phi + k*En_r*En_z + j*En_r + i]
 *   = E_z(R_i, phi_k, z_j, t_m)   [V/m]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int E_3DST_init_offload(E_3DST_offload_data* offload_data, real** offload_array) {

    /* Spline initialization. */
    int err = 0;
    int E_size = offload_data->n_r * offload_data->n_z
      * offload_data->n_phi * offload_data->n_t;

    /* Allocate enough space to store three 3D arrays and one 2D array */
    real* coeff_array = (real*) malloc( (3*NSIZE_COMP4D*E_size)*sizeof(real));
    real* E_r   = &(coeff_array[0*E_size*NSIZE_COMP4D]);
    real* E_phi = &(coeff_array[1*E_size*NSIZE_COMP4D]);
    real* E_z   = &(coeff_array[2*E_size*NSIZE_COMP4D]);

    err += interp4Dcomp_init_coeff(
        E_r, *offload_array + 0*E_size,
        offload_data->n_r, offload_data->n_phi,
        offload_data->n_z, offload_data->n_t,
        NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
        offload_data->r_min,   offload_data->r_max,
        offload_data->phi_min, offload_data->phi_max,
        offload_data->z_min,   offload_data->z_max,
        offload_data->t_min,   offload_data->t_max);

    err += interp4Dcomp_init_coeff(
        E_phi, *offload_array + 1*E_size,
        offload_data->n_r, offload_data->n_phi,
        offload_data->n_z, offload_data->n_t,
        NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
        offload_data->r_min,   offload_data->r_max,
        offload_data->phi_min, offload_data->phi_max,
        offload_data->z_min,   offload_data->z_max,
        offload_data->t_min,   offload_data->t_max);

    err += interp4Dcomp_init_coeff(
        E_z, *offload_array + 2*E_size,
        offload_data->n_r, offload_data->n_phi,
        offload_data->n_z, offload_data->n_t,
        NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
        offload_data->r_min,   offload_data->r_max,
        offload_data->phi_min, offload_data->phi_max,
        offload_data->z_min,   offload_data->z_max,
        offload_data->t_min,   offload_data->t_max);

    if(err) {
      print_err("Error: Failed to initialize splines.\n");
      return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);

    *offload_array = coeff_array;
        offload_data->offload_array_length = NSIZE_COMP4D*E_size*3;

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nTime-dependent 3D electric field (B_3DST)\n");
    print_out(VERBOSE_IO, "Grid: nR = %4.d Rmin = %3.3f Rmax = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f zmax = %3.3f\n",
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max);
    print_out(VERBOSE_IO, "     nphi = %4.d phimin = %3.3f phimax = %3.3f\n",
              offload_data->n_phi,
              math_rad2deg(offload_data->phi_min),
              math_rad2deg(offload_data->phi_max));
    print_out(VERBOSE_IO, "     ntime = %4.d tmin = %3.3f tmax = %3.3f\n",
              offload_data->n_t,
              offload_data->t_min, offload_data->t_max);

    return err;


}

/**
 * @brief Free offload array and reset parameters.
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void E_3DST_free_offload(E_3DST_offload_data* offload_data, real** offload_array) {
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
void E_3DST_init(E_3DST_data* Edata, E_3DST_offload_data* offload_data,
               real* offload_array) {

    int E_size = NSIZE_COMP4D * offload_data->n_r
        * offload_data->n_z * offload_data->n_phi * offload_data->n_t;

    /* Initialize spline structs from the coefficients */
    interp4Dcomp_init_spline(&Edata->E_r, &(offload_array[0*E_size]),
                             offload_data->n_r,
                             offload_data->n_phi,
                             offload_data->n_z,
                             offload_data->n_t,
                             NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->phi_min,
                             offload_data->phi_max,
                             offload_data->z_min,
                             offload_data->z_max,
                             offload_data->t_min,
                             offload_data->t_max);

    interp4Dcomp_init_spline(&Edata->E_phi, &(offload_array[1*E_size]),
                             offload_data->n_r,
                             offload_data->n_phi,
                             offload_data->n_z,
                             offload_data->n_t,
                             NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->phi_min,
                             offload_data->phi_max,
                             offload_data->z_min,
                             offload_data->z_max,
                             offload_data->t_min,
                             offload_data->t_max);

    interp4Dcomp_init_spline(&Edata->E_z, &(offload_array[2*E_size]),
                             offload_data->n_r,
                             offload_data->n_phi,
                             offload_data->n_z,
                             offload_data->n_t,
                             NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->phi_min,
                             offload_data->phi_max,
                             offload_data->z_min,
                             offload_data->z_max,
                             offload_data->t_min,
                             offload_data->t_max);

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
a5err E_3DST_eval_E(real E[3], real r, real phi, real z, real t,
                   E_3DST_data* Edata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += interp4Dcomp_eval_f(&E[0], &Edata->E_r, r, phi, z, t);
    interperr += interp4Dcomp_eval_f(&E[1], &Edata->E_phi, r, phi, z, t);
    interperr += interp4Dcomp_eval_f(&E[2], &Edata->E_z, r, phi, z, t);

    /* Test for E field interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_3DST );}

    return err;
}
