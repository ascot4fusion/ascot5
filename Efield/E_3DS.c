/**
 * @file E_3DS.c @brief 3D electric field with trilinear interpolation
 *
 * This module represents an electric field where data is given in \f$R\phi z\f$-
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
 * \f$N=0\f$.
 *
 * The splines may either have compact or explicit forms which is toggled by
 * INTERP_SPL_EXPL in ascot5.h. Compact forms require 1/8 th of memory (in 3D)
 * but require more floating point operations.
 *
 * @see E_field.c linint3D.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "E_3DS.h"
#include "../spline/interp3D.h"
#include "../spline/interp3Dcomp.h"
#include "../spline/interp3Dexpl.h"

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
    
    /* Fill rest of the offload data struct */
    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
	/ (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
	/ (offload_data->n_z - 1);
    offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
	/ (offload_data->n_phi - 1);

    /* Spline initialization. Use spline structs for temporary storage */
    int err = 0;
    int E_size = offload_data->n_r * offload_data->n_z
      * offload_data->n_phi;

    interp3D_data E_r;
    interp3D_data E_phi;
    interp3D_data E_z;

#if INTERP_SPL_EXPL
    err += interp3Dexpl_init(
			     &E_r, *offload_array + 0*E_size,
			     offload_data->n_r, offload_data->n_phi, offload_data->n_z,
			     offload_data->r_min, offload_data->r_max,
			     offload_data->r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->z_min, offload_data->z_max,
			     offload_data->z_grid);

    err += interp3Dexpl_init(
			     &E_phi, *offload_array + 1*E_size,
                             offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                             offload_data->r_min, offload_data->r_max,
                             offload_data->r_grid,
                             offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                             offload_data->z_min, offload_data->z_max,
                             offload_data->z_grid);

    err += interp3Dexpl_init(
			     &E_z, *offload_array + 2*E_size,
                             offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                             offload_data->r_min, offload_data->r_max,
                             offload_data->r_grid,
                             offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                             offload_data->z_min, offload_data->z_max,
                             offload_data->z_grid);

    /* The data is now presented with splines, each data point has
     * 64 spline coefficients in 3D */
    E_size *= 64;

#else
    err += interp3Dcomp_init(
                             &E_r, *offload_array + 0*E_size,
                             offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                             offload_data->r_min, offload_data->r_max,
                             offload_data->r_grid,
                             offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                             offload_data->z_min, offload_data->z_max,
                             offload_data->z_grid);

    err += interp3Dcomp_init(
                             &E_phi, *offload_array + 1*E_size,
                             offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                             offload_data->r_min, offload_data->r_max,
                             offload_data->r_grid,
                             offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                             offload_data->z_min, offload_data->z_max,
                             offload_data->z_grid);

    err += interp3Dcomp_init(
                             &E_z, *offload_array + 2*E_size,
                             offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                             offload_data->r_min, offload_data->r_max,
                             offload_data->r_grid,
                             offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                             offload_data->z_min, offload_data->z_max,
                             offload_data->z_grid);

    /* The data is now presented with splines, each data point has
     * 8 spline coefficients in 3D */
    E_size *= 8;

#endif
    if(err) {
      print_err("Error: Failed to initialize splines.\n");
      return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);

    offload_data->offload_array_length = 3*E_size;
    *offload_array =
      (real*) malloc( offload_data->offload_array_length * sizeof(real) );

    for(int i = 0; i < E_size; i++) {
      (*offload_array)[0*E_size + i] = E_r.c[i];
      (*offload_array)[1*E_size + i] = E_phi.c[i];
      (*offload_array)[2*E_size + i] = E_z.c[i];
    }

    interp3D_free(&E_r);
    interp3D_free(&E_phi);
    interp3D_free(&E_z);

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\3D electric field, tricubic interpolation (E_3DS)\n");
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
void E_3DS_free_offload(E_3DS_offload_data* offload_data, real** offload_array) {
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

    int E_size = offload_data->n_r * offload_data->n_z
      * offload_data->n_phi ;
#if INTERP_SPL_EXPL
    E_size   *= 64;
#else
    E_size   *= 8;
#endif

    /* Copy parameters and assign pointers to offload array to initialize the
       spline structs */
    Edata->E_r.n_r        = offload_data->n_r;
    Edata->E_r.r_min      = offload_data->r_min;
    Edata->E_r.r_max      = offload_data->r_max;
    Edata->E_r.r_grid     = offload_data->r_grid;
    Edata->E_r.n_z        = offload_data->n_z;
    Edata->E_r.z_min      = offload_data->z_min;
    Edata->E_r.z_max      = offload_data->z_max;
    Edata->E_r.z_grid     = offload_data->z_grid;
    Edata->E_r.n_phi      = offload_data->n_phi;
    Edata->E_r.phi_min    = offload_data->phi_min;
    Edata->E_r.phi_max    = offload_data->phi_max;
    Edata->E_r.phi_grid   = offload_data->phi_grid;
    Edata->E_r.c          = &(offload_array[0*E_size]);

    Edata->E_phi.n_r      = offload_data->n_r;
    Edata->E_phi.r_min    = offload_data->r_min;
    Edata->E_phi.r_max    = offload_data->r_max;
    Edata->E_phi.r_grid   = offload_data->r_grid;
    Edata->E_phi.n_z      = offload_data->n_z;
    Edata->E_phi.z_min    = offload_data->z_min;
    Edata->E_phi.z_max    = offload_data->z_max;
    Edata->E_phi.z_grid   = offload_data->z_grid;
    Edata->E_phi.n_phi    = offload_data->n_phi;
    Edata->E_phi.phi_min  = offload_data->phi_min;
    Edata->E_phi.phi_max  = offload_data->phi_max;
    Edata->E_phi.phi_grid = offload_data->phi_grid;
    Edata->E_phi.c        = &(offload_array[1*E_size]);

    Edata->E_z.n_r        = offload_data->n_r;
    Edata->E_z.r_min      = offload_data->r_min;
    Edata->E_z.r_max      = offload_data->r_max;
    Edata->E_z.r_grid     = offload_data->r_grid;
    Edata->E_z.n_z        = offload_data->n_z;
    Edata->E_z.z_min      = offload_data->z_min;
    Edata->E_z.z_max      = offload_data->z_max;
    Edata->E_z.z_grid     = offload_data->z_grid;
    Edata->E_z.n_phi      = offload_data->n_phi;
    Edata->E_z.phi_min    = offload_data->phi_min;
    Edata->E_z.phi_max    = offload_data->phi_max;
    Edata->E_z.phi_grid   = offload_data->phi_grid;
    Edata->E_z.c          = &(offload_array[2*E_size]);

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
#if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_B(&E[0], &Edata->E_r, r, phi, z);
    interperr += interp3Dexpl_eval_B(&E[1], &Edata->E_phi, r, phi, z);
    interperr += interp3Dexpl_eval_B(&E[2], &Edata->E_z, r, phi, z);
#else
    interperr += interp3Dcomp_eval_B(&E[0], &Edata->E_r, r, phi, z);
    interperr += interp3Dcomp_eval_B(&E[1], &Edata->E_phi, r, phi, z);
    interperr += interp3Dcomp_eval_B(&E[2], &Edata->E_z, r, phi, z);
#endif

    /* Test for E field interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_3DS );}

    return err;
}


