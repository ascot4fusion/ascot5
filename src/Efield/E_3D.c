/**
 * @file E_3D.c @brief 3D electric field with trilinear interpolation
 *
 * This module represents an electric field where data is given in \f$R\phi z\f$-
 * grid from which it is interpolated with trilinear splines.
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
 * @see E_field.c linint3D.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "E_3D.h"
#include "../linint/linint.h"

/**
 * @brief Initialize electric field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct.
 *
 * The offload data struct must have the following fields initialized:
 * - E_3D_offload_data.n_r
 * - E_3D_offload_data.n_z
 * - E_3D_offload_data.r_min
 * - E_3D_offload_data.r_max
 * - E_3D_offload_data.z_min
 * - E_3D_offload_data.z_max
 * - E_3D_offload_data.n_phi
 * - E_3D_offload_data.phi_min
 * - E_3D_offload_data.phi_max
 *
 * E_3D_offload_data.offload_array_length is set here.
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int E_3D_init_offload(E_3D_offload_data* offload_data, real** offload_array) {

    /* Fill rest of the offload data struct */
    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
        / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
        / (offload_data->n_z - 1);
    offload_data->phi_grid = (offload_data->phi_max - offload_data->phi_min)
        / (offload_data->n_phi - 1);

    offload_data->offload_array_length =
        3*offload_data->n_r * offload_data->n_phi * offload_data->n_z;

    print_out(VERBOSE_IO, "\n3D electric field, trilinear interpolation (E_3D)\n");
    print_out(VERBOSE_IO, "Grid: nR = %4.d Rmin = %3.3f Rmax = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f zmax = %3.3f\n",
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max);
    print_out(VERBOSE_IO, "     nphi = %4.d phimin = %3.3f phimax = %3.3f\n",
              offload_data->n_phi,
              offload_data->phi_min, offload_data->phi_max);

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
void E_3D_free_offload(E_3D_offload_data* offload_data, real** offload_array) {
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
 * Any initialization that requires any computations must have been done already
 * when the offload struct was initialized. The linear interpolation
 * can be done here because it does not require any computation.
 *
 * The offload array must contain the following data:
 * - offload_array[                     z*En_r*En_z + j*En_r + i]
 *   = E_R(R_i, phi_z, z_j)   [V/m]
 * - offload_array[  En_r*En_z*En_phi + z*En_r*En_z + j*En_r + i]
 *   = E_phi(R_i, phi_z, z_j)   [V/m]
 * - offload_array[2*En_r*En_z*En_phi + z*En_r*En_z + j*En_r + i]
 *   = E_z(R_i, phi_z, z_j)   [V/m]
 *
 * @param EData pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void E_3D_init(E_3D_data* Edata, E_3D_offload_data* offload_data,
               real* offload_array) {
    int E_size = offload_data->n_r * offload_data->n_z
        * offload_data->n_phi;

    linint3D_init(&Edata->E_r, &offload_array[E_size*0],
                  offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
                  offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                  offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    linint3D_init(&Edata->E_phi, &offload_array[E_size*1],
                  offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
                  offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                  offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    linint3D_init(&Edata->E_z, &offload_array[E_size*2],
                  offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
                  offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                  offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    return;
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
a5err E_3D_eval_E(real E[3], real r, real phi, real z,
                   E_3D_data* Edata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint3D_eval_f(&E[0], &Edata->E_r, r, phi, z);
    interperr += linint3D_eval_f(&E[1], &Edata->E_phi, r, phi, z);
    interperr += linint3D_eval_f(&E[2], &Edata->E_z, r, phi, z);


    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_3D);}

    return err;
}
