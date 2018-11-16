/**
 * @file N0_ST.c @brief Stellarator neutral density with trilinear interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../consts.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "N0_ST.h"
#include "../linint/linint3D.h" /* for 3D interpolation routines */

int N0_ST_init_offload(N0_ST_offload_data* offload_data,
                        real** offload_array) {

    /* We need to use stellarator symmetry here.
     * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
     * The data is expected to include half a period.
     */
    real* temp_n0 = *offload_array;
    int n0_size = offload_data->n_r * offload_data->n_z * (2*(offload_data->n_phi - 1));
    *offload_array = (real*) malloc(n0_size * sizeof(real));
    offload_data->offload_array_length = n0_size;

    int i_phi;
    int i_z;
    int i_r;
    int temp_ind, off_ind, sym_ind;
    for (i_phi = 0; i_phi < offload_data->n_phi; i_phi++) {
        for (i_z = 0; i_z < offload_data->n_z; i_z++) {
            for (i_r = 0; i_r < offload_data->n_r; i_r++) {
                /* Stellarator symmetry: i_phi        <=> 2*(n_phi-1)-i_phi
                 *                       i_z          <=> n_z-i_z-1
                 * So, a point at:      (i_r, i_phi, i_z)  <=> (i_r, 2*(n_phi-1)-i_phi, n_z-i_z-1)
                 *
                 * temp_n0 data is in the format: (i_r, i_phi, i_z) = temp_n0(i_z*n_phi*n_r + i_phi*n_r + i_r)
                 * offload_array data -"-"-"-"-  : (i_r, i_phi, i_z) = (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ]
                 * => (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ] = temp_n0(i_z*n_phi*n_r + i_phi*n_r + i_r);
                 * The values are: Sym[B_r, B_phi, B_z] = [-B_r, B_phi, B_z]
                 */

                /* Index of data point in temp arrays */
                temp_ind = i_z*offload_data->n_phi*offload_data->n_r + i_phi*offload_data->n_r + i_r;

                /* Index of data point in offload_array and corresponding stel.-symmetric index  */
                off_ind = i_phi*offload_data->n_z*offload_data->n_r + i_z*offload_data->n_r + i_r;
                sym_ind = (2*(offload_data->n_phi - 1) - i_phi)*offload_data->n_z*offload_data->n_r
                    + (offload_data->n_z - i_z - 1)*offload_data->n_r + i_r;

                (*offload_array)[off_ind] =  temp_n0[temp_ind];
                if (i_phi != 0) {
                    (*offload_array)[sym_ind] = -temp_n0[temp_ind];
                }
            }
        }
    }
    /* Phi data is now for one toroidal period */
    offload_data->n_phi = 2*(offload_data->n_phi - 1);
    offload_data->phi_max = 2*offload_data->phi_max;

    print_out(VERBOSE_IO, "\n3D neutral density (N0_3D)\n");
    print_out(VERBOSE_IO, "Grid: nR = %4.d Rmin = %3.3f Rmax = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f zmax = %3.3f\n",
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max);
    print_out(VERBOSE_IO, "     nphi = %4.d phimin = %3.3f phimax = %3.3f\n",
              offload_data->n_phi,
              offload_data->phi_min, offload_data->phi_max);

    free(temp_n0);

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
void N0_ST_free_offload(N0_ST_offload_data* offload_data,
                        real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize magnetic field data struct on target
 *
 * This function copies the magnetic field parameters from the offload struct
 * to the struct on target and sets the magnetic field data pointers to
 * correct offsets in the offload array.
 *
 * @param BData pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
int N0_ST_init(N0_ST_data* ndata, N0_ST_offload_data* offload_data,
               real* offload_array) {
    int err = 0;

    ndata->periods = offload_data->periods;

    err += linint3D_init(
        &ndata->n0, offload_array,
        offload_data->n_r, offload_data->n_phi, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    return err;
}

/**
 * @brief Evaluate neutral density
 *
 * This function evaluates the neutral density n0 at the given coordinates using
 * trilinear interpolation on the stellarator neutral density data.
 *
 * @param n0 n0 value will be stored in n0[0]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param ndata pointer to neutral density data struct
 *
 */
a5err N0_ST_eval_n0(real n0[], real r, real phi, real z,
                   N0_ST_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    phi = fmod(phi, 2*CONST_PI/ndata->periods);
    /* if(phi < ndata->n0.phi_min) { */
    if(phi < 0) {
        phi += 2*CONST_PI/ndata->periods;
    }
    if(phi > CONST_PI/ndata->periods) {
        /* Stellarator-symmetric mirroring */
        phi = 2*CONST_PI/ndata->periods - phi;
        z = -z;
    }

    interperr += linint3D_eval(&n0[0], &ndata->n0, r, phi, z);

    if(interperr) {err = error_raise( ERR_OUTSIDE_N0DATA, __LINE__ );}

    return err;
}
