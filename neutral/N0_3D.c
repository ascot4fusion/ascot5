/**
 * @file N0_3D.c @brief 3D neutral density with trilinear interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "N0_3D.h"
#include "../linint/linint3D.h" /* for 3D interpolation routines */

void N0_3D_init_offload(N0_3D_offload_data* offload_data, real** offload_array) {
    /* In hdf5 directory */
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void N0_3D_free_offload(N0_3D_offload_data* offload_data, real** offload_array) {
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
int N0_3D_init(N0_3D_data* ndata, N0_3D_offload_data* offload_data,
               real* offload_array) {
    int err = 0;
    
    err += linint3D_init(&ndata->n0, offload_array,
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
 * trilinear interpolation on the 3D neutral density data.
 *
 * @param n0 n0 value will be stored in n0[0]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param ndata pointer to neutral density data struct
 *
 */
a5err N0_3D_eval_n0(real n0[], real r, real phi, real z,
                   N0_3D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint3D_eval(&n0[0], &ndata->n0, r, phi, z);

    if(interperr) {err = error_raise( ERR_OUTSIDE_N0DATA, __LINE__ );}

    return err;
}

a5err N0_3D_eval_n0_SIMD(int i, real n0[NSIMD], real r, real phi, real z,
			  N0_3D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += linint3D_eval_SIMD(i, n0, &ndata->n0, r, phi, z);

    if(interperr) {err = error_raise( ERR_OUTSIDE_N0DATA, __LINE__ );}

    return err;
}

