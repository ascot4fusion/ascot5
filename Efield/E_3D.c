x/**
 * @file E_3D.c @brief 3D electric field with trilinear interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "E_3D.h"
#include "../linint/linint3D.h" /* for 3D interpolation routines */

void E_3D_init_offload(E_3D_offload_data* offload_data, real** offload_array) {
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
 * @param BData pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
int N0_3D_init(E_3D_data* Edata, E_3D_offload_data* offload_data,
               real* offload_array) {
    int err = 0;
    
    err += linint3D_init(&ndata->E_r, offload_array,
                         offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                         offload_data->r_min, offload_data->r_max, offload_data->r_grid,
                         offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                         offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    err += linint3D_init(&ndata->E_phi, offload_array,
                         offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                         offload_data->r_min, offload_data->r_max, offload_data->r_grid,
                         offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                         offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += linint3D_init(&ndata->E_z, offload_array,
                         offload_data->n_r, offload_data->n_phi, offload_data->n_z,
                         offload_data->r_min, offload_data->r_max, offload_data->r_grid,
                         offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
                         offload_data->z_min, offload_data->z_max, offload_data->z_grid);


    return err;
}

/**
 * @brief Evaluate electric field
 *
 * This function evaluates the electric field at the given coordinates using
 * trilinear interpolation on the 3D electric field data.
 *
 * @param n0 n0 value will be stored in n0[0]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param ndata pointer to neutral density data struct
 *
 */
a5err E_3D_eval_E(real E[], real r, real phi, real z,
                   E_3D_data* Edata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += linint3D_eval(&E[0], &Edata->E_r, r, phi, z);
    interperr += linint3D_eval(&E[1], &Edata->E_phi, r, phi, z);
    interperr += linint3D_eval(&E[2], &Edata->E_z, r, phi, z);


    if(interperr) {err = error_raise( ERR_OUTSIDE_N0DATA, __LINE__ );}

    return err;
}

a5err E_3D_eval_E_SIMD(int i, real E[3][NSIMD], real r, real phi, real z,
			  E_3D_data* Edata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += linint3D_eval_SIMD(i, E[0], &Edata->E_r, r, phi, z);
    interperr += linint3D_eval_SIMD(i, E[1], &Edata->E_phi, r, phi, z);
    interperr += linint3D_eval_SIMD(i, E[2], &Edata->E_z, r, phi, z);

    if(interperr) {err = error_raise( ERR_OUTSIDE_N0DATA, __LINE__ );}

    return err;
}

