/**
 * @file B_ST.c
 * @brief Stellarator magnetic field with cubic interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "math.h"
#include "ascot5.h"
#include "B_ST.h"
#include "B_3D.h" /* for 3D interpolation routines */
#include "B_2D.h" /* for 2D interpolation routines */
#include "hdf5_bfield.h"

/**
 * @brief Load magnetic field data and prepare parameters
 *
 * This function reads the magnetic field data from the input.h5 file,
 * fills the offload struct with parameters and 
 * allocates and fills the offload array.
 *
 * @todo Error checking
 * @todo Move reading the file to ascot4_interface
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_ST_init_offload(B_ST_offload_data* offload_data, real** offload_array) {

    hdf5_bfield_init_offload_ST(offload_data, offload_array);

}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_ST_free_offload(B_ST_offload_data* offload_data, real** offload_array) {
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
void B_ST_init(B_ST_data* Bdata, B_ST_offload_data* offload_data,
               real* offload_array) {
    Bdata->n_r = offload_data->n_r;
    Bdata->n_z = offload_data->n_z;
    Bdata->n_phi = offload_data->n_phi;
    Bdata->periods = offload_data->periods;
    Bdata->r_min = offload_data->r_min;
    Bdata->r_max = offload_data->r_max;
    Bdata->r_grid = offload_data->r_grid;
    Bdata->z_min = offload_data->z_min;
    Bdata->z_max = offload_data->z_max;
    Bdata->z_grid = offload_data->z_grid;
    Bdata->phi_min = offload_data->phi_min;
    Bdata->phi_max = offload_data->phi_max;
    Bdata->phi_grid = offload_data->phi_grid;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;
    Bdata->B_r = &offload_array[0];
    Bdata->B_phi = &offload_array[(offload_data->n_phi+4)*offload_data->n_r*offload_data->n_z];
    Bdata->B_z = &offload_array[2*(offload_data->n_phi+4)*offload_data->n_r*offload_data->n_z];
    Bdata->s = &offload_array[3*(offload_data->n_phi+4)*offload_data->n_r*offload_data->n_z];
}

/**
 * @brief Evaluate magnetic field
 *
 * This function evaluates the magnetic field at the given coordinates using
 * tricubic interpolation on the 3D magnetic field data. This is a SIMD
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param B array where magnetic field values will be stored (Br -> B[0][i],
 *          Bphi -> B[1][i], Bz -> B[2][i]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 */
void B_ST_eval_B(real B[], real r, real phi, real z, B_ST_data* Bdata) {
    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }
    int i_r = (int) floor((r - Bdata->r_min)
                    / ((Bdata->r_max - Bdata->r_min)
                       / (Bdata->n_r-1)));
    if(i_r < 0 || i_r >= Bdata->n_r) {
        i_r = 0;
    }

    int i_phi = (int) floor((phi - Bdata->phi_min)
                    / ((Bdata->phi_max - Bdata->phi_min)
                       / (Bdata->n_phi+3)));
    if(i_phi < 0 || i_phi >= Bdata->n_phi+4) {
        i_phi = 0;
    }

    int i_z = (int) floor((z - Bdata->z_min)
                    / ((Bdata->z_max - Bdata->z_min)
                       / (Bdata->n_z-1)));
    if(i_z < 0 || i_z >= Bdata->n_z) {
        i_z = 0;
    }

    real t_r = (r - (Bdata->r_min + i_r * Bdata->r_grid)) / Bdata->r_grid;
    real t_phi = (phi - (Bdata->phi_min + i_phi * Bdata->phi_grid)) / Bdata->phi_grid;
    real t_z = (z - (Bdata->z_min + i_z * Bdata->z_grid)) / Bdata->z_grid;

    B[0] = B_3D_tricubic(t_r, t_phi, t_z, i_r, i_phi, i_z, Bdata->n_z, Bdata->n_r, Bdata->B_r);
    B[1] = B_3D_tricubic(t_r, t_phi, t_z, i_r, i_phi, i_z, Bdata->n_z, Bdata->n_r, Bdata->B_phi);
    B[2] = B_3D_tricubic(t_r, t_phi, t_z, i_r, i_phi, i_z, Bdata->n_z, Bdata->n_r, Bdata->B_z);
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * This function evaluates the poloidal flux psi at the given coordinates using
 * tricubic interpolation on the stellarator 3D s data. This is a SIMD
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param psi psi value will be stored in psi[0][i]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Change to a scalar elemental function and compare performance
 */
void B_ST_eval_psi(real psi[], real r, real phi, real z,
                   B_ST_data* Bdata)
{
    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }
    int i_r = (int) floor((r - Bdata->r_min)
                    / ((Bdata->r_max - Bdata->r_min)
                       / (Bdata->n_r-1)));
    if(i_r < 0 || i_r >= Bdata->n_r) {
        i_r = 0;
    }

    int i_phi = (int) floor((phi - Bdata->phi_min)
                    / ((Bdata->phi_max - Bdata->phi_min)
                       / (Bdata->n_phi+3)));
    if(i_phi < 0 || i_phi >= Bdata->n_phi+4) {
        i_phi = 0;
    }

    int i_z = (int) floor((z - Bdata->z_min)
                    / ((Bdata->z_max - Bdata->z_min)
                       / (Bdata->n_z-1)));
    if(i_z < 0 || i_z >= Bdata->n_z) {
        i_z = 0;
    }

    real t_r = (r - (Bdata->r_min + i_r * Bdata->r_grid)) / Bdata->r_grid;
    real t_phi = (phi - (Bdata->phi_min + i_phi * Bdata->phi_grid)) / Bdata->phi_grid;
    real t_z = (z - (Bdata->z_min + i_z * Bdata->z_grid)) / Bdata->z_grid;

    psi[0] = B_3D_tricubic(t_r, t_phi, t_z, i_r, i_phi, i_z, Bdata->n_z, Bdata->n_r, Bdata->s);
}

/**
 * @brief Evaluate radial coordinate rho
 *
 * This function evaluates the radial coordinate rho at the given psi value.
 * This is a SIMD function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param rho rho value will be stored in rho[0][i]
 * @param psi poloidal flux value 
 * @param Bdata pointer to magnetic field data struct
 *
 */
void B_ST_eval_rho(real rho[], real psi, B_ST_data* Bdata) {

    rho[0] = sqrt(psi);
    
}


/**
 * @brief Evaluate magnetic field and derivatives
 *
 * This function evaluates the magnetic field and it's derivatives at the given 
 * coordinates using bicubic interpolation on the stellarator magnetic field data.
 * This is a SIMD function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param B_dB array where magnetic field values will be stored (Br -> B[0][i],
 *          dBr/dr -> B[1][i], dBr/dphi -> B[2][i], dBr/dz -> B[3][i],
 *          Bphi -> B[4][i], dBphi/dr -> B[5][i], dBphi/dphi -> B[6][i],
 *          dBphi/dz -> B[7][i], Bz -> B[8][i], dBz/dr -> B[9][i],
 *          dBz/dphi -> B[10][i], dBz/dz -> B[11][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 */
void B_ST_eval_B_dB(real B_dB[], real r, real phi, real z, B_ST_data* Bdata) {
    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    int i_r = (int) floor((r - Bdata->r_min)
                    / ((Bdata->r_max - Bdata->r_min)
                       / (Bdata->n_r-1)));
    if(i_r < 0 || i_r >= Bdata->n_r) {
        i_r = 0;
    }

    int i_phi = (int) floor((phi - Bdata->phi_min)
                    / ((Bdata->phi_max - Bdata->phi_min)
                       / (Bdata->n_phi+3)));
    if(i_phi < 0 || i_phi >= Bdata->n_phi+4) {
        i_phi = 0;
    }

    int i_z = (int) floor((z - Bdata->z_min)
                    / ((Bdata->z_max - Bdata->z_min)
                       / (Bdata->n_z-1)));
    if(i_z < 0 || i_z >= Bdata->n_z) {
        i_z = 0;
    }

    real t_r = (r - (Bdata->r_min + i_r * Bdata->r_grid)) / Bdata->r_grid;
    real t_phi = (phi - (Bdata->phi_min + i_phi * Bdata->phi_grid)) / Bdata->phi_grid;
    real t_z = (z - (Bdata->z_min + i_z * Bdata->z_grid)) / Bdata->z_grid;

    B_3D_tricubic_derivs(&B_dB[0], t_r, t_phi, t_z, i_r, i_phi, i_z,
                         Bdata->n_r, Bdata->n_z, Bdata->r_grid, Bdata->phi_grid,
                         Bdata->z_grid, Bdata->B_r);
    B_3D_tricubic_derivs(&B_dB[4], t_r, t_phi, t_z, i_r, i_phi, i_z,
                         Bdata->n_r, Bdata->n_z, Bdata->r_grid, Bdata->phi_grid,
                         Bdata->z_grid, Bdata->B_phi);
    B_3D_tricubic_derivs(&B_dB[8], t_r, t_phi, t_z, i_r, i_phi, i_z,
                         Bdata->n_r, Bdata->n_z, Bdata->r_grid, Bdata->phi_grid,
                         Bdata->z_grid, Bdata->B_z);
}

real B_ST_get_axis_r(B_ST_data* Bdata) {
    // 3D magnetic axis not implemented yet
    return 0;
}

real B_ST_get_axis_z(B_ST_data* Bdata) {
    // 3D magnetic axis not implemented yet
    return 0;
}
