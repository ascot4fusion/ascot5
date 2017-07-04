/**
 * @file B_2DS.c
 * @brief 2D magnetic field with cubic spline interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "B_2DS.h"
#include "splinePatrik/interp2D.h" /* for 2D interpolation routines */

/**
 * @brief Load magnetic field data and prepare parameters
 *
 * This function reads the magnetic field data from input.magn_bkg and 
 * input.magn_header files, fills the offload struct with parameters and 
 * allocates and fills the offload array.
 *
 * @todo Error checking
 * @todo Move reading the file to ascot4_interface
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_2DS_init_offload(B_2DS_offload_data* offload_data, real** offload_array) {
    int i;
    FILE* f = fopen("input.magn_bkg", "r");

    /* Skip first line */
    fscanf(f, "%*d %*d %*d %*d %*d");

    /* Read r and z parameters */
    fscanf(f, "%lf %lf %d", &(offload_data->r_min), &(offload_data->r_max),
                            &(offload_data->n_r));
    fscanf(f, "%lf %lf %d", &(offload_data->z_min), &(offload_data->z_max),
                            &(offload_data->n_z));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);

    /* Allocate offload_array; psi and each component (r,phi,z) is
     * size n_r*n_z */
    int B_size = offload_data->n_r * offload_data->n_z;
    *offload_array = (real*) malloc(4 * B_size * sizeof(real));
    offload_data->offload_array_length = 4 * B_size;
    
    /* Skip phimaps */
    fscanf(f, "%*f %*f");

    /* Read psi */
    for(i = 0; i < B_size; i++) {
        fscanf(f, "%lf", &(*offload_array)[i]);
    }

    /* Read B_r */
    for(i = 0; i < B_size; i++) {
        fscanf(f, "%lf", &(*offload_array)[i+B_size]);
    }

    /* Read B_phi */
    for(i = 0; i < B_size; i++)
        fscanf(f, "%lf", &(*offload_array)[i+2*B_size]);

    /* Read B_z */
    for(i = 0; i < B_size; i++) {
        fscanf(f, "%lf", &(*offload_array)[i+3*B_size]);
    }

    fclose(f);

    /* Read rho parameters from input.magn_header */
    f = fopen("input.magn_header", "r");

    /* Skip first four lines */
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');

    /* Read the first two values; These are the poloidal flux (psi) values at
     * magnetic axis and at x point (that is, separatrix). */
    fscanf(f, "%lf %lf", &(offload_data->psi0), &(offload_data->psi1));

    /* Read magnetic axis r and z coordinates */
    while(fgetc(f) != '\n');
    fscanf(f, "%lf", &(offload_data->axis_r));
    while(fgetc(f) != '\n');
    fscanf(f, "%lf", &(offload_data->axis_z));

    fclose(f);
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_2DS_free_offload(B_2DS_offload_data* offload_data, real** offload_array) {
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
void B_2DS_init(B_2DS_data* Bdata, B_2DS_offload_data* offload_data,
               real* offload_array) {
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;
    /* Spline initialization and storage. */
    interp2D_init(&Bdata->psi, offload_array,
    		  offload_data->n_r, offload_data->n_z,
    		  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
    		  offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    interp2D_init(&Bdata->B_r, offload_array+offload_data->n_z*offload_data->n_r,
    		  offload_data->n_r, offload_data->n_z,
    		  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
    		  offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    interp2D_init(&Bdata->B_phi, offload_array+2*offload_data->n_z*offload_data->n_r,
    		  offload_data->n_r, offload_data->n_z,
    		  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
    		  offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    interp2D_init(&Bdata->B_z, offload_array+3*offload_data->n_z*offload_data->n_r,
    		  offload_data->n_r, offload_data->n_z,
    		  offload_data->r_min, offload_data->r_max, offload_data->r_grid,
		  offload_data->z_min, offload_data->z_max, offload_data->z_grid);
}

/**
 * @brief Evaluate magnetic field
 *
 * This function evaluates the magnetic field at the given coordinates using
 * bicubic interpolation on the 2D magnetic field data. This is a SIMD
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
void B_2DS_eval_B(real B[], real r, real phi, real z, B_2DS_data* Bdata) {
    interp2D_eval_B(&B[0], &Bdata->B_r, r, z);
    interp2D_eval_B(&B[1], &Bdata->B_phi, r, z);
    interp2D_eval_B(&B[2], &Bdata->B_z, r, z);
    #ifndef NOPSI
    real psi_dpsi[4];
    B_2DS_eval_psi_dpsi(psi_dpsi, r, phi, z, Bdata);
    B[0] = B[0]+1/(2*math_pi)*-psi_dpsi[2]/r;
    B[2] = B[2]+1/(2*math_pi)*psi_dpsi[1]/r;
    #endif
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * This function evaluates the poloidal flux psi at the given coordinates using
 * bicubic interpolation on the 2D magnetic field data. This is a SIMD
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
void B_2DS_eval_psi(real psi[], real r, real phi, real z, B_2DS_data* Bdata)
{
    interp2D_eval_B(&psi[0], &Bdata->psi, r, z);
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * This function evaluates the poloidal flux psi and its derivatives at the given coordinates using
 * bicubic interpolation on the 2D magnetic field data. This is a SIMD
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param psi psi values (psi   -> psi_dpsi[0][i]    dpsi/dr -> psi_dpsi[1][i]
 *        dpsi/dphi -> psi_dpsi[2][i]    dpsi/dz -> psi_dpsi[3][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 */
void B_2DS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_2DS_data* Bdata) {
    real psi_dpsi_temp[6];
    interp2D_eval_dB(psi_dpsi_temp, &Bdata->psi, r, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];
}

/**
 * @brief Evaluate radial coordinate rho
 *
 * This function evaluates the radial coordinate rho at the given psi value
 * using linear interpolation. This is a SIMD function, so the values are 
 * placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param rho rho value will be stored in rho[0][i]
 * @param psi poloidal flux value 
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Change to a scalar elemental function and compare performance
 */
void B_2DS_eval_rho(real rho[], real psi, B_2DS_data* Bdata) {
    if(psi - Bdata->psi0 < 0) {
        rho[0] = 0;
    }
    else {
        rho[0] = sqrt((psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0));
    }
}

/**
 * @brief Evaluate radial coordinate rho and its derivatives
 *
 * This function evaluates the radial coordinate rho and its derivatives
 * at the given coordinates using bicubic interpolation on the 
 * 2D magnetic field data. This is a SIMD 
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param rho rho values (rho   -> rho_drho[0][i]    drho/dr -> rho_drho[1][i]
 *        drho/dphi -> rho_drho[2][i]    drho/dz -> rho_drho[3][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 */
void B_2DS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_2DS_data* Bdata) {
    real rho;
    B_2DS_eval_psi_dpsi(rho_drho, r, phi, z, Bdata);
    /* Convert: rho = sqrt(psi), drho = dpsi/(2 * sqrt(psi))
     * Note that rho_drho[2] = 1/R * drho/dphi, because of cylindrical gradient
     */
    rho = sqrt(rho_drho[0]);
    rho_drho[0] = rho;
    rho_drho[1] = rho_drho[1] / (2*rho);
    rho_drho[2] = rho_drho[2] / (2*rho);
    rho_drho[3] = rho_drho[3] / (2*rho);
}

/**
 * @brief Evaluate magnetic field and derivatives
 *
 * This function evaluates the magnetic field and it's derivatives at the given 
 * coordinates using bicubic interpolation on the 2D magnetic field data. This 
 * is a SIMD function, so the values are placed in an NSIMD length struct.
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
void B_2DS_eval_B_dB(real B_dB[], real r, real phi, real z, B_2DS_data* Bdata) {
    real B_dB_temp[6];
    interp2D_eval_dB(B_dB_temp, &Bdata->B_r, r, z);
    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = 0;
    B_dB[3] = B_dB_temp[2];
    interp2D_eval_dB(B_dB_temp, &Bdata->B_phi, r, z);
    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = 0;
    B_dB[7] = B_dB_temp[2];
    interp2D_eval_dB(B_dB_temp, &Bdata->B_z, r, z);
    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = 0;
    B_dB[11] = B_dB_temp[2];
    #ifndef NOPSI
    real psi_dpsi[6];
    interp2D_eval_dB(psi_dpsi, &Bdata->psi, r, z);
    B_dB[0] = B_dB[0]+1/(2*math_pi)*-psi_dpsi[2]/r;
    B_dB[1] = B_dB[1] = 1/(2*math_pi*r)*(1/r*psi_dpsi[2]-psi_dpsi[5]);
    B_dB[3] = -1/(2*math_pi*r)*psi_dpsi[4];
    B_dB[8] = B_dB[0]+1/(2*math_pi)*psi_dpsi[1]/r;
    B_dB[9] = B_dB[1] = 1/(2*math_pi*r)*(-1/r*psi_dpsi[1]*psi_dpsi[3]);
    B_dB[11] = 1/(2*math_pi*r)*psi_dpsi[5];
    #endif
}

/**
 * @brief Interpolate magnetic field value and derivatives using bicubic 
 *        interpolation
 *
 * This function interpolates the magnetic field and it's derivatives for a
 * single component at the given point using 2D bicubic interpolation with   
 * Catmull-Rom splines and their partial derivatives.
 *
 * Different components can be evaluated with successive calls to the function, 
 * giving a pointer to the corresponding location in the full B_dB array as the 
 * second parameter (i.e. B_dB[0] for Br, B_dB[4] for Bphi, B_dB[8] for Bz).
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param B_dB_component array where magnetic field values and derivatives for 
 *             a single component will be stored (value->B_dB_component[0][i], 
 *             d/dr -> B_dB_component[1][i], d/dphi -> B_dB_component[2][i],
 *             d/dz -> B_dB_component[3][i]
 * @param t_r r parameter in the given cell, t_r = [0,1)
 * @param t_z z parameter in the given cell, t_z = [0,1)
 * @param i_r r axis index of the cell
 * @param i_z z axis index of the cell
 * @param n_r length of the data grid on r axis
 * @param r_grid r grid interval 
 * @param z_grid z grid interval
 * @param B pointer to magnetic field component data; can be B_r, B_phi or B_z
 *          from Bdata struct
 */

real B_2DS_get_axis_r(B_2DS_data* Bdata) {
    return Bdata->axis_r;
}

real B_2DS_get_axis_z(B_2DS_data* Bdata) {
    return Bdata->axis_z;
}
