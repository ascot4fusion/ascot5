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
//TODO Rewrite all of this for stellarator HDF5 files
void B_ST_init_offload(B_ST_offload_data* offload_data, real** offload_array) {
    FILE* f = fopen("input.magn_bkg", "r");

    /* Read phi parameters */
    fscanf(f, "%*d %*d %d %*d %*d", &(offload_data->n_phi));

    /* Read r and z parameters */
    fscanf(f, "%lf %lf %d", &(offload_data->r_min), &(offload_data->r_max),
                            &(offload_data->n_r));
    fscanf(f, "%lf %lf %d", &(offload_data->z_min), &(offload_data->z_max),
                            &(offload_data->n_z));

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                           / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                           / (offload_data->n_z - 1);
    offload_data->phi_grid = 2*math_pi / (offload_data->n_phi);

    /* phi array starts from -0.5*phi_grid and ends at 0.5*phi_grid + 2pi! */
    offload_data->phi_min = -1.5 * offload_data->phi_grid;
    offload_data->phi_max = 1.5 * offload_data->phi_grid + 2*math_pi;

    /* Allocate offload_array */
    int psi_size = offload_data->n_r*offload_data->n_z;
    int B_size = offload_data->n_r*offload_data->n_z*(offload_data->n_phi+4);
    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    /* Skip phimaps */
    fscanf(f, "%*f %*f");

    /* Read psi */
    int i;
    for(i = 0; i < psi_size; i++) {
        fscanf(f, "%lf", &(*offload_array)[i + 3*B_size]);
    }

    /* Calculate 2D components of poloidal field */
    real* eq_B_r = (real*) malloc(psi_size*sizeof(real));
    real* eq_B_z = (real*) malloc(psi_size*sizeof(real));

    for(i = 0; i < offload_data->n_r; i++) {
        int j;
        for(j = 0; j < offload_data->n_z; j++) {
            real psi[4];
            B_2D_bicubic_derivs(psi, 0, 0, i, j, offload_data->n_r, offload_data->r_grid, offload_data->z_grid, &(*offload_array)[3*B_size]);
            eq_B_r[j*offload_data->n_r + i] = 1/(2*math_pi) * -psi[3] / (offload_data->r_min + i*offload_data->r_grid);
            eq_B_z[j*offload_data->n_r + i] = 1/(2*math_pi) * psi[1] / (offload_data->r_min + i*offload_data->r_grid);
        }
    }

    real* temp_B_r = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_phi = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));
    real* temp_B_z = (real*) malloc(psi_size*offload_data->n_phi*sizeof(real));

    /* Read B_r */
    for(i = 0; i < psi_size*offload_data->n_phi; i++) {
        fscanf(f, "%lf", &temp_B_r[i]);
    }

    /* Read B_phi */
    for(i = 0; i < psi_size*offload_data->n_phi; i++) {
        fscanf(f, "%lf", &temp_B_phi[i]);
    }

    /* Read B_z */
    for(i = 0; i < psi_size*offload_data->n_phi; i++) {
        fscanf(f, "%lf", &temp_B_z[i]);
    }

    /* permute the phi and z dimensions */
    for(i = 0; i < offload_data->n_phi; i++) {
        int j;
        for(j = 0; j < offload_data->n_z; j++) {
            int k;
            for(k = 0; k < offload_data->n_r; k++) {
               (*offload_array)[2*psi_size+i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k] =
                        eq_B_r[j*offload_data->n_r + k] +
                        temp_B_r[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
               (*offload_array)[2*psi_size+i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k + B_size] =
                        temp_B_phi[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
               (*offload_array)[2*psi_size+i*offload_data->n_z*offload_data->n_r
                          + j*offload_data->n_r + k + 2*B_size] =
                        eq_B_z[j*offload_data->n_r + k] +
                        temp_B_z[j*offload_data->n_phi*offload_data->n_r
                          + i*offload_data->n_r + k];
            }
        }
    }
    
    free(eq_B_r);
    free(eq_B_z);
    fclose(f);

    /* Copy two phi slices into opposite ends for each field
     * component */
    memcpy(&(*offload_array)[0],
        &(*offload_array)[offload_data->n_phi*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[B_size],
        &(*offload_array)[offload_data->n_phi*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[2*B_size],
        &(*offload_array)[offload_data->n_phi*psi_size+2*B_size],
        psi_size * sizeof(real));

    memcpy(&(*offload_array)[psi_size],
        &(*offload_array)[(offload_data->n_phi+1)*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[psi_size + B_size],
        &(*offload_array)[(offload_data->n_phi+1)*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[psi_size + 2*B_size],
        &(*offload_array)[(offload_data->n_phi+1)*psi_size+2*B_size],
        psi_size * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi+2)*psi_size],
        &(*offload_array)[2*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+2)*psi_size+B_size],
        &(*offload_array)[2*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+2)*psi_size+2*B_size],
        &(*offload_array)[2*psi_size+2*B_size],
        psi_size * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi+3)*psi_size],
        &(*offload_array)[3*psi_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+3)*psi_size+B_size],
        &(*offload_array)[3*psi_size+B_size],
        psi_size * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi+3)*psi_size+2*B_size],
        &(*offload_array)[3*psi_size+2*B_size],
        psi_size * sizeof(real));

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
    Bdata->r_min = offload_data->r_min;
    Bdata->r_max = offload_data->r_max;
    Bdata->r_grid = offload_data->r_grid;
    Bdata->z_min = offload_data->z_min;
    Bdata->z_max = offload_data->z_max;
    Bdata->z_grid = offload_data->z_grid;
    Bdata->phi_min = offload_data->phi_min;
    Bdata->phi_max = offload_data->phi_max;
    Bdata->phi_grid = offload_data->phi_grid;
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;
    Bdata->B_r = &offload_array[0];
    Bdata->B_phi = &offload_array[(offload_data->n_phi+4)*offload_data->n_r*offload_data->n_z];
    Bdata->B_z = &offload_array[2*(offload_data->n_phi+4)*offload_data->n_r*offload_data->n_z];
    Bdata->psi = &offload_array[3*(offload_data->n_phi+4)*offload_data->n_r*offload_data->n_z];
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
//TODO
void B_ST_eval_B(real B[], real r, real phi, real z, B_ST_data* Bdata) {
    phi = fmod(phi, 2*math_pi);
    if(phi < 0) {
        phi += 2*math_pi;
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
    real t_phi = (phi - (i_phi * Bdata->phi_grid)) / Bdata->phi_grid;
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
    int i_r = (int) floor((r - Bdata->r_min)
                    / ((Bdata->r_max - Bdata->r_min)
                       / (Bdata->n_r-1)));
    if(i_r < 0 || i_r >= Bdata->n_r) {
        i_r = 0;
    }

    int i_z = (int) floor((z - Bdata->z_min)
                    / ((Bdata->z_max - Bdata->z_min)
                       / (Bdata->n_z-1)));
    if(i_z < 0 || i_z >= Bdata->n_z) {
        i_z = 0;
    }

    real t_r = (r - (Bdata->r_min + i_r * Bdata->r_grid)) / Bdata->r_grid;
    real t_z = (z - (Bdata->z_min + i_z * Bdata->z_grid)) / Bdata->z_grid;

    psi[0] = B_2D_bicubic(t_r, t_z, i_r, i_z, Bdata->n_r, Bdata->psi);
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
//TODO 3D interpolation of rho using S
void B_ST_eval_rho(real rho[], real psi, B_ST_data* Bdata) {
    if(psi - Bdata->psi0 < 0) {
        rho[0] = 0;
    }
    else {
        rho[0] = sqrt((psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0));
    }
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
//TODO
void B_ST_eval_B_dB(real B_dB[], real r, real phi, real z, B_ST_data* Bdata) {
    phi = fmod(phi, 2*math_pi);
    if(phi < 0) {
        phi += 2*math_pi;
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
    real t_phi = (phi - (i_phi * Bdata->phi_grid)) / Bdata->phi_grid;
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
