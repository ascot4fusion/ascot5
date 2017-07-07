/**
 * @file B_2D.c
 * @brief 2D magnetic field with cubic interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"
#include "B_2Dspl.h"

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
void B_2Dspl_init_offload(B_2D_offload_data* offload_data, real** offload_array) {
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
void B_2Dspl_free_offload(B_2D_offload_data* offload_data, real** offload_array) {
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
void B_2Dspl_init(B_2D_data* Bdata, B_2D_offload_data* offload_data,
               real* offload_array) {
    Bdata->n_r = offload_data->n_r;
    Bdata->n_z = offload_data->n_z;
    Bdata->r_min = offload_data->r_min;
    Bdata->r_max = offload_data->r_max;
    Bdata->r_grid = offload_data->r_grid;
    Bdata->z_min = offload_data->z_min;
    Bdata->z_max = offload_data->z_max;
    Bdata->z_grid = offload_data->z_grid;
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;
    Bdata->psi = &offload_array[0];
    Bdata->B_r = &offload_array[offload_data->n_r*offload_data->n_z];
    Bdata->B_phi = &offload_array[2*offload_data->n_r*offload_data->n_z];
    Bdata->B_z = &offload_array[3*offload_data->n_r*offload_data->n_z];
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
void B_2Dspl_eval_B(real B[], real r, real phi, real z, B_2D_data* Bdata) {
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

    B[0] = B_2D_bicubic(t_r, t_z, i_r, i_z, Bdata->n_r, Bdata->B_r);
    B[1] = B_2D_bicubic(t_r, t_z, i_r, i_z, Bdata->n_r, Bdata->B_phi);
    B[2] = B_2D_bicubic(t_r, t_z, i_r, i_z, Bdata->n_r, Bdata->B_z); 
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
void B_2Dspl_eval_psi(real psi[], real r, real phi, real z, B_2D_data* Bdata)
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
void B_2Dspl_eval_rho(real rho[], real psi, B_2D_data* Bdata) {
    if(psi - Bdata->psi0 < 0) {
        rho[0] = 0;
    }
    else {
        rho[0] = sqrt((psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0));
    }
}

/**
 * @brief Interpolate magnetic field value using bicubic interpolation
 *
 * This function interpolates the magnetic field at the given point using
 * 2D bicubic interpolation with Catmull-Rom splines.
 *
 * @param t_r r parameter in the given cell, t_r = [0,1)
 * @param t_z z parameter in the given cell, t_z = [0,1)
 * @param i_r r axis index of the cell
 * @param i_z z axis index of the cell
 * @param n_r length of the data grid on r axis
 * @param B pointer to magnetic field component data; can be B_r, B_phi or B_z
 *          from Bdata struct
 */
real B_2Dspl_bicubic(real t_r, real t_z, int i_r, int i_z, int n_r, real* B) {
    /* Value of the field at the corners of the grid square to be interpolated
     * (1,1)-(2,2) and the adjacent squares. */
    real B00 = B[(i_z-1)*n_r + i_r-1];
    real B10 = B[i_z*n_r + i_r-1];
    real B20 = B[(i_z+1)*n_r + i_r-1];
    real B30 = B[(i_z+2)*n_r + i_r-1];

    real B01 = B[(i_z-1)*n_r + i_r];
    real B11 = B[i_z*n_r + i_r];
    real B21 = B[(i_z+1)*n_r + i_r];
    real B31 = B[(i_z+2)*n_r + i_r];

    real B02 = B[(i_z-1)*n_r + i_r+1];
    real B12 = B[i_z*n_r + i_r+1];
    real B22 = B[(i_z+1)*n_r + i_r+1]; 
    real B32 = B[(i_z+2)*n_r + i_r+1];

    real B03 = B[(i_z-1)*n_r + i_r+2];
    real B13 = B[i_z*n_r + i_r+2];
    real B23 = B[(i_z+1)*n_r + i_r+2];
    real B33 = B[(i_z+2)*n_r + i_r+2];

    /* Interpolate the values along the r axis */
    real t_r2 = t_r*t_r;
    real t_r3 = t_r2*t_r;
    real t_r_coeff_0 = 0.5*(-t_r3 + 2*t_r2 - t_r);
    real t_r_coeff_1 = 0.5*(3*t_r3 - 5*t_r2 + 2);
    real t_r_coeff_2 = 0.5*(-3*t_r3 + 4*t_r2 + t_r);
    real t_r_coeff_3 = 0.5*(t_r3 - t_r2);

    real B0t =   B00 * (t_r_coeff_0)
               + B01 * (t_r_coeff_1)
               + B02 * (t_r_coeff_2)
               + B03 * (t_r_coeff_3);
    real B1t =   B10 * (t_r_coeff_0)
               + B11 * (t_r_coeff_1)
               + B12 * (t_r_coeff_2)
               + B13 * (t_r_coeff_3);
    real B2t =   B20 * (t_r_coeff_0)
               + B21 * (t_r_coeff_1)
               + B22 * (t_r_coeff_2)
               + B23 * (t_r_coeff_3);
    real B3t =   B30 * (t_r_coeff_0)
               + B31 * (t_r_coeff_1)
               + B32 * (t_r_coeff_2)
               + B33 * (t_r_coeff_3);

    /* Interpolate the values along the z axis using the r axis interpolated
     * values */
    real t_z2 = t_z*t_z;
    real t_z3 = t_z2*t_z;

    real t_z_coeff_0 = 0.5*(-t_z3 + 2*t_z2 - t_z);
    real t_z_coeff_1 = 0.5*(3*t_z3 - 5*t_z2 + 2);
    real t_z_coeff_2 = 0.5*(-3*t_z3 + 4*t_z2 + t_z);
    real t_z_coeff_3 = 0.5*(t_z3 - t_z2);

    real Btt =   B0t * (t_z_coeff_0)
               + B1t * (t_z_coeff_1)
               + B2t * (t_z_coeff_2)
               + B3t * (t_z_coeff_3);

    return Btt;
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
void B_2Dspl_eval_B_dB(real B_dB[], real r, real phi, real z, B_2D_data* Bdata) {
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

    B_2D_bicubic_derivs(&B_dB[0], t_r, t_z, i_r, i_z, Bdata->n_r,
                        Bdata->r_grid, Bdata->z_grid, Bdata->B_r);
    B_2D_bicubic_derivs(&B_dB[4], t_r, t_z, i_r, i_z, Bdata->n_r,
                        Bdata->r_grid, Bdata->z_grid, Bdata->B_phi);
    B_2D_bicubic_derivs(&B_dB[8], t_r, t_z, i_r, i_z, Bdata->n_r,
                        Bdata->r_grid, Bdata->z_grid, Bdata->B_z);
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
void B_2Dsspl_bicubic_derivs(real B_dB_component[], real t_r,
                         real t_z, int i_r, int i_z, int n_r, real r_grid,
                         real z_grid, real* B) {
    /* Value of the field at the corners of the grid square to be interpolated
     * (1,1)-(2,2) and the adjacent squares. */
    real B00 = B[(i_z-1)*n_r + i_r-1];
    real B10 = B[i_z*n_r + i_r-1];
    real B20 = B[(i_z+1)*n_r + i_r-1];
    real B30 = B[(i_z+2)*n_r + i_r-1];

    real B01 = B[(i_z-1)*n_r + i_r];
    real B11 = B[i_z*n_r + i_r];
    real B21 = B[(i_z+1)*n_r + i_r];
    real B31 = B[(i_z+2)*n_r + i_r];

    real B02 = B[(i_z-1)*n_r + i_r+1];
    real B12 = B[i_z*n_r + i_r+1];
    real B22 = B[(i_z+1)*n_r + i_r+1]; 
    real B32 = B[(i_z+2)*n_r + i_r+1];

    real B03 = B[(i_z-1)*n_r + i_r+2];
    real B13 = B[i_z*n_r + i_r+2];
    real B23 = B[(i_z+1)*n_r + i_r+2];
    real B33 = B[(i_z+2)*n_r + i_r+2];

    /* Interpolate the values along the r axis */
    real t_r2 = t_r*t_r;
    real t_r3 = t_r2*t_r;
    real t_r_coeff_0 = 0.5*(-t_r3 + 2*t_r2 - t_r);
    real t_r_coeff_1 = 0.5*(3*t_r3 - 5*t_r2 + 2);
    real t_r_coeff_2 = 0.5*(-3*t_r3 + 4*t_r2 + t_r);
    real t_r_coeff_3 = 0.5*(t_r3 - t_r2);

    real B0t =   B00 * (t_r_coeff_0)
               + B01 * (t_r_coeff_1)
               + B02 * (t_r_coeff_2)
               + B03 * (t_r_coeff_3);
    real B1t =   B10 * (t_r_coeff_0)
               + B11 * (t_r_coeff_1)
               + B12 * (t_r_coeff_2)
               + B13 * (t_r_coeff_3);
    real B2t =   B20 * (t_r_coeff_0)
               + B21 * (t_r_coeff_1)
               + B22 * (t_r_coeff_2)
               + B23 * (t_r_coeff_3);
    real B3t =   B30 * (t_r_coeff_0)
               + B31 * (t_r_coeff_1)
               + B32 * (t_r_coeff_2)
               + B33 * (t_r_coeff_3);

    /* Interpolate the values along the z axis using the r axis interpolated
     * values */
    real t_z2 = t_z*t_z;
    real t_z3 = t_z2*t_z;

    real t_z_coeff_0 = 0.5*(-t_z3 + 2*t_z2 - t_z);
    real t_z_coeff_1 = 0.5*(3*t_z3 - 5*t_z2 + 2);
    real t_z_coeff_2 = 0.5*(-3*t_z3 + 4*t_z2 + t_z);
    real t_z_coeff_3 = 0.5*(t_z3 - t_z2);

    real Btt =   B0t * (t_z_coeff_0)
               + B1t * (t_z_coeff_1)
               + B2t * (t_z_coeff_2)
               + B3t * (t_z_coeff_3);

    /* d/dr coefficients */
    real dt_r_coeff_0 = 0.5*(-3*t_r2 + 4*t_r - 1);
    real dt_r_coeff_1 = 0.5*(9*t_r2 - 10*t_r);
    real dt_r_coeff_2 = 0.5*(-9*t_r2 + 8*t_r + 1);
    real dt_r_coeff_3 = 0.5*(3*t_r2 - 2*t_r);

    /* d/dr interpolated derivatives */
    real dB0t =   B00 * (dt_r_coeff_0)
                + B01 * (dt_r_coeff_1)
                + B02 * (dt_r_coeff_2)
                + B03 * (dt_r_coeff_3);
    real dB1t =   B10 * (dt_r_coeff_0)
                + B11 * (dt_r_coeff_1)
                + B12 * (dt_r_coeff_2)
                + B13 * (dt_r_coeff_3);
    real dB2t =   B20 * (dt_r_coeff_0)
                + B21 * (dt_r_coeff_1)
                + B22 * (dt_r_coeff_2)
                + B23 * (dt_r_coeff_3);
    real dB3t =   B30 * (dt_r_coeff_0)
                + B31 * (dt_r_coeff_1)
                + B32 * (dt_r_coeff_2)
                + B33 * (dt_r_coeff_3);

    /* d/dr using the interpolated derivatives */
    real dr_Btt =   dB0t * (t_z_coeff_0)
                  + dB1t * (t_z_coeff_1)
                  + dB2t * (t_z_coeff_2)
                  + dB3t * (t_z_coeff_3);

    /* d/dz coefficients */
    real dt_z_coeff_0 = 0.5*(-3*t_z2 + 4*t_z - 1);
    real dt_z_coeff_1 = 0.5*(9*t_z2 - 10*t_z);
    real dt_z_coeff_2 = 0.5*(-9*t_z2 + 8*t_z + 1);
    real dt_z_coeff_3 = 0.5*(3*t_z2 - 2*t_z);

    /* d/dz using the interpolated values */
    real dz_Btt =   B0t * (dt_z_coeff_0)
                  + B1t * (dt_z_coeff_1)
                  + B2t * (dt_z_coeff_2)
                  + B3t * (dt_z_coeff_3);

    B_dB_component[0] = Btt;
    B_dB_component[1] = dr_Btt/r_grid;
    B_dB_component[2] = 0;
    B_dB_component[3] = dz_Btt/z_grid;
}

real B_2Dspl_get_axis_r(B_2D_data* Bdata) {
    return Bdata->axis_r;
}

real B_2Dspl_get_axis_z(B_2D_data* Bdata) {
    return Bdata->axis_z;
}
