/**
 * @file B_3D.c
 * @brief 3D magnetic field with cubic interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "math.h"
#include "ascot5.h"
#include "B_3D.h"
#include "B_2D.h" /* for 2D interpolation routines */
#include "B_GS.h"

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
void B_3D_init_offload(B_3D_offload_data* offload_data, real** offload_array) {
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
 * @brief Create 3D field from BGS field for testing
 */
void B_3D_init_offload_dummy(B_3D_offload_data* offload_data, real** offload_array) {
    offload_data->r_min = 3;
    offload_data->r_max = 9;
    offload_data->n_r = 100;
    offload_data->z_min = -5;
    offload_data->z_max = 5;
    offload_data->n_z = 150;
    offload_data->n_phi = 360;

    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
                               / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
                                   / (offload_data->n_z - 1);
    offload_data->phi_grid = 2*math_pi / (offload_data->n_phi - 1);

    int psi_size = offload_data->n_r*offload_data->n_z;
    int B_size = offload_data->n_r*offload_data->n_z*(offload_data->n_phi+2);
    *offload_array = (real*) malloc((psi_size + 3 * B_size) * sizeof(real));
    offload_data->offload_array_length = psi_size + 3 * B_size;

    B_GS_offload_data B_offload_data;
    real* B_offload_array;
    B_GS_init_offload(&B_offload_data, &B_offload_array);
    B_GS_data B_data;
    B_GS_init(&B_data, &B_offload_data, B_offload_array);

    int i, j, k;

    for(i = 1; i < (offload_data->n_phi+1); i++) {
        for(j = 0; j < offload_data->n_z; j++) {
            for(k = 0; k < offload_data->n_r; k++) {
                real r = offload_data->r_min + k * offload_data->r_grid;
                real z = offload_data->z_min + j * offload_data->z_grid;
                real phi = (i-1) * offload_data->phi_grid;
                real B[3];
                B_GS_eval_B(B, r, phi, z, &B_data);
                (*offload_array)[i*offload_data->n_z*offload_data->n_r
                              + j*offload_data->n_r + k] = B[0];
                (*offload_array)[i*offload_data->n_z*offload_data->n_r
                              + j*offload_data->n_r + k + B_size] = B[1];
                (*offload_array)[i*offload_data->n_z*offload_data->n_r
                              + j*offload_data->n_r + k + 2*B_size] = B[2];
            }
        }
    }

    memcpy(&(*offload_array)[0],
        &(*offload_array)[(offload_data->n_phi-2)*offload_data->n_r*offload_data->n_z],
        offload_data->n_r * offload_data->n_z * sizeof(real));
    memcpy(&(*offload_array)[B_size],
 &(*offload_array)[(offload_data->n_phi-2)*offload_data->n_r*offload_data->n_z+B_size],
           offload_data->n_r * offload_data->n_z * sizeof(real));
    memcpy(&(*offload_array)[2*B_size],
 &(*offload_array)[(offload_data->n_phi-2)*offload_data->n_r*offload_data->n_z+2*B_size],
           offload_data->n_r * offload_data->n_z * sizeof(real));

    memcpy(&(*offload_array)[(offload_data->n_phi)*offload_data->n_r*offload_data->n_z],
        &(*offload_array)[offload_data->n_r*offload_data->n_z],
        offload_data->n_r * offload_data->n_z * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi)*offload_data->n_r*offload_data->n_z+B_size],
        &(*offload_array)[offload_data->n_r*offload_data->n_z+B_size],
        offload_data->n_r * offload_data->n_z * sizeof(real));
    memcpy(&(*offload_array)[(offload_data->n_phi)*offload_data->n_r*offload_data->n_z+2*B_size],
        &(*offload_array)[offload_data->n_r*offload_data->n_z+2*B_size],
        offload_data->n_r * offload_data->n_z * sizeof(real));

    offload_data->psi0 = -0.0365;
    offload_data->psi1 = 0.0;

    for(j = 0; j < offload_data->n_z; j++) {
        for(k = 0; k < offload_data->n_r; k++) {
            real r = offload_data->r_min + k * offload_data->r_grid;
            real z = offload_data->z_min + j * offload_data->z_grid;
            real psi[1];
            B_GS_eval_psi(psi, r, 0.0, z, &B_data);
            (*offload_array)[j*offload_data->n_r + k + 3*B_size] = psi[0];
        }
    }
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3D_free_offload(B_3D_offload_data* offload_data, real** offload_array) {
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
void B_3D_init(B_3D_data* Bdata, B_3D_offload_data* offload_data,
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
 * bicubic interpolation on the 3D magnetic field data. This is a SIMD
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
void B_3D_eval_B(real B[], real r, real phi, real z, B_3D_data* Bdata) {
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
 * bicubic interpolation on the 3D magnetic field data. This is a SIMD
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
void B_3D_eval_psi(real psi[], real r, real phi, real z,
                   B_3D_data* Bdata)
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
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * This function evaluates the poloidal flux psi at the given coordinates using
 * tricubic interpolation on the 3D magnetic field data. This is a SIMD
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
 * @todo Change to a scalar elemental function and compare performance
 */
void B_3D_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z,
                   B_3D_data* Bdata)
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

    B_2D_bicubic_derivs(psi_dpsi, t_r, t_z, i_r, i_z, Bdata->n_r,
                        Bdata->r_grid, Bdata->z_grid, Bdata->psi);
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
void B_3D_eval_rho(real rho[], real psi, B_3D_data* Bdata) {
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
 * at the given coordinates using tricubic interpolation on the 
 * 3D magnetic field data. This is a SIMD function, so the values are
 * placed in an NSIMD length struct.
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
void B_3D_eval_rho_drho(real rho_drho[], real r, real phi, real z,
                    B_3D_data* Bdata) {
    real rho;
    B_3D_eval_psi_dpsi(rho_drho, r, phi, z, Bdata);
    /* Convert: rho = sqrt(psi), drho = dpsi/(2 * sqrt(psi))
     * Note that rho_drho[2] = 1/R * drho/dphi, because of cylindrical gradient
     */
    rho = sqrt(rho_drho[0]);
    rho_drho[0] = rho;
    rho_drho[1] = rho_drho[1] / (2*rho);
    rho_drho[2] = rho_drho[2] / (2*rho);
    rho_drho[3] = rho_drho[3] / (2*rho);
}

real B_3D_tricubic(real t_r, real t_phi, real t_z, int i_r, int i_phi, int i_z, int n_z, int n_r, real* B) {
    real B000 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B010 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r-1];
    real B020 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B030 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B001 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B011 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r];
    real B021 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B031 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B002 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B012 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r+1];
    real B022 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B032 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B003 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B013 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r+2];
    real B023 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B033 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real B100 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B110 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r-1];
    real B120 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B130 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B101 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B111 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r];
    real B121 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B131 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B102 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B112 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r+1];
    real B122 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B132 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B103 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B113 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r+2];
    real B123 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B133 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real B200 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B210 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r-1];
    real B220 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B230 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B201 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B211 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r];
    real B221 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B231 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B202 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B212 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r+1];
    real B222 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B232 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B203 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B213 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r+2];
    real B223 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B233 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real B300 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B310 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r-1];
    real B320 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B330 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B301 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B311 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r];
    real B321 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B331 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B302 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B312 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r+1];
    real B322 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B332 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B303 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B313 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r+2];
    real B323 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B333 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real t_phi2 = t_phi*t_phi;
    real t_phi3 = t_phi2*t_phi;

    real t_phi_coeff_0 = 0.5*(-t_phi3 + 2*t_phi2 - t_phi);
    real t_phi_coeff_1 = 0.5*(3*t_phi3 - 5*t_phi2 + 2);
    real t_phi_coeff_2 = 0.5*(-3*t_phi3 + 4*t_phi2 + t_phi);
    real t_phi_coeff_3 = 0.5*(t_phi3 - t_phi2);
    
    real t_r2 = t_r*t_r;
    real t_r3 = t_r2*t_r;
    real t_r_coeff_0 = 0.5*(-t_r3 + 2*t_r2 - t_r);
    real t_r_coeff_1 = 0.5*(3*t_r3 - 5*t_r2 + 2);
    real t_r_coeff_2 = 0.5*(-3*t_r3 + 4*t_r2 + t_r);
    real t_r_coeff_3 = 0.5*(t_r3 - t_r2);

    real t_z2 = t_z*t_z;
    real t_z3 = t_z2*t_z;

    real t_z_coeff_0 = 0.5*(-t_z3 + 2*t_z2 - t_z);
    real t_z_coeff_1 = 0.5*(3*t_z3 - 5*t_z2 + 2);
    real t_z_coeff_2 = 0.5*(-3*t_z3 + 4*t_z2 + t_z);
    real t_z_coeff_3 = 0.5*(t_z3 - t_z2);

    real B00t =  B000 * (t_r_coeff_0)
                + B001 * (t_r_coeff_1)
                + B002 * (t_r_coeff_2)
                + B003 * (t_r_coeff_3);
    real B01t =  B010 * (t_r_coeff_0)
                + B011 * (t_r_coeff_1)
                + B012 * (t_r_coeff_2)
                + B013 * (t_r_coeff_3);
    real B02t =  B020 * (t_r_coeff_0)
                + B021 * (t_r_coeff_1)
                + B022 * (t_r_coeff_2)
                + B023 * (t_r_coeff_3);
    real B03t =  B030 * (t_r_coeff_0)
                + B031 * (t_r_coeff_1)
                + B032 * (t_r_coeff_2)
                + B033 * (t_r_coeff_3);
    real B0tt =  B00t * (t_z_coeff_0)
                + B01t * (t_z_coeff_1)
                + B02t * (t_z_coeff_2)
                + B03t * (t_z_coeff_3);

    real B10t =  B100 * (t_r_coeff_0)
                + B101 * (t_r_coeff_1)
                + B102 * (t_r_coeff_2)
                + B103 * (t_r_coeff_3);
    real B11t =  B110 * (t_r_coeff_0)
                + B111 * (t_r_coeff_1)
                + B112 * (t_r_coeff_2)
                + B113 * (t_r_coeff_3);
    real B12t =  B120 * (t_r_coeff_0)
                + B121 * (t_r_coeff_1)
                + B122 * (t_r_coeff_2)
                + B123 * (t_r_coeff_3);
    real B13t =  B130 * (t_r_coeff_0)
                + B131 * (t_r_coeff_1)
                + B132 * (t_r_coeff_2)
                + B133 * (t_r_coeff_3);
    real B1tt =  B10t * (t_z_coeff_0)
                + B11t * (t_z_coeff_1)
                + B12t * (t_z_coeff_2)
                + B13t * (t_z_coeff_3);

    real B20t =  B200 * (t_r_coeff_0)
                + B201 * (t_r_coeff_1)
                + B202 * (t_r_coeff_2)
                + B203 * (t_r_coeff_3);
    real B21t =  B210 * (t_r_coeff_0)
                + B211 * (t_r_coeff_1)
                + B212 * (t_r_coeff_2)
                + B213 * (t_r_coeff_3);
    real B22t =  B220 * (t_r_coeff_0)
                + B221 * (t_r_coeff_1)
                + B222 * (t_r_coeff_2)
                + B223 * (t_r_coeff_3);
    real B23t =  B230 * (t_r_coeff_0)
                + B231 * (t_r_coeff_1)
                + B232 * (t_r_coeff_2)
                + B233 * (t_r_coeff_3);
    real B2tt =  B20t * (t_z_coeff_0)
                + B21t * (t_z_coeff_1)
                + B22t * (t_z_coeff_2)
                + B23t * (t_z_coeff_3);

    real B30t =  B300 * (t_r_coeff_0)
                + B301 * (t_r_coeff_1)
                + B302 * (t_r_coeff_2)
                + B303 * (t_r_coeff_3);
    real B31t =  B310 * (t_r_coeff_0)
                + B311 * (t_r_coeff_1)
                + B312 * (t_r_coeff_2)
                + B313 * (t_r_coeff_3);
    real B32t =  B320 * (t_r_coeff_0)
                + B321 * (t_r_coeff_1)
                + B322 * (t_r_coeff_2)
                + B323 * (t_r_coeff_3);
    real B33t =  B330 * (t_r_coeff_0)
                + B331 * (t_r_coeff_1)
                + B332 * (t_r_coeff_2)
                + B333 * (t_r_coeff_3);
    real B3tt =  B30t * (t_z_coeff_0)
                + B31t * (t_z_coeff_1)
                + B32t * (t_z_coeff_2)
                + B33t * (t_z_coeff_3);

    real Bttt =  B0tt * (t_phi_coeff_0)
                + B1tt * (t_phi_coeff_1)
                + B2tt * (t_phi_coeff_2)
                + B3tt * (t_phi_coeff_3);

    return Bttt;
}

/**
 * @brief Evaluate magnetic field and derivatives
 *
 * This function evaluates the magnetic field and it's derivatives at the given 
 * coordinates using bicubic interpolation on the 3D magnetic field data. This 
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
void B_3D_eval_B_dB(real B_dB[], real r, real phi, real z, B_3D_data* Bdata) {
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

void B_3D_tricubic_derivs(real B_dB_component[], real t_r,
                         real t_phi, real t_z, int i_r, int i_phi, int i_z,
                         int n_r, int n_z, real r_grid, real phi_grid,
                         real z_grid, real* B) {
    real B000 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B010 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r-1];
    real B020 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B030 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B001 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B011 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r];
    real B021 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B031 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B002 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B012 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r+1];
    real B022 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B032 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B003 = B[(i_phi-1)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B013 = B[(i_phi-1)*n_z*n_r + i_z*n_r + i_r+2];
    real B023 = B[(i_phi-1)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B033 = B[(i_phi-1)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real B100 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B110 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r-1];
    real B120 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B130 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B101 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B111 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r];
    real B121 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B131 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B102 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B112 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r+1];
    real B122 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B132 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B103 = B[(i_phi)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B113 = B[(i_phi)*n_z*n_r + i_z*n_r + i_r+2];
    real B123 = B[(i_phi)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B133 = B[(i_phi)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real B200 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B210 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r-1];
    real B220 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B230 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B201 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B211 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r];
    real B221 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B231 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B202 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B212 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r+1];
    real B222 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B232 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B203 = B[(i_phi+1)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B213 = B[(i_phi+1)*n_z*n_r + i_z*n_r + i_r+2];
    real B223 = B[(i_phi+1)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B233 = B[(i_phi+1)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real B300 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r-1];
    real B310 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r-1];
    real B320 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r-1];
    real B330 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r-1];

    real B301 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r];
    real B311 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r];
    real B321 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r];
    real B331 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r];

    real B302 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r+1];
    real B312 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r+1];
    real B322 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r+1]; 
    real B332 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r+1];

    real B303 = B[(i_phi+2)*n_z*n_r + (i_z-1)*n_r + i_r+2];
    real B313 = B[(i_phi+2)*n_z*n_r + i_z*n_r + i_r+2];
    real B323 = B[(i_phi+2)*n_z*n_r + (i_z+1)*n_r + i_r+2];
    real B333 = B[(i_phi+2)*n_z*n_r + (i_z+2)*n_r + i_r+2];

    real t_phi2 = t_phi*t_phi;
    real t_phi3 = t_phi2*t_phi;

    real t_phi_coeff_0 = 0.5*(-t_phi3 + 2*t_phi2 - t_phi);
    real t_phi_coeff_1 = 0.5*(3*t_phi3 - 5*t_phi2 + 2);
    real t_phi_coeff_2 = 0.5*(-3*t_phi3 + 4*t_phi2 + t_phi);
    real t_phi_coeff_3 = 0.5*(t_phi3 - t_phi2);

    real t_r2 = t_r*t_r;
    real t_r3 = t_r2*t_r;
    real t_r_coeff_0 = 0.5*(-t_r3 + 2*t_r2 - t_r);
    real t_r_coeff_1 = 0.5*(3*t_r3 - 5*t_r2 + 2);
    real t_r_coeff_2 = 0.5*(-3*t_r3 + 4*t_r2 + t_r);
    real t_r_coeff_3 = 0.5*(t_r3 - t_r2);

    real t_z2 = t_z*t_z;
    real t_z3 = t_z2*t_z;

    real t_z_coeff_0 = 0.5*(-t_z3 + 2*t_z2 - t_z);
    real t_z_coeff_1 = 0.5*(3*t_z3 - 5*t_z2 + 2);
    real t_z_coeff_2 = 0.5*(-3*t_z3 + 4*t_z2 + t_z);
    real t_z_coeff_3 = 0.5*(t_z3 - t_z2);

    real B00t =  B000 * (t_r_coeff_0)
                + B001 * (t_r_coeff_1)
                + B002 * (t_r_coeff_2)
                + B003 * (t_r_coeff_3);
    real B01t =  B010 * (t_r_coeff_0)
                + B011 * (t_r_coeff_1)
                + B012 * (t_r_coeff_2)
                + B013 * (t_r_coeff_3);
    real B02t =  B020 * (t_r_coeff_0)
                + B021 * (t_r_coeff_1)
                + B022 * (t_r_coeff_2)
                + B023 * (t_r_coeff_3);
    real B03t =  B030 * (t_r_coeff_0)
                + B031 * (t_r_coeff_1)
                + B032 * (t_r_coeff_2)
                + B033 * (t_r_coeff_3);
    real B0tt =  B00t * (t_z_coeff_0)
                + B01t * (t_z_coeff_1)
                + B02t * (t_z_coeff_2)
                + B03t * (t_z_coeff_3);

    real B10t =  B100 * (t_r_coeff_0)
                + B101 * (t_r_coeff_1)
                + B102 * (t_r_coeff_2)
                + B103 * (t_r_coeff_3);
    real B11t =  B110 * (t_r_coeff_0)
                + B111 * (t_r_coeff_1)
                + B112 * (t_r_coeff_2)
                + B113 * (t_r_coeff_3);
    real B12t =  B120 * (t_r_coeff_0)
                + B121 * (t_r_coeff_1)
                + B122 * (t_r_coeff_2)
                + B123 * (t_r_coeff_3);
    real B13t =  B130 * (t_r_coeff_0)
                + B131 * (t_r_coeff_1)
                + B132 * (t_r_coeff_2)
                + B133 * (t_r_coeff_3);
    real B1tt =  B10t * (t_z_coeff_0)
                + B11t * (t_z_coeff_1)
                + B12t * (t_z_coeff_2)
                + B13t * (t_z_coeff_3);

    real B20t =  B200 * (t_r_coeff_0)
                + B201 * (t_r_coeff_1)
                + B202 * (t_r_coeff_2)
                + B203 * (t_r_coeff_3);
    real B21t =  B210 * (t_r_coeff_0)
                + B211 * (t_r_coeff_1)
                + B212 * (t_r_coeff_2)
                + B213 * (t_r_coeff_3);
    real B22t =  B220 * (t_r_coeff_0)
                + B221 * (t_r_coeff_1)
                + B222 * (t_r_coeff_2)
                + B223 * (t_r_coeff_3);
    real B23t =  B230 * (t_r_coeff_0)
                + B231 * (t_r_coeff_1)
                + B232 * (t_r_coeff_2)
                + B233 * (t_r_coeff_3);
    real B2tt =  B20t * (t_z_coeff_0)
                + B21t * (t_z_coeff_1)
                + B22t * (t_z_coeff_2)
                + B23t * (t_z_coeff_3);

    real B30t =  B300 * (t_r_coeff_0)
                + B301 * (t_r_coeff_1)
                + B302 * (t_r_coeff_2)
                + B303 * (t_r_coeff_3);
    real B31t =  B310 * (t_r_coeff_0)
                + B311 * (t_r_coeff_1)
                + B312 * (t_r_coeff_2)
                + B313 * (t_r_coeff_3);
    real B32t =  B320 * (t_r_coeff_0)
                + B321 * (t_r_coeff_1)
                + B322 * (t_r_coeff_2)
                + B323 * (t_r_coeff_3);
    real B33t =  B330 * (t_r_coeff_0)
                + B331 * (t_r_coeff_1)
                + B332 * (t_r_coeff_2)
                + B333 * (t_r_coeff_3);
    real B3tt =  B30t * (t_z_coeff_0)
                + B31t * (t_z_coeff_1)
                + B32t * (t_z_coeff_2)
                + B33t * (t_z_coeff_3);

    real Bttt =  B0tt * (t_phi_coeff_0)
                + B1tt * (t_phi_coeff_1)
                + B2tt * (t_phi_coeff_2)
                + B3tt * (t_phi_coeff_3);

    /* d/dr */
    real dt_r_coeff_0 = 0.5*(-3*t_r2 + 4*t_r - 1);
    real dt_r_coeff_1 = 0.5*(9*t_r2 - 10*t_r);
    real dt_r_coeff_2 = 0.5*(-9*t_r2 + 8*t_r + 1);
    real dt_r_coeff_3 = 0.5*(3*t_r2 - 2*t_r);

    real dr_B00t = B000 * (dt_r_coeff_0)
                 + B001 * (dt_r_coeff_1)
                 + B002 * (dt_r_coeff_2)
                 + B003 * (dt_r_coeff_3);
    real dr_B01t = B010 * (dt_r_coeff_0)
                 + B011 * (dt_r_coeff_1)
                 + B012 * (dt_r_coeff_2)
                 + B013 * (dt_r_coeff_3);
    real dr_B02t = B020 * (dt_r_coeff_0)
                 + B021 * (dt_r_coeff_1)
                 + B022 * (dt_r_coeff_2)
                 + B023 * (dt_r_coeff_3);
    real dr_B03t = B030 * (dt_r_coeff_0)
                 + B031 * (dt_r_coeff_1)
                 + B032 * (dt_r_coeff_2)
                 + B033 * (dt_r_coeff_3);
    real dr_B0tt = dr_B00t * (t_z_coeff_0)
                 + dr_B01t * (t_z_coeff_1)
                 + dr_B02t * (t_z_coeff_2)
                 + dr_B03t * (t_z_coeff_3);

    real dr_B10t = B100 * (dt_r_coeff_0)
                 + B101 * (dt_r_coeff_1)
                 + B102 * (dt_r_coeff_2)
                 + B103 * (dt_r_coeff_3);
    real dr_B11t = B110 * (dt_r_coeff_0)
                 + B111 * (dt_r_coeff_1)
                 + B112 * (dt_r_coeff_2)
                 + B113 * (dt_r_coeff_3);
    real dr_B12t = B120 * (dt_r_coeff_0)
                 + B121 * (dt_r_coeff_1)
                 + B122 * (dt_r_coeff_2)
                 + B123 * (dt_r_coeff_3);
    real dr_B13t = B130 * (dt_r_coeff_0)
                 + B131 * (dt_r_coeff_1)
                 + B132 * (dt_r_coeff_2)
                 + B133 * (dt_r_coeff_3);
    real dr_B1tt = dr_B10t * (t_z_coeff_0)
                 + dr_B11t * (t_z_coeff_1)
                 + dr_B12t * (t_z_coeff_2)
                 + dr_B13t * (t_z_coeff_3);

    real dr_B20t = B200 * (dt_r_coeff_0)
                 + B201 * (dt_r_coeff_1)
                 + B202 * (dt_r_coeff_2)
                 + B203 * (dt_r_coeff_3);
    real dr_B21t = B210 * (dt_r_coeff_0)
                 + B211 * (dt_r_coeff_1)
                 + B212 * (dt_r_coeff_2)
                 + B213 * (dt_r_coeff_3);
    real dr_B22t = B220 * (dt_r_coeff_0)
                 + B221 * (dt_r_coeff_1)
                 + B222 * (dt_r_coeff_2)
                 + B223 * (dt_r_coeff_3);
    real dr_B23t = B230 * (dt_r_coeff_0)
                 + B231 * (dt_r_coeff_1)
                 + B232 * (dt_r_coeff_2)
                 + B233 * (dt_r_coeff_3);
    real dr_B2tt = dr_B20t * (t_z_coeff_0)
                 + dr_B21t * (t_z_coeff_1)
                 + dr_B22t * (t_z_coeff_2)
                 + dr_B23t * (t_z_coeff_3);

    real dr_B30t = B300 * (dt_r_coeff_0)
                 + B301 * (dt_r_coeff_1)
                 + B302 * (dt_r_coeff_2)
                 + B303 * (dt_r_coeff_3);
    real dr_B31t = B310 * (dt_r_coeff_0)
                 + B311 * (dt_r_coeff_1)
                 + B312 * (dt_r_coeff_2)
                 + B313 * (dt_r_coeff_3);
    real dr_B32t = B320 * (dt_r_coeff_0)
                 + B321 * (dt_r_coeff_1)
                 + B322 * (dt_r_coeff_2)
                 + B323 * (dt_r_coeff_3);
    real dr_B33t = B330 * (dt_r_coeff_0)
                 + B331 * (dt_r_coeff_1)
                 + B332 * (dt_r_coeff_2)
                 + B333 * (dt_r_coeff_3);
    real dr_B3tt = dr_B30t * (t_z_coeff_0)
                 + dr_B31t * (t_z_coeff_1)
                 + dr_B32t * (t_z_coeff_2)
                 + dr_B33t * (t_z_coeff_3);

    real dr_Bttt =  dr_B0tt * (t_phi_coeff_0)
                  + dr_B1tt * (t_phi_coeff_1)
                  + dr_B2tt * (t_phi_coeff_2)
                  + dr_B3tt * (t_phi_coeff_3);

    /* d/dz */
    real dt_z_coeff_0 = 0.5*(-3*t_z2 + 4*t_z - 1);
    real dt_z_coeff_1 = 0.5*(9*t_z2 - 10*t_z);
    real dt_z_coeff_2 = 0.5*(-9*t_z2 + 8*t_z + 1);
    real dt_z_coeff_3 = 0.5*(3*t_z2 - 2*t_z);

    /* d/dz using the interpolated values */
    real dz_B0tt =   B00t * (dt_z_coeff_0)
                   + B01t * (dt_z_coeff_1)
                   + B02t * (dt_z_coeff_2)
                   + B03t * (dt_z_coeff_3);
    real dz_B1tt =   B10t * (dt_z_coeff_0)
                   + B11t * (dt_z_coeff_1)
                   + B12t * (dt_z_coeff_2)
                   + B13t * (dt_z_coeff_3);
    real dz_B2tt =   B20t * (dt_z_coeff_0)
                   + B21t * (dt_z_coeff_1)
                   + B22t * (dt_z_coeff_2)
                   + B23t * (dt_z_coeff_3);
    real dz_B3tt =   B30t * (dt_z_coeff_0)
                   + B31t * (dt_z_coeff_1)
                   + B32t * (dt_z_coeff_2)
                   + B33t * (dt_z_coeff_3);

    real dz_Bttt =  dz_B0tt * (t_phi_coeff_0)
                  + dz_B1tt * (t_phi_coeff_1)
                  + dz_B2tt * (t_phi_coeff_2)
                  + dz_B3tt * (t_phi_coeff_3);

    /* d/dphi */
    real dt_phi_coeff_0 = 0.5*(-3*t_phi2 + 4*t_phi - 1);
    real dt_phi_coeff_1 = 0.5*(9*t_phi2 - 10*t_phi);
    real dt_phi_coeff_2 = 0.5*(-9*t_phi2 + 8*t_phi + 1);
    real dt_phi_coeff_3 = 0.5*(3*t_phi2 - 2*t_phi);

    real dphi_Bttt =  B0tt * (dt_phi_coeff_0)
                    + B1tt * (dt_phi_coeff_1)
                    + B2tt * (dt_phi_coeff_2)
                    + B3tt * (dt_phi_coeff_3);


    B_dB_component[0] = Bttt;
    B_dB_component[1] = dr_Bttt/r_grid;
    B_dB_component[2] = dphi_Bttt/phi_grid;
    B_dB_component[3] = dz_Bttt/z_grid;
}

real B_3D_get_axis_r(B_3D_data* Bdata) {
    return Bdata->axis_r;
}

real B_3D_get_axis_z(B_3D_data* Bdata) {
    return Bdata->axis_z;
}
