/**
 * @file B_3DS.h
 * @brief Header file for B_3DS.c
 */
#ifndef B_3DS_H
#define B_3DS_H
#include "../ascot5.h"
#include "../spline/interp2D.h" /* for 2D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */

/**
 * @brief 3D magnetic field parameters that will be offloaded to target
 */
typedef struct {
    int n_r;                    /**< number of r grid points */
    int n_z;                    /**< number of z grid points */
    int n_phi;                  /**< number of phi grid points */
    real r_min;                 /**< minimum r coordinate in the grid */
    real r_max;                 /**< maximum r coordinate in the grid */
    real r_grid;                /**< r grid interval (r_max-r_min)/(n_r-1) */
    real z_min;                 /**< minimum z coordinate in the grid */
    real z_max;                 /**< maximum z coordinate in the grid */
    real z_grid;                /**< z grid interval (z_max-z_min)/(n_z-1) */
    real phi_min;               /**< minimum phi coordinate in the grid */
    real phi_max;               /**< maximum phi coordinate in the grid */
    real phi_grid;              /**< phi grid interval 2pi/(n_phi-1) */
    real psi0;                  /**< sqrt(psi) value at magnetic axis */
    real psi1;                  /**< sqrt(psi) value at separatrix */
    real axis_r;                /**< r coordinate of magnetic axis */
    real axis_z;                /**< z coordinate of magnetic axis */
    int offload_array_length;   /**< number of elements in offload_array */
} B_3DS_offload_data;

/**
 * @brief 3D magnetic field parameters on the target
 */
typedef struct {
    real psi0;              /**< sqrt(psi) value at magnetic axis */
    real psi1;              /**< sqrt(psi) value at separatrix */
    real axis_r;            /**< r coordinate of magnetic axis */
    real axis_z;            /**< z coordinate of magnetic axis */
    interp2D_data psi;     /**< pointer to start of psi interpolation data struct */
    interp3D_data B_r;     /**< pointer to start of B_r interpolation data struct */
    interp3D_data B_phi;   /**< pointer to start of B_phi interpolation data struct */
    interp3D_data B_z;     /**< pointer to start of B_z interpolation data struct */
} B_3DS_data;

void B_3DS_init_offload(B_3DS_offload_data* offload_data, real** offload_array);
void B_3DS_free_offload(B_3DS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_3DS_init(B_3DS_data* Bdata, B_3DS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
void B_3DS_eval_B(real B[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_3DS_eval_psi(real psi[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_3DS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_3DS_eval_rho(real rho[], real psi, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_3DS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(B)
real B_3DS_tricubic(real t_r, real t_phi, real t_z, int i_r, int i_phi, int i_z, int n_z, int n_r, real* B);
#pragma omp declare simd uniform(Bdata)
void B_3DS_eval_B_dB(real B_dB[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(B)
void B_3DS_tricubic_derivs(real B_dB_component[], real t_r, real t_phi, real t_z,
                          int i_r, int i_phi, int i_z, int n_r, int n_z,
                          real r_grid, real phi_grid, real z_grid, real* B);
#pragma omp declare simd uniform(B)
real B_3DS_get_axis_r(B_3DS_data* Bdata);
#pragma omp declare simd uniform(B)
real B_3DS_get_axis_z(B_3DS_data* Bdata);
#pragma omp end declare target   
#endif
