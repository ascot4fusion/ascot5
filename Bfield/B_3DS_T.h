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
 * @brief 3D magnetic field parameters on the host
 */
typedef struct {
    int n_time;
    real time[10];

    int psigrid_n_r;            /**< number of r grid points in psi data */
    int psigrid_n_z;            /**< number of z grid points in psi data */
    int psigrid_n_phi;          /**< number of phi grid points in psi data */
    real psigrid_r_min;         /**< minimum r coordinate in the grid in psi data */
    real psigrid_r_max;         /**< maximum r coordinate in the grid in psi data */
    real psigrid_r_grid;        /**< r grid interval (r_max-r_min)/(n_r-1) in psi data */
    real psigrid_z_min;         /**< minimum z coordinate in the grid in psi data */
    real psigrid_z_max;         /**< maximum z coordinate in the grid in psi data */
    real psigrid_z_grid;        /**< z grid interval (z_max-z_min)/(n_z-1) in psi data */

    int Bgrid_n_r;              /**< number of r grid points in B data */
    int Bgrid_n_z;              /**< number of z grid points in B data */
    int n_phi;                  /**< number of phi grid points in B data */
    real Bgrid_r_min;           /**< minimum r coordinate in the grid in B data */
    real Bgrid_r_max;           /**< maximum r coordinate in the grid in B data */
    real Bgrid_r_grid;          /**< r grid interval (r_max-r_min)/(n_r-1) in B data */
    real Bgrid_z_min;           /**< minimum z coordinate in the grid in B data */
    real Bgrid_z_max;           /**< maximum z coordinate in the grid in B data */
    real Bgrid_z_grid;          /**< z grid interval (z_max-z_min)/(n_z-1) in B data */
    real phi_min;               /**< minimum phi coordinate in the grid in B data */
    real phi_max;               /**< maximum phi coordinate in the grid in B data */
    real phi_grid;              /**< phi grid interval 2pi/(n_phi-1) in B data */
    
    real psi0[10];                  /**< sqrt(psi) value at magnetic axis */
    real psi1[10];                  /**< sqrt(psi) value at separatrix */
    real axis_r[10];                /**< r coordinate of magnetic axis */
    real axis_z[10];                /**< z coordinate of magnetic axis */
    int offload_array_length;   /**< number of elements in offload_array */
} B_3DS_T_offload_data;

/**
 * @brief 3D magnetic field parameters on the target
 */
typedef struct {
    int n_time;
    real time;
    real psi0;              /**< sqrt(psi) value at magnetic axis */
    real psi1;              /**< sqrt(psi) value at separatrix */
    real axis_r;            /**< r coordinate of magnetic axis */
    real axis_z;            /**< z coordinate of magnetic axis */
    interp2D_data psi;     /**< pointer to start of psi interpolation data struct */
    interp3D_data B_r;     /**< pointer to start of B_r interpolation data struct */
    interp3D_data B_phi;   /**< pointer to start of B_phi interpolation data struct */
    interp3D_data B_z;     /**< pointer to start of B_z interpolation data struct */
} B_3DS_T_data;

void B_3DS_init_offload(B_3DS_offload_data* offload_data, real** offload_array);
void B_3DS_free_offload(B_3DS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
int B_3DS_init(B_3DS_data* Bdata, B_3DS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_psi(real psi[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd linear(i) uniform(psi, Bdata)
a5err B_3DS_eval_psi_SIMD(int i, real psi[NSIMD], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_rho(real rho[], real psi, B_3DS_data* Bdata);
#pragma omp declare simd linear(i) uniform(rho, Bdata)
a5err B_3DS_eval_rho_SIMD(int i, real rho[NSIMD], real psi, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_B(real B[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd linear(i) uniform(B, Bdata)
a5err B_3DS_eval_B_SIMD(int i, real B[3][NSIMD], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_eval_B_dB(real B_dB[], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd linear(i) uniform(B_dB, Bdata)
a5err B_3DS_eval_B_dB_SIMD(int i, real B_dB[12][NSIMD], real r, real phi, real z, B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DS_get_axis_r(B_3DS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DS_get_axis_z(B_3DS_data* Bdata);
#pragma omp end declare target   
#endif
