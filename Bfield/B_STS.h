/**
 * @file B_STS.h
 * @brief Header file for B_STS.c
 */
#ifndef B_STS_H
#define B_STS_H
#include "../ascot5.h"
#include "../linint/linint1D.h" /* for 1D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */

/**
 * @brief stellarator magnetic field parameters that will be offloaded to target
 */
typedef struct {
    int n_r;                    /**< number of r grid points */
    int n_z;                    /**< number of z grid points */
    int n_phi;                  /**< number of phi grid points */
    int periods;                /**< number of toroidal periods */
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
    int n_axis;                 /**< number of phi grid points for magnetic axis */
    real axis_min;              /**< minimum phi coordinate in the magnetic axis grid */
    real axis_max;              /**< maximum phi coordinate in the magnetic axis grid */
    real axis_grid;             /**< phi grid interval 2pi/(n_phi-1) */
    int offload_array_length;   /**< number of elements in offload_array */
} B_STS_offload_data;

/**
 * @brief stellarator magnetic field parameters on the target
 */
typedef struct {
    real psi0;                  /**< sqrt(psi) value at magnetic axis */
    real psi1;                  /**< sqrt(psi) value at separatrix */
    real periods;               /**< number of toroidal periods */
    linint1D_data axis_r;       /**< r coordinate of magnetic axis */
    linint1D_data axis_z;       /**< z coordinate of magnetic axis */
    interp3D_data psi;          /**< pointer to start of psi interpolation data struct */
    interp3D_data B_r;          /**< pointer to start of B_r interpolation data struct */
    interp3D_data B_phi;        /**< pointer to start of B_phi interpolation data struct */
    interp3D_data B_z;          /**< pointer to start of B_z interpolation data struct */
} B_STS_data;

void B_STS_init_offload(B_STS_offload_data* offload_data, real** offload_array);
void B_STS_free_offload(B_STS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
int B_STS_init(B_STS_data* Bdata, B_STS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_psi(real psi[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd linear(i) uniform(psi, Bdata)  
a5err B_STS_eval_psi_SIMD(int i, real psi[NSIMD], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_rho(real rho[], real psi, B_STS_data* Bdata);
#pragma omp declare simd linear(i) uniform(rho, Bdata)  
a5err B_STS_eval_rho_SIMD(int i, real rho[NSIMD], real psi, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_eval_B(real B[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd linear(i) uniform(B, Bdata)
a5err B_STS_eval_B_SIMD(int i, real B[3][NSIMD], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_STS_eval_B_dB(real B_dB[], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd linear(i) uniform(B_dB, Bdata)  
a5err B_STS_eval_B_dB_SIMD(int i, real B_dB[12][NSIMD], real r, real phi, real z, B_STS_data* Bdata);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_get_axis_r(real axis_r[], B_STS_data* Bdata, real phi);
#pragma omp declare simd uniform(Bdata)  
a5err B_STS_get_axis_z(real axis_r[], B_STS_data* Bdata, real phi);
#pragma omp end declare target   
#endif
