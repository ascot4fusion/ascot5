/**
 * @file B_2DS.h
 * @brief Header file for B_2DS.c
 */
#ifndef B_2DS_H
#define B_2DS_H
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp2D.h" /* for 2D interpolation routines */

/**
 * @brief 2D magnetic field parameters that will be offloaded to target
 */
typedef struct {
    int n_r;                    /**< number of r grid points */
    int n_z;                    /**< number of z grid points */
    real r_min;                 /**< minimum r coordinate in the grid */
    real r_max;                 /**< maximum r coordinate in the grid */
    real r_grid;                /**< r grid interval (r_max-r_min)/(n_r-1) */
    real z_min;                 /**< minimum z coordinate in the grid */
    real z_max;                 /**< maximum z coordinate in the grid */
    real z_grid;                /**< z grid interval (z_max-z_min)/(n_r-1) */
    real psi0;                  /**< sqrt(psi) value at magnetic axis */
    real psi1;                  /**< sqrt(psi) value at separatrix */
    real axis_r;                /**< r coordinate of magnetic axis */
    real axis_z;                /**< z coordinate of magnetic axis */
    int offload_array_length;   /**< number of elements in offload_array */
} B_2DS_offload_data;

/**
 * @brief 2D magnetic field parameters on the target
 */
typedef struct {
    real psi0;             /**< sqrt(psi) value at magnetic axis */
    real psi1;             /**< sqrt(psi) value at separatrix */
    real axis_r;           /**< r coordinate of magnetic axis */
    real axis_z;           /**< z coordinate of magnetic axis */
    interp2D_data psi;     /**< pointer to start of psi interpolation data struct */
    interp2D_data B_r;     /**< pointer to start of B_r interpolation data struct */
    interp2D_data B_phi;   /**< pointer to start of B_phi interpolation data struct */
    interp2D_data B_z;     /**< pointer to start of B_z interpolation data struct */
} B_2DS_data;

void B_2DS_init_offload(B_2DS_offload_data* offload_data, real** offload_array);
void B_2DS_free_offload(B_2DS_offload_data* offload_data, real** offload_array);

#pragma omp declare target
int B_2DS_init(B_2DS_data* Bdata, B_2DS_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_psi(real psi[], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_psi_SIMD(int i, real psi[NSIMD], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_rho(real rho[], real psi, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_rho_SIMD(int i, real rho[NSIMD], real psi, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_rho_drho(real rho[], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_B(real B[], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_B_SIMD(int i, real B[3][NSIMD], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_B_dB(real B_dB[], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
a5err B_2DS_eval_B_dB_SIMD(int i, real B_dB[12][NSIMD], real r, real phi, real z, B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
real B_2DS_get_axis_r(B_2DS_data* Bdata);
#pragma omp declare simd uniform(Bdata) simdlen(8)
real B_2DS_get_axis_z(B_2DS_data* Bdata);
#pragma omp end declare target   
#endif
