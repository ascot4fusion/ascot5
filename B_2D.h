/**
 * @file B_2D.h
 * @brief Header file for B_2D.c
 */
#ifndef B_2D_H
#define B_2D_H
#include "ascot5.h"

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
} B_2D_offload_data;

/**
 * @brief 2D magnetic field parameters on the target
 */
typedef struct {
    int n_r;        /**< number of r grid points */
    int n_z;        /**< number of z grid points */
    real r_min;     /**< minimum r coordinate in the grid */
    real r_max;     /**< maximum r coordinate in the grid */
    real r_grid;    /**< r grid interval (r_max-r_min)/(n_r-1) */
    real z_min;     /**< minimum z coordinate in the grid */
    real z_max;     /**< maximum z coordinate in the grid */
    real z_grid;    /**< z grid interval (z_max-z_min)/(n_r-1) */
    real psi0;      /**< sqrt(psi) value at magnetic axis */
    real psi1;      /**< sqrt(psi) value at separatrix */
    real axis_r;    /**< r coordinate of magnetic axis */
    real axis_z;    /**< z coordinate of magnetic axis */
    real* psi;      /**< pointer to start of psi data in offload_array */
    real* B_r;      /**< pointer to start of B_r data in offload_array */
    real* B_phi;    /**< pointer to start of B_phi data in offload_array */
    real* B_z;      /**< pointer to start of B_z data in offload_array */
} B_2D_data;

void B_2D_init_offload(B_2D_offload_data* offload_data, real** offload_array);
void B_2D_free_offload(B_2D_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_2D_init(B_2D_data* Bdata, B_2D_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
void B_2D_eval_B(real B[], real r, real phi, real z, B_2D_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_2D_eval_psi(real psi[], real r, real phi, real z, B_2D_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_2D_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_2D_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_2D_eval_rho(real rho[], real psi, B_2D_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_2D_eval_rho_drho(real rho[], real r, real phi, real z, B_2D_data* Bdata);
#pragma omp declare simd uniform(B)
real B_2D_bicubic(real t_r, real t_z, int i_r, int i_z, int n_r, real* B);
#pragma omp declare simd uniform(Bdata)
void B_2D_eval_B_dB(real B_dB[], real r, real phi, real z, B_2D_data* Bdata);
#pragma omp declare simd uniform(B)
void B_2D_bicubic_derivs(real B_dB_component[], real t_r, real t_z, int i_r,
                         int i_z, int n_r, real r_grid, real z_grid, real* B);
real B_2D_get_axis_r(B_2D_data* Bdata);
real B_2D_get_axis_z(B_2D_data* Bdata);
#pragma omp end declare target   
#endif
