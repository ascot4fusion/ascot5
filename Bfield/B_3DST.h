/**
 * @file B_3DS_T.h
 * @brief Header file for B_3DS_T.c
 */
#ifndef B_3DS_H
#define B_3DS_H

#define N_MAX_TIME_SLICE 10

#include "../ascot5.h"
#include "../spline/interp2D.h" /* for 2D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */
#include "B_3DS.h"
/**
 * @brief 3D magnetic field parameters on the host
 */
typedef struct {
    int n_time;
    real time[N_MAX_TIME_SLICE];

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
    
    real psi0[N_MAX_TIME_SLICE];                  /**< sqrt(psi) value at magnetic axis */
    real psi1[N_MAX_TIME_SLICE];                  /**< sqrt(psi) value at separatrix */
    real axis_r[N_MAX_TIME_SLICE];                /**< r coordinate of magnetic axis */
    real axis_z[N_MAX_TIME_SLICE];                /**< z coordinate of magnetic axis */
    int offload_array_length;   /**< number of elements in offload_array */
} B_3DS_T_offload_data;

/**
 * @brief 3D magnetic field parameters on the target
 */
typedef struct {
    int n_time;                 /**< number of time slices in B & psi data */
    real time[N_MAX_TIME_SLICE];/**< values of time for each time slice */
    B_3DS_data Bslice[N_MAX_TIME_SLICE];/**< B3DS grids for each time slice */
} B_3DS_T_data;

void B_3DS_T_init_offload(B_3DS_T_offload_data* offload_data, real** offload_array);
void B_3DS_T_free_offload(B_3DS_T_offload_data* offload_data, real** offload_array);

#pragma omp declare target
int B_3DS_T_init(B_3DS_T_data* Bdata, B_3DS_T_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_T_eval_psi(real psi[], real r, real phi, real z, real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_T_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_T_eval_rho(real rho[], real psi, real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_T_eval_rho_drho(real rho_drho[], real r, real phi, real z, real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_T_eval_B(real B[], real r, real phi, real z, real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
a5err B_3DS_T_eval_B_dB(real B_dB[], real r, real phi, real z, real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DS_T_get_axis_r(real time, B_3DS_T_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_3DS_T_get_axis_z(real time, B_3DS_T_data* Bdata);
#pragma omp end declare target   
#endif
