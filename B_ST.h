/**
 * @file B_ST.h
 * @brief Header file for B_ST.c
 */
#ifndef B_ST_H
#define B_ST_H
#include "ascot5.h"

/**
 * @brief Stellarator magnetic field parameters that will be offloaded to target
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
    real phi_grid;              /**< phi grid interval 2pio/(n_phi-1) */
    real axis_r;                /**< r coordinate of magnetic axis */
    real axis_z;                /**< z coordinate of magnetic axis */
    int offload_array_length;   /**< number of elements in offload_array */
} B_ST_offload_data;

/**
 * @brief Stellarator magnetic field parameters on the target
 */
typedef struct {
    int n_r;        /**< number of r grid points */
    int n_z;        /**< number of z grid points */
    int n_phi;      /**< number of phi grid points */
    int periods;    /**< number of toroidal periods */
    real r_min;     /**< minimum r coordinate in the grid */
    real r_max;     /**< maximum r coordinate in the grid */
    real r_grid;    /**< r grid interval (r_max-r_min)/(n_r-1) */
    real z_min;     /**< minimum z coordinate in the grid */
    real z_max;     /**< maximum z coordinate in the grid */
    real z_grid;    /**< z grid interval (z_max-z_min)/(n_r-1) */
    real phi_min;   /**< minimum phi coordinate in the grid */
    real phi_max;   /**< maximum phi coordinate in the grid */
    real phi_grid;  /**< phi grid interval 2pi/(n_phi-1) */
    real axis_r;    /**< r coordinate of magnetic axis */
    real axis_z;    /**< z coordinate of magnetic axis */
    real* s;        /**< pointer to start of psi data in offload_array */
    real* B_r;      /**< pointer to start of B_r data in offload_array */
    real* B_phi;    /**< pointer to start of B_phi data in offload_array */
    real* B_z;      /**< pointer to start of B_z data in offload_array */
} B_ST_data;

void B_ST_init_offload(B_ST_offload_data* offload_data, real** offload_array);
void B_ST_free_offload(B_ST_offload_data* offload_data, real** offload_array);

#pragma omp declare target
void B_ST_init(B_ST_data* Bdata, B_ST_offload_data* offload_data,
               real* offload_array);
#pragma omp declare simd uniform(Bdata)
void B_ST_eval_B(real B[], real r, real phi, real z, B_ST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_ST_eval_psi(real psi[], real r, real phi, real z, B_ST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_ST_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_ST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_ST_eval_rho(real rho[], real psi, B_ST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_ST_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_ST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
void B_ST_eval_B_dB(real B_dB[], real r, real phi, real z, B_ST_data* Bdata);
#pragma omp declare simd uniform(Bdata)
real B_ST_get_axis_r(B_ST_data* Bdata);
real B_ST_get_axis_z(B_ST_data* Bdata);
#pragma omp end declare target   
#endif
