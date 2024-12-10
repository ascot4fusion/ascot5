/**
 * @file mhd_nonstat.h
 * @brief Header file for mhd_nonstat.c
 */
#ifndef MHD_NONSTAT_H
#define MHD_NONSTAT_H

#include "../ascot5.h"
#include "../offload.h"
#include "../error.h"
#include "../boozer.h"
#include "../spline/interp.h"
#include "../B_field.h"

/**
 * @brief MHD parameters
 */
typedef struct {
    int n_modes;        /**< Number of modes            */
    real rho_min;       /**< rho grid minimum value     */
    real rho_max;       /**< rho grid maximum value     */
    int* nmode;         /**< Toroidal mode numbers      */
    int* mmode;         /**< Poloidal mode numbers      */
    real* amplitude_nm; /**< Amplitude of each mode     */
    real* omega_nm;     /**< Toroidal frequency [rad/s] */
    real* phase_nm;     /**< Phase of each mode [rad]   */

    /**
     * @brief 2D splines (rho,time) for each mode's magnetic eigenfunction
     */
    interp2D_data* alpha_nm;
    /**
     * @brief 2D splines (rho,time) for each mode's electric eigenfunction
     */
    interp2D_data* phi_nm;
} mhd_nonstat_data;

int mhd_nonstat_init(mhd_nonstat_data* data, int nmode, int nrho, int ntime,
                     real rhomin, real rhomax, real tmin, real tmax,
                     int* moden, int* modem, real* amplitude_nm,
                     real* omega_nm, real* phase_nm, real* alpha, real* phi);
void mhd_nonstat_free(mhd_nonstat_data* data);
void mhd_nonstat_offload(mhd_nonstat_data* data);
GPU_DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, includemode)
a5err mhd_nonstat_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
                       int includemode, boozer_data* boozerdata,
                       mhd_nonstat_data* mhddata, B_field_data* Bdata);
DECLARE_TARGET_END
GPU_DECLARE_TARGET_SIMD_UNIFORM(boozerdata, mhddata, Bdata, pertonly,	\
                                 includemode)
a5err mhd_nonstat_perturbations(real pert_field[7], real r, real phi, real z,
                                real t, int pertonly, int includemode,
                                boozer_data* boozerdata,
                                mhd_nonstat_data* mhddata, B_field_data* Bdata);
DECLARE_TARGET_END

#endif
