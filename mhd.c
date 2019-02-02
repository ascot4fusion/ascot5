/**
 * @file mhd.c
 * @brief Module for evaluating MHD parameters.
 */
#include <stdlib.h>
#include "ascot5.h"
#include "print.h"
#include "error.h"
#include "mhd.h"
#include "boozer.h"
#include "spline/interp.h"

/**
 * @brief Load MHD data and prepare parameters for offload.
 *
 * This function fills the MHD offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * It is assumed that the offload_data struct is completely filled before
 * calling this function (except for the offload_array_length and psigrid).
 * Furthermore, offload array should contain following data:
 *
 * - offload_array[j*npsi + i] : alpha(mode_j, psi_i).
 * - offload_array[n_modes*npsi + j*npsi + i] : phi(mode_j, psi_i).
 *
 * 1D splines are constructed here and stored to offload array which is
 * reallocated.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 */
int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array) {

    /* Allocate array for storing 2D spline coefficients for two quantities
     * for each mode                                                         */
    real* coeff_array = (real*)malloc(2 * NSIZE_COMP2D * offload_data->n_modes
                                      * offload_data->npsi * offload_data->ntime
                                      * sizeof(real));

    /* Go through all modes, and evaluate and store coefficients for each */
    int err      = 0;
    int datasize = offload_data->npsi * offload_data->ntime;
    int n_modes  = offload_data->n_modes;
    for(int j=0; j<offload_data->n_modes; j++) {

        /* alpha_nm */
        err += interp2Dcomp_init_coeff(
            &coeff_array[NSIZE_COMP2D * datasize * j],
            &(*offload_array)[j*datasize],
            offload_data->npsi, offload_data->ntime,
            NATURALBC, NATURALBC,
            offload_data->psi_min,
            offload_data->psi_max,
            offload_data->t_min,
            offload_data->t_max);

        /* omega_nm */
        err += interp2Dcomp_init_coeff(
            &coeff_array[NSIZE_COMP2D * datasize * (n_modes + j)],
            &(*offload_array)[j*datasize],
            offload_data->npsi,
            offload_data->ntime,
            NATURALBC, NATURALBC,
            offload_data->psi_min,
            offload_data->psi_max,
            offload_data->t_min,
            offload_data->t_max);
    }

    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = 2 * NSIZE_COMP2D
        * offload_data->n_modes * offload_data->npsi * offload_data->ntime;

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nMHD input\n");
    print_out(VERBOSE_IO, "Grid: npsi = %4.d psimin = %3.3f psimax = %3.3f\n",
              offload_data->npsi,
              offload_data->psi_min, offload_data->psi_max);
    print_out(VERBOSE_IO, "      ntime = %4.d tmin = %3.3f tmax = %3.3f\n",
              offload_data->ntime,
              offload_data->t_min, offload_data->t_max);

    print_out(VERBOSE_IO, "\nModes:\n");
    for(int j=0; j<n_modes; j++) {
        print_out(VERBOSE_IO,
                  "(n,m) = (%2.d,%2.d) Amplitude = %3.3g Frequency = %3.3g\n",
                  offload_data->nmode[j], offload_data->mmode[j],
                  offload_data->amplitude_nm[j], offload_data->omega_nm[j]);
    }

    return err;
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void mhd_free_offload(mhd_offload_data* offload_data,
                          real** offload_array) {
    free(*offload_array);
}

/**
 * @brief Initialize MHD data struct on target
 *
 * @param mhddata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void mhd_init(mhd_data* mhddata, mhd_offload_data* offload_data,
              real* offload_array) {

    mhddata->n_modes = offload_data->n_modes;

    int n_modes  = offload_data->n_modes;
    int datasize = NSIZE_COMP2D * offload_data->npsi * offload_data->ntime;

    for(int j=0; j<mhddata->n_modes; j++) {
        mhddata->nmode[j]        = offload_data->nmode[j];
        mhddata->mmode[j]        = offload_data->mmode[j];
        mhddata->amplitude_nm[j] = offload_data->amplitude_nm[j];
        mhddata->omega_nm[j]     = offload_data->omega_nm[j];

        interp2Dcomp_init_spline(&(mhddata->alpha_nm[j]),
                                 &(offload_array[j*datasize]),
                                 offload_data->npsi, offload_data->ntime,
                                 NATURALBC, NATURALBC,
                                 offload_data->psi_min, offload_data->psi_max,
                                 offload_data->t_min,   offload_data->t_max);

        interp2Dcomp_init_spline(&(mhddata->phi_nm[j]),
                                 &(offload_array[(n_modes + j)*datasize]),
                                 offload_data->npsi, offload_data->ntime,
                                 NATURALBC, NATURALBC,
                                 offload_data->psi_min, offload_data->psi_max,
                                 offload_data->t_min,   offload_data->t_max);
    }
}

/**
 * @brief Evaluate the needed quantities from MHD mode for orbit following, i.e. alpha, phi, grad alpha, grad phi, partial t alpha, partial t phi
 * @todo This is just a dummy.
 * @return Non-zero a5err value if evaluation failed, zero otherwise 
 * The values are stored in the given array as:
 * - mhd_dmhd[0] = alpha 
 * - mhd_dmhd[1] = dalpha/dt
 * - mhd_dmhd[2] = grad alpha, r component
 * - mhd_dmhd[3] = grad alpha, phi component
 * - mhd_dmhd[4] = grad alpha, z component
 * - mhd_dmhd[5] = phi
 * - mhd_dmhd[6] = dphi/dt
 * - mhd_dmhd[7] = grad phi, r component
 * - mhd_dmhd[8] = grad phi, phi component
 * - mhd_dmhd[9] = grad phi, z component
 */
a5err mhd_eval(real mhd_dmhd[10], real phase, real r, real phi, real z, real t, boozer_data* boozerdata, mhd_data* MHDdata) {
    a5err err = 0;
    real ptz_dptz[12];
    boozer_eval_gradients(ptz_dptz, r, phi, z, boozerdata);
    int i; 
    for(i = 0; i <  MHDdata->n_modes; i++){
	/*get interpolated values */
        real a_da[6];
        interp2Dcomp_eval_df(a_da, &(MHDdata->alpha_nm[i]),r,t);

        real phi_dphi[6];
	interp2Dcomp_eval_df(phi_dphi,&(MHDdata->alpha_nm[i]),r,t);
        
        /*this is used frequently, so define here*/
        real mhdarg = (MHDdata->nmode[i])* ptz_dptz[8]-(MHDdata->mmode[i])*ptz_dptz[4]-(MHDdata->omega_nm[i])*t +phase;

        /*sum over modes to get alpha, phi */
	/*possible normalization errors*/
        
        mhd_dmhd[0] += a_da[0]*(MHDdata->amplitude_nm[i])*sin(mhdarg);
        mhd_dmhd[5] += phi_dphi[0]*(MHDdata->amplitude_nm[i])*sin(mhdarg);

	/*time derivs */

	mhd_dmhd[1] += -1*a_da[0]*(MHDdata->amplitude_nm[i])*(MHDdata->omega_nm[i])*cos(mhdarg);
	mhd_dmhd[6] += -1*phi_dphi[0]*(MHDdata->amplitude_nm[i])*(MHDdata->omega_nm[i])*cos(mhdarg);
	/*following code could be written better*/
	/*r component of gradients */

	mhd_dmhd[2] += (MHDdata->amplitude_nm[i])*(a_da[1]*ptz_dptz[1]*sin(mhdarg) + -1*a_da[0]*(MHDdata->mmode[i])*ptz_dptz[5]*cos(mhdarg) + a_da[0]*(MHDdata->nmode[i])*ptz_dptz[9]*cos(mhdarg));


	mhd_dmhd[7] += (MHDdata->amplitude_nm[i])*(phi_dphi[1]*ptz_dptz[1]*sin(mhdarg) + -1*phi_dphi[0]*(MHDdata->mmode[i])*ptz_dptz[5]*cos(mhdarg) + phi_dphi[0]*(MHDdata->nmode[i])*ptz_dptz[9]*cos(mhdarg));

	/*phi component of gradients */

	mhd_dmhd[3] += (MHDdata->amplitude_nm[i])*(a_da[1]*ptz_dptz[2]*sin(mhdarg) + -1*a_da[0]*(MHDdata->mmode[i])*ptz_dptz[6]*cos(mhdarg) + a_da[0]*(MHDdata->nmode[i])*ptz_dptz[10]*cos(mhdarg));

	mhd_dmhd[8] += (MHDdata->amplitude_nm[i])*(phi_dphi[1]*ptz_dptz[2]*sin(mhdarg) + -1*phi_dphi[0]*(MHDdata->mmode[i])*ptz_dptz[6]*cos(mhdarg) + phi_dphi[0]*(MHDdata->nmode[i])*ptz_dptz[10]*cos(mhdarg));

	/*z component of gradients */

	mhd_dmhd[4] += (MHDdata->amplitude_nm[i])*(a_da[1]*ptz_dptz[3]*sin(mhdarg) + -1*a_da[0]*(MHDdata->mmode[i])*ptz_dptz[7]*cos(mhdarg) + a_da[0]*(MHDdata->nmode[i])*ptz_dptz[11]*cos(mhdarg)); 
 
	mhd_dmhd[9] += (MHDdata->amplitude_nm[i])*(phi_dphi[1]*ptz_dptz[3]*sin(mhdarg) + -1*phi_dphi[0]*(MHDdata->mmode[i])*ptz_dptz[7]*cos(mhdarg) + phi_dphi[0]*(MHDdata->nmode[i])*ptz_dptz[11]*cos(mhdarg));
     }
    return err;
}
