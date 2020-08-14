/**
 * @file mhd_stat.c
 * @brief MHD module for stationary amplitudes (eigenmodes).
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "../print.h"
#include "../error.h"
#include "../boozer.h"
#include "../spline/interp.h"
#include "../B_field.h"
#include "../math.h"
#include "mhd_stat.h"

/**
 * @brief Load MHD data and prepare parameters for offload.
 *
 * This function fills the MHD offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * It is assumed that the offload_data struct is completely filled before
 * calling this function (except for the offload_array_length and rhogrid).
 * Furthermore, offload array should contain following data:
 *
 * - offload_array[j*nrho + i] : alpha(mode_j, rho_i).
 * - offload_array[n_modes*nrho + j*nrho + i] : phi(mode_j, rho_i).
 *
 * 1D splines are constructed here and stored to offload array which is
 * reallocated.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 */
int mhd_stat_init_offload(mhd_stat_offload_data* offload_data,
                          real** offload_array) {

    /* Allocate array for storing 1D spline coefficients for two quantities
     * for each mode                                                         */
    real* coeff_array = (real*)malloc(2 * NSIZE_COMP1D * offload_data->n_modes
                                      * offload_data->nrho * sizeof(real));

    /* Go through all modes, and evaluate and store coefficients for each */
    int err      = 0;
    int datasize = offload_data->nrho;
    int n_modes  = offload_data->n_modes;
    for(int j=0; j<offload_data->n_modes; j++) {

        /* alpha_nm */
        err += interp1Dcomp_init_coeff(
            &coeff_array[NSIZE_COMP1D * datasize * j],
            &(*offload_array)[j*datasize],
            offload_data->nrho,
            NATURALBC,
            offload_data->rho_min,
            offload_data->rho_max);

        /* phi_nm */
        err += interp1Dcomp_init_coeff(
            &coeff_array[NSIZE_COMP1D * datasize * (n_modes + j)],
            &(*offload_array)[(n_modes + j)*datasize],
            offload_data->nrho,
            NATURALBC,
            offload_data->rho_min,
            offload_data->rho_max);
    }

    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = 2 * NSIZE_COMP1D
        * offload_data->n_modes * offload_data->nrho;

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nMHD (stationary) input\n");
    print_out(VERBOSE_IO, "Grid: nrho = %4.d rhomin = %3.3f rhomax = %3.3f\n",
              offload_data->nrho,
              offload_data->rho_min, offload_data->rho_max);

    print_out(VERBOSE_IO, "\nModes:\n");
    for(int j=0; j<n_modes; j++) {
        print_out(VERBOSE_IO,
                  "(n,m) = (%2.d,%2.d) Amplitude = %3.3g Frequency = %3.3g"
                  " Phase = %3.3g\n",
                  offload_data->nmode[j], offload_data->mmode[j],
                  offload_data->amplitude_nm[j], offload_data->omega_nm[j],
                  offload_data->phase_nm[j]);
    }

    return err;
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void mhd_stat_free_offload(mhd_stat_offload_data* offload_data,
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
void mhd_stat_init(mhd_stat_data* mhddata, mhd_stat_offload_data* offload_data,
                   real* offload_array) {

    mhddata->n_modes = offload_data->n_modes;
    mhddata->rho_min = offload_data->rho_min;
    mhddata->rho_max = offload_data->rho_max;

    int n_modes  = offload_data->n_modes;
    int datasize = NSIZE_COMP1D * offload_data->nrho;

    for(int j=0; j<mhddata->n_modes; j++) {
        mhddata->nmode[j]        = offload_data->nmode[j];
        mhddata->mmode[j]        = offload_data->mmode[j];
        mhddata->amplitude_nm[j] = offload_data->amplitude_nm[j];
        mhddata->omega_nm[j]     = offload_data->omega_nm[j];
        mhddata->phase_nm[j]     = offload_data->phase_nm[j];

        interp1Dcomp_init_spline(&(mhddata->alpha_nm[j]),
                                 &(offload_array[j*datasize]),
                                 offload_data->nrho,
                                 NATURALBC,
                                 offload_data->rho_min, offload_data->rho_max);

        interp1Dcomp_init_spline(&(mhddata->phi_nm[j]),
                                 &(offload_array[(n_modes + j)*datasize]),
                                 offload_data->nrho,
                                 NATURALBC,
                                 offload_data->rho_min, offload_data->rho_max);

    }
}

/**
 * @brief Evaluate the needed quantities from MHD mode for orbit following
 *
 * The quantities to be evaluated are alpha, phi, grad alpha, grad phi,
 * partial t alpha, partial t phi
 *
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
 *
 * @param mhd_dmhd
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_stat_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
                    boozer_data* boozerdata, mhd_stat_data* mhddata) {

    a5err err = 0;

    real ptz[12];
    int isinside;
    if(!err) {
        err = boozer_eval_psithetazeta(ptz, &isinside, r, phi, z, boozerdata);
    }
    real rho[2];
    if(!err) {
        err = boozer_eval_rho_drho(rho, ptz[0], boozerdata);
    }

    int iterations = mhddata->n_modes;

    /* Initialize values */
    for(int i=0; i<10; i++) {
        mhd_dmhd[i] = 0;
    }

    /* Check that we are within MHD grid */
    if(rho[0] < mhddata->rho_min || rho[0] > mhddata->rho_max) {
        isinside = 0;
    }

    /* Skip evaluation if evaluation failed or point outside the grid. */
    if(err || !isinside) {
        iterations = 0;
    }

    int interperr = 0;
    for(int i = 0; i < iterations; i++){
        /* Get interpolated values */
        real a_da[3], phi_dphi[3];
        interperr += interp1Dcomp_eval_df(a_da, &(mhddata->alpha_nm[i]),
                                          rho[0]);
        interperr += interp1Dcomp_eval_df(phi_dphi, &(mhddata->phi_nm[i]),
                                          rho[0]);

        /* The interpolation returns dx/drho but we require dx/dpsi.
         * The second order derivatives are not needed anywhere */
        a_da[1]     *= rho[1];
        phi_dphi[1] *= rho[1];

        /* These are used frequently, so store them in separate variables */
        real mhdarg = mhddata->nmode[i] * ptz[8]
                    - mhddata->mmode[i] * ptz[4]
                    - mhddata->omega_nm[i] * t
                    + mhddata->phase_nm[i];
        real sinmhd = sin(mhdarg);
        real cosmhd = cos(mhdarg);

        /* Sum over modes to get alpha, phi */
        mhd_dmhd[0] +=     a_da[0] * mhddata->amplitude_nm[i] * cosmhd;
        mhd_dmhd[5] += phi_dphi[0] * mhddata->amplitude_nm[i] * cosmhd;

        /* Time derivatives */
        mhd_dmhd[1] +=     a_da[0] * mhddata->amplitude_nm[i]
                                   * mhddata->omega_nm[i] * sinmhd;
        mhd_dmhd[6] += phi_dphi[0] * mhddata->amplitude_nm[i]
                                   * mhddata->omega_nm[i] * sinmhd;

        /* R component of gradients */
        mhd_dmhd[2] += mhddata->amplitude_nm[i]
            * (  a_da[1] * ptz[1] * cosmhd
               + a_da[0] * mhddata->mmode[i] * ptz[5] * sinmhd
               - a_da[0] * mhddata->nmode[i] * ptz[9] * sinmhd);
        mhd_dmhd[7] += mhddata->amplitude_nm[i]
            * (   phi_dphi[1] * ptz[1] * cosmhd
                + phi_dphi[0] * mhddata->mmode[i] * ptz[5] * sinmhd
                - phi_dphi[0] * mhddata->nmode[i] * ptz[9] * sinmhd);

        /* phi component of gradients */
        mhd_dmhd[3] += (1/r) * mhddata->amplitude_nm[i]
            * (  a_da[1] * ptz[2] * cosmhd
               + a_da[0] * mhddata->mmode[i] * ptz[6]  * sinmhd
               - a_da[0] * mhddata->nmode[i] * ptz[10] * sinmhd);
        mhd_dmhd[8] += (1/r) * mhddata->amplitude_nm[i]
            * (   phi_dphi[1] * ptz[2] * cosmhd
                + phi_dphi[0] * mhddata->mmode[i] * ptz[6]  * sinmhd
                - phi_dphi[0] * mhddata->nmode[i] * ptz[10] * sinmhd);

        /* z component of gradients */
        mhd_dmhd[4] += mhddata->amplitude_nm[i]
            * (   a_da[1] * ptz[3] * cosmhd
                + a_da[0] * mhddata->mmode[i] * ptz[7]  * sinmhd
                - a_da[0] * mhddata->nmode[i] * ptz[11] * sinmhd);
        mhd_dmhd[9] += mhddata->amplitude_nm[i]
            * (   phi_dphi[1] * ptz[3] * cosmhd
                + phi_dphi[0] * mhddata->mmode[i] * ptz[7]  * sinmhd
                - phi_dphi[0] * mhddata->nmode[i] * ptz[11] * sinmhd);
    }

    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_MHD );
    }
    return err;
}

/**
 * @brief Evaluate perturbed fields Btilde, Etilde and potential Phi explicitly
 *
 * The values are stored in the given array as
 * - pert_field[0] = BtildeR
 * - pert_field[1] = BtildePhi
 * - pert_field[2] = BtildeZ
 * - pert_field[3] = EtildeR
 * - pert_field[4] = EtildePhi
 * - pert_field[5] = EtildeZ
 * - pert_field[6] = Phi
 *
 * Only the perturbation values for the magnetic field are returned if
 * pertonly=1, otherwise, the total perturbed field is returned. This is done to
 * avoid double evaluation of the magnetic field e.g. in field line tracing.
 * For electric field only the perturbation component is returned always.
 *
 * @param pert_field perturbation field components
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param pertonly flag whether to return the whole field or only perturbation
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_stat_perturbations(real pert_field[7], real r, real phi, real z,
                             real t, int pertonly, boozer_data* boozerdata,
                             mhd_stat_data* mhddata, B_field_data* Bdata) {
    a5err err = 0;
    real mhd_dmhd[10];
    if(!err) {
        err = mhd_stat_eval(mhd_dmhd, r, phi, z, t, boozerdata, mhddata);
    }
    /*  see example of curl evaluation in step_gc_rk4.c, ydot_gc*/
    real B_dB[15];
    if(!err) {
        err = B_field_eval_B_dB(B_dB, r, phi, z, t, Bdata);
    }

    if(!err) {
        real B[3];
        B[0] = B_dB[0];
        B[1] = B_dB[4];
        B[2] = B_dB[8];

        real curlB[3];
        curlB[0] = B_dB[10]/r - B_dB[7];
        curlB[1] = B_dB[3] - B_dB[9];
        curlB[2] = (B[1] - B_dB[2])/r + B_dB[5];

        real gradalpha[3];
        gradalpha[0] = mhd_dmhd[2];
        gradalpha[1] = mhd_dmhd[3];
        gradalpha[2] = mhd_dmhd[4];

        real gradalphacrossB[3];

        math_cross(gradalpha, B, gradalphacrossB);

        pert_field[0] = mhd_dmhd[0]*curlB[0] + gradalphacrossB[0];
        pert_field[1] = mhd_dmhd[0]*curlB[1] + gradalphacrossB[1];
        pert_field[2] = mhd_dmhd[0]*curlB[2] + gradalphacrossB[2];

        pert_field[3] = -mhd_dmhd[7] - B[0]*mhd_dmhd[1];
        pert_field[4] = -mhd_dmhd[8] - B[1]*mhd_dmhd[1];
        pert_field[5] = -mhd_dmhd[9] - B[2]*mhd_dmhd[1];
        pert_field[6] = mhd_dmhd[5];

        if(!pertonly) {
            pert_field[0] += B[0];
            pert_field[1] += B[1];
            pert_field[2] += B[2];
        }
    }

    return err;
}
