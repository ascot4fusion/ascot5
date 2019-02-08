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
#include "B_field.h"
#include "math.h"

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
 * @param phase
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_eval(real mhd_dmhd[10], real phase, real r, real phi, real z, real t,
               boozer_data* boozerdata, mhd_data* mhddata,
               B_field_data* Bdata) {

    a5err err = 0;
    real psi_dpsi[4];
    real thetazeta[12];
    B_field_eval_psi_dpsi(psi_dpsi, r, phi, z, Bdata);
    boozer_eval_thetazeta(thetazeta, r, phi, z, psi_dpsi,boozerdata);

    for(int i = 0; i <  mhddata->n_modes; i++){
        /*get interpolated values */
        real a_da[6];
        interp2Dcomp_eval_df(a_da, &(mhddata->alpha_nm[i]),r,t);

        real phi_dphi[6];
        interp2Dcomp_eval_df(phi_dphi,&(mhddata->phi_nm[i]),r,t);

        /* These are used frequently, so store them in separate variables */
        real mhdarg = mhddata->nmode[i] * thetazeta[4]
                    - mhddata->mmode[i] * thetazeta[0]
                    - mhddata->omega_nm[i] * t + phase;
        real sinmhd = sin(mhd);
        real cosmhd = cos(mhd);

        /*sum over modes to get alpha, phi */
        /*possible normalization errors*/

        mhd_dmhd[0] +=     a_da[0] * mhddata->amplitude_nm[i] * sinmhd;
        mhd_dmhd[5] += phi_dphi[0] * mhddata->amplitude_nm[i] * sinmhd;

        /* Time derivatives */
        mhd_dmhd[1] +=     - a_da[0] * mhddata->amplitude_nm[i]
                                     * mhddata->omega_nm[i] * cosmhd
                           + a_da[2] * mhddata->amplitude_nm[i] * sinmhd;

        mhd_dmhd[6] += - phi_dphi[0] * mhddata->amplitude_nm[i]
                                     * mhddata->omega_nm[i] * cosmhd
                       + phi_dphi[2] * mhddata->amplitude_nm[i] * sinmhd;

        /*following code could be written better*/
        /*r component of gradients */

        mhd_dmhd[2] += mhddata->amplitude_nm[i]
            * (  a_da[1] * psi_dpsi[1] * sinmhd
               - a_da[0] * mhddata->mmode[i] * thetazeta[1] * cosmhd
               + a_da[0] * mhddata->nmode[i] * thetazeta[5] * cosmhd);


        mhd_dmhd[7] += mhddata->amplitude_nm[i]
            * (   phi_dphi[1] * psi_dpsi[1] * sinmhd
                - phi_dphi[0] * mhddata->mmode[i] * thetazeta[1] * cosmhd
                + phi_dphi[0] * mhddata->nmode[i] * thetazeta[5] * cosmhd);

        /*phi component of gradients */

        mhd_dmhd[3] += (1/r) * mhddata->amplitude_nm[i]
            * (  a_da[1] * psi_dpsi[2] * sinmhd
               - a_da[0] * mhddata->mmode[i] * thetazeta[2] * cosmhd
               + a_da[0] * mhddata->nmode[i] * thetazeta[6] * cosmhd);

        mhd_dmhd[8] += (1/r) * mhddata->amplitude_nm[i]
            * (   phi_dphi[1] * psi_dpsi[2] * sinmhd
                - phi_dphi[0] * mhddata->mmode[i] * thetazeta[2] * cosmhd
                + phi_dphi[0] * mhddata->nmode[i] * thetazeta[6] * cosmhd);

        /*z component of gradients */

        mhd_dmhd[4] += mhddata->amplitude_nm[i]
            * (   a_da[1] * psi_dpsi[3] * sinmhd
                - a_da[0] * mhddata->mmode[i] * thetazeta[3] * cosmhd
                + a_da[0] * mhddata->nmode[i] * thetazeta[7] * cosmhd);

        mhd_dmhd[9] += mhddata->amplitude_nm[i]
            * (   phi_dphi[1] * psi_dpsi[3] * sinmhd
                - phi_dphi[0] * mhddata->mmode[i] * thetazeta[3] * cosmhd
                + phi_dphi[0] * mhddata->nmode[i] * thetazeta[7] * cosmhd);
     }
    return err;
}

/**
 * @brief Evaluate mhd perturbed fields Btilde, Etilde for full orbit
 *
 * The values are stored in the given array as
 * - pert_field[0] = BtildeR
 * - pert_field[1] = BtildePhi
 * - pert_field[2] = BtildeZ
 * - pert_field[3] = EtildeR
 * - pert_field[4] = EtildePhi
 * - pert_field[5] = EtildeZ
 *
 * @param pert_field perturbation field components
 * @param phase
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_perturbations(real pert_field[6], real phase, real r, real phi,
                        real z, real t, boozer_data* boozerdata,
                        mhd_data* mhddata, B_field_data* Bdata){
    a5err err = 0;
    real mhd_dmhd[10];
    mhd_eval(mhd_dmhd, phase, r, phi, z, t, boozerdata, mhddata,Bdata);

    /*  see example of curl evaluation in step_gc_rk4.c, ydot_gc*/
    real B_dB[12];
    B_field_eval_B_dB(B_dB, r, phi, z, Bdata);

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

    pert_field[3] = -1*mhd_dmhd[7] + -1*B[0]*mhd_dmhd[1];
    pert_field[4] = -1*mhd_dmhd[8] + -1*B[1]*mhd_dmhd[1];
    pert_field[5] = -1*mhd_dmhd[9] + -1*B[2]*mhd_dmhd[1];

    return err;
}
