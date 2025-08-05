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
#include "../mhd.h"
#include "mhd_stat.h"

/**
 * @brief Load MHD data
 *
 * @param offload_data pointer to the data struct
 *
 * @return zero if initialization succeeded.
 */
int mhd_stat_init(mhd_stat_data* data, int nmode, int nrho,
                  real rhomin, real rhomax, int* moden, int* modem,
                  real* amplitude_nm, real* omega_nm, real* phase_nm,
                  real* alpha, real* phi) {

    int err = 0;
    data->n_modes = nmode;
    data->rho_min = rhomin;
    data->rho_max = rhomax;
    data->nmode = (int*) malloc(nmode * sizeof(int));
    data->mmode = (int*) malloc(nmode * sizeof(int));
    data->omega_nm = (real*) malloc(nmode * sizeof(real));
    data->phase_nm = (real*) malloc(nmode * sizeof(real));
    data->amplitude_nm = (real*) malloc(nmode * sizeof(real));
    data->phi_nm = (interp1D_data*) malloc(nmode * sizeof(interp1D_data));
    data->alpha_nm = (interp1D_data*) malloc(nmode * sizeof(interp1D_data));
    for(int i = 0; i < nmode; i++) {
        data->nmode[i] = moden[i];
        data->mmode[i] = modem[i];
        data->omega_nm[i] = omega_nm[i];
        data->phase_nm[i] = phase_nm[i];
        data->amplitude_nm[i] = amplitude_nm[i];

        err = interp1Dcomp_setup(&data->alpha_nm[i], &alpha[i*nrho],
                                 nrho, NATURALBC, rhomin, rhomax);
        if(err) {
            print_err("Error: Failed to initialize splines.\n");
            return err;
        }
        err = interp1Dcomp_setup(&data->phi_nm[i], &phi[i*nrho],
                                 nrho, NATURALBC, rhomin, rhomax);
        if(err) {
            print_err("Error: Failed to initialize splines.\n");
            return err;
        }
    }

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nMHD (stationary) input\n");
    print_out(VERBOSE_IO, "Grid: nrho = %4.d rhomin = %3.3f rhomax = %3.3f\n",
              nrho, data->rho_min, data->rho_max);

    print_out(VERBOSE_IO, "\nModes:\n");
    for(int i = 0; i < nmode; i++) {
        print_out(VERBOSE_IO,
                  "(n,m) = (%2.d,%2.d) Amplitude = %3.3g Frequency = %3.3g"
                  " Phase = %3.3g\n",
                  data->nmode[i], data->mmode[i], data->amplitude_nm[i],
                  data->omega_nm[i], data->phase_nm[i]);
    }

    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void mhd_stat_free(mhd_stat_data* data) {
    for(int i = 0; i < data->n_modes; i++) {
        free(data->phi_nm[i].c);
        free(data->alpha_nm[i].c);
    }
    free(data->phi_nm);
    free(data->alpha_nm);
    free(data->nmode);
    free(data->mmode);
    free(data->phase_nm);
    free(data->omega_nm);
    free(data->amplitude_nm);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void mhd_stat_offload(mhd_stat_data* data) {
    //TODO: Implement
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
 * @param includemode mode number to include or MHD_INCLUDE_ALL
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_stat_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
                    int includemode, boozer_data* boozerdata,
                    mhd_stat_data* mhddata, B_field_data* Bdata) {

    a5err err = 0;

    real ptz[12];
    int isinside;
    if(!err) {
        err = boozer_eval_psithetazeta(ptz, &isinside, r, phi, z, Bdata,
                                       boozerdata);
    }
    real rho[2];
    if(!err && isinside) {
        err = B_field_eval_rho(rho, ptz[0], Bdata);
    }

    /* Initialize values */
    for(int i=0; i<10; i++) {
        mhd_dmhd[i] = 0;
    }

    int interperr = 0;
    for(int i = 0; i < mhddata->n_modes; i++){
        if( includemode != MHD_INCLUDE_ALL && includemode != i ) { continue; }
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

    /* Omit evaluation if point outside the boozer or mhd grid. */
    if(!isinside || interperr) {
        for(int i=0; i<10; i++) {
            mhd_dmhd[i] = 0;
        }
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
 * @param includemode mode number to include or MHD_INCLUDE_ALL
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_stat_perturbations(real pert_field[7], real r, real phi, real z,
                             real t, int pertonly, int includemode,
                             boozer_data* boozerdata, mhd_stat_data* mhddata,
                             B_field_data* Bdata) {
    a5err err = 0;
    real mhd_dmhd[10];
    if(!err) {
        err = mhd_stat_eval(mhd_dmhd, r, phi, z, t, includemode, boozerdata,
                            mhddata, Bdata);
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

/**
 * @brief Evaluate perturbed fields Btilde, Etilde and potential Phi explicitly,
 * and its time derivatives.
 *
 * The values are stored in the given array as
 * - pert_field[0] = BtildeR
 * - pert_field[1] = BtildePhi
 * - pert_field[2] = BtildeZ
 * - pert_field[3] = EtildeR
 * - pert_field[4] = EtildePhi
 * - pert_field[5] = EtildeZ
 * - pert_field[6] = Phi
 * - pert_field[7] = dBtildeR/dt
 * - pert_field[8] = dBtildePhi/dt
 * - pert_field[9] = dBtildeZ/dt
 * - pert_field[10] = dEtildeR/dt
 * - pert_field[11] = dEtildePhi/dt
 * - pert_field[12] = dEtildeZ/dt
 * - pert_field[13] = dPhi/dt
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
 * @param includemode mode number to include or MHD_INCLUDE_ALL
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_stat_eval_perturbations_dt(real pert_field[14], real r, real phi, real z,
                                    real t, int pertonly, int includemode,
                                    boozer_data* boozerdata, mhd_stat_data* mhddata,
                                    B_field_data* Bdata){
    a5err err = 0;
    real ptz[12];
    real gradalpha[3];
    real gradphipot[3];
    real gradalpha_dt[3];
    real gradphipot_dt[3];
    real phipot, phipot_dt;
    real alpha, alpha_dt, alpha_dt2;
    int isinside;
    if(!err) {
        err = boozer_eval_psithetazeta(ptz, &isinside, r, phi, z, Bdata,
                                       boozerdata);
    }
    real rho[2];
    if(!err && isinside) {
        err = B_field_eval_rho(rho, ptz[0], Bdata);
    }
    
    real B_dB[15];
    if(!err) {
        err = B_field_eval_B_dB(B_dB, r, phi, z, t, Bdata);
    }

    /* Get interpolated values */
    real a_da[3], phi_dphi[3];
    int interperr = 0;
    interperr += interp1Dcomp_eval_df(a_da, &(mhddata->alpha_nm[includemode]),
                                        rho[0]);
    interperr += interp1Dcomp_eval_df(phi_dphi, &(mhddata->phi_nm[includemode]),
                                        rho[0]);
    if(interperr) {
        a_da[0] = 0;
        a_da[1] = 0;
        a_da[2] = 0;
        phi_dphi[0] = 0;
        phi_dphi[1] = 0;
        phi_dphi[2] = 0;
    }

    /* The interpolation returns dx/drho but we require dx/dpsi.
        * The second order derivatives are not needed anywhere */
    a_da[1]     *= rho[1];
    phi_dphi[1] *= rho[1];

    /* These are used frequently, so store them in separate variables */
    real mhdarg = mhddata->nmode[includemode] * ptz[8]
        - mhddata->mmode[includemode] * ptz[4]
        - mhddata->omega_nm[includemode] * t
        + mhddata->phase_nm[includemode];
    real sinmhd = sin(mhdarg);
    real cosmhd = cos(mhdarg);

    /* Getting the potentials */
    alpha  =     a_da[0] * mhddata->amplitude_nm[includemode] * cosmhd;
    phipot = phi_dphi[0] * mhddata->amplitude_nm[includemode] * cosmhd;

    /* Time derivatives */
    alpha_dt =     a_da[0] * mhddata->amplitude_nm[includemode]
        * mhddata->omega_nm[includemode] * sinmhd;
    phipot_dt = phi_dphi[0] * mhddata->amplitude_nm[includemode]
        * mhddata->omega_nm[includemode] * sinmhd;

    /* Second time derivative: the sign comes from the -omega*t term in the sin() */
    real omega_nm2 = mhddata->omega_nm[includemode] * mhddata->omega_nm[includemode];
    alpha_dt2 =  -  a_da[0] * mhddata->amplitude_nm[includemode]
        * omega_nm2 * cosmhd;

    /* R component of gradients */
    gradalpha[0] = mhddata->amplitude_nm[includemode]
        * (  a_da[1] * ptz[1] * cosmhd
            + a_da[0] * mhddata->mmode[includemode] * ptz[5] * sinmhd
            - a_da[0] * mhddata->nmode[includemode] * ptz[9] * sinmhd);
    gradphipot[0] = mhddata->amplitude_nm[includemode]
        * (   phi_dphi[1] * ptz[1] * cosmhd
            + phi_dphi[0] * mhddata->mmode[includemode] * ptz[5] * sinmhd
            - phi_dphi[0] * mhddata->nmode[includemode] * ptz[9] * sinmhd);

    /* phi component of gradients */
    gradalpha[1] = (1/r) * mhddata->amplitude_nm[includemode]
        * (  a_da[1] * ptz[2] * cosmhd
            + a_da[0] * mhddata->mmode[includemode] * ptz[6]  * sinmhd
            - a_da[0] * mhddata->nmode[includemode] * ptz[10] * sinmhd);
    gradphipot[1] = (1/r) * mhddata->amplitude_nm[includemode]
        * (   phi_dphi[1] * ptz[2] * cosmhd
            + phi_dphi[0] * mhddata->mmode[includemode] * ptz[6]  * sinmhd
            - phi_dphi[0] * mhddata->nmode[includemode] * ptz[10] * sinmhd);

    /* z component of gradients */
    gradalpha[2] = mhddata->amplitude_nm[includemode]
        * (   a_da[1] * ptz[3] * cosmhd
            + a_da[0] * mhddata->mmode[includemode] * ptz[7]  * sinmhd
            - a_da[0] * mhddata->nmode[includemode] * ptz[11] * sinmhd);
    gradphipot[2] = mhddata->amplitude_nm[includemode]
        * (   phi_dphi[1] * ptz[3] * cosmhd
            + phi_dphi[0] * mhddata->mmode[includemode] * ptz[7]  * sinmhd
            - phi_dphi[0] * mhddata->nmode[includemode] * ptz[11] * sinmhd);

    /* Time derivatives of gradients */
    /* R component of gradients */
    gradalpha_dt[0] = mhddata->amplitude_nm[includemode] * mhddata->omega_nm[includemode]
        * (  a_da[1] * ptz[1] * sinmhd
            - a_da[0] * mhddata->mmode[includemode] * ptz[5] * cosmhd
            + a_da[0] * mhddata->nmode[includemode] * ptz[9] * cosmhd);

    gradphipot_dt[0] = mhddata->amplitude_nm[includemode] * mhddata->omega_nm[includemode]
        * (   phi_dphi[1] * ptz[1] * sinmhd
            - phi_dphi[0] * mhddata->mmode[includemode] * ptz[5] * cosmhd
            + phi_dphi[0] * mhddata->nmode[includemode] * ptz[9] * cosmhd);

    /* phi component of gradients */
    gradalpha_dt[1] = (1/r) * mhddata->amplitude_nm[includemode] * mhddata->omega_nm[includemode]
        * (  a_da[1] * ptz[2] * sinmhd
            - a_da[0] * mhddata->mmode[includemode] * ptz[6]  * cosmhd
            + a_da[0] * mhddata->nmode[includemode] * ptz[10] * cosmhd);

    gradphipot_dt[1] = (1/r) * mhddata->amplitude_nm[includemode] * mhddata->omega_nm[includemode]
        * (   phi_dphi[1] * ptz[2] * sinmhd
            - phi_dphi[0] * mhddata->mmode[includemode] * ptz[6]  * cosmhd
            + phi_dphi[0] * mhddata->nmode[includemode] * ptz[10] * cosmhd);

    /* z component of gradients */
    gradalpha_dt[2] = mhddata->amplitude_nm[includemode] * mhddata->omega_nm[includemode]
        * (   a_da[1] * ptz[3] * sinmhd
            - a_da[0] * mhddata->mmode[includemode] * ptz[7]  * cosmhd
            + a_da[0] * mhddata->nmode[includemode] * ptz[11] * cosmhd);

    gradphipot_dt[2] = mhddata->amplitude_nm[includemode] * mhddata->omega_nm[includemode]
        * (   phi_dphi[1] * ptz[3] * sinmhd
            - phi_dphi[0] * mhddata->mmode[includemode] * ptz[7]  * cosmhd
            + phi_dphi[0] * mhddata->nmode[includemode] * ptz[11] * cosmhd);
    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[4];
    B[2] = B_dB[8];

    real curlB[3];
    curlB[0] = B_dB[10]/r - B_dB[7];
    curlB[1] = B_dB[3] - B_dB[9];
    curlB[2] = (B[1] - B_dB[2])/r + B_dB[5];

    real gradalphacrossB[3];
    math_cross(gradalpha, B, gradalphacrossB);

    real gradalphadt_x_B[3];
    math_cross(gradalpha_dt, B, gradalphadt_x_B);

    /** Fields */
    pert_field[0] = alpha*curlB[0] + gradalphacrossB[0];
    pert_field[1] = alpha*curlB[1] + gradalphacrossB[1];
    pert_field[2] = alpha*curlB[2] + gradalphacrossB[2];

    pert_field[3] = -gradphipot[0] - B[0]*alpha_dt;
    pert_field[4] = -gradphipot[1] - B[1]*alpha_dt;
    pert_field[5] = -gradphipot[2] - B[2]*alpha_dt;
    pert_field[6] = phipot;

    /** Instantaneous field evolution */
    pert_field[7] = alpha_dt * curlB[0] + gradalphadt_x_B[0];
    pert_field[8] = alpha_dt * curlB[1] + gradalphadt_x_B[1];
    pert_field[9] = alpha_dt * curlB[2] + gradalphadt_x_B[2];
    pert_field[10] = - gradphipot_dt[0] - B[0]*alpha_dt2;
    pert_field[11] = - gradphipot_dt[1] - B[1]*alpha_dt2;
    pert_field[12] = - gradphipot_dt[2] - B[2]*alpha_dt2;
    pert_field[13] = phipot_dt;

    if(!pertonly) {
        pert_field[0] += B[0];
        pert_field[1] += B[1];
        pert_field[2] += B[2];
    }

    return err;                                    
}

/**
 * @brief Evaluate the potentials alpha and phi at a given psi value
 *
 * The values are stored in the given array as:
 * - alpha_nm[0] = alpha
 * - phi_nm[0] = phi
 *
 * @param alpha alpha value at the position where psi is evaluated
 * @param phi phi value at the position where psi is evaluated
 * @param psi psi value for which the potentials are evaluated
 * @param mode mode index to evaluate
 * @param mhddata pointer to mhd data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_stat_eval_potentials(real *alpha, real *phi, real psi, int mode,
                                mhd_stat_data* mhddata) {
    a5err err = 0;

    if(psi )

    err += interp1Dcomp_eval_f(alpha, &(mhddata->alpha_nm[mode]), psi);
    err += interp1Dcomp_eval_f(phi, &(mhddata->phi_nm[mode]), psi);

    return err;
}
