/**
 * Implements mhd_dynamic.h.
 */
#include "mhd_dynamic.h"
#include "defines.h"
#include "bfield.h"
#include "boozer.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include "mhd.h"
#include <stdlib.h>

/**
 * @brief Load MHD data
 *
 * @param data pointer to the data struct
 * @param nmode number of modes
 * @param nrho number of rho grid points
 * @param ntime number of time grid points
 * @param rhomin minimum rho value in the grid
 * @param rhomax maximum rho value in the grid
 * @param tmin minimum time value in the grid
 * @param tmax maximum time value in the grid
 * @param moden toroidal mode numbers
 * @param modem poloidal mode numbers
 * @param amplitude_nm amplitude of each mode
 * @param omega_nm toroidal frequency of each mode [rad/s]
 * @param phase_nm phase of each mode [rad]
 * @param alpha magnetic perturbation eigenfunctions
 * @param phi electric perturbation eigenfunctions
 *
 * @return zero if initialization succeeded.
 */
int MhdDynamic_init(
    MhdDynamic *mhd, size_t n, size_t nrho, size_t ntime, int nmode[n],
    int mmode[n], real rholim[2], real tlim[2], real amplitude[n],
    real omega[n], real phase[n], real alpha[n * nrho * ntime],
    real phi[n * nrho * ntime])
{

    int err = 0;
    mhd->n = n;
    mhd->nmode = (int *)malloc(n * sizeof(int));
    mhd->mmode = (int *)malloc(n * sizeof(int));
    mhd->omega = (real *)malloc(n * sizeof(real));
    mhd->phase = (real *)malloc(n * sizeof(real));
    mhd->amplitude = (real *)malloc(n * sizeof(real));
    mhd->phi = (Spline2D *)malloc(n * sizeof(Spline2D));
    mhd->alpha = (Spline2D *)malloc(n * sizeof(Spline2D));
    for (size_t i = 0; i < n; i++)
    {
        mhd->nmode[i] = nmode[i];
        mhd->mmode[i] = mmode[i];
        mhd->omega[i] = omega[i];
        mhd->phase[i] = phase[i];
        mhd->amplitude[i] = amplitude[i];

        err += interp2Dcomp_setup(
            &mhd->alpha[i], &alpha[i * nrho * ntime], nrho, ntime,
            NATURALBC, NATURALBC, rholim[0], rholim[1], tlim[0], tlim[1]);
        err += interp2Dcomp_setup(
            &mhd->phi[i], &phi[i * nrho * ntime], nrho, ntime, NATURALBC,
            NATURALBC, rholim[0], rholim[1], tlim[0], tlim[1]);
    }
    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void MhdDynamic_free(MhdDynamic *mhd)
{
    for (size_t i = 0; i < mhd->n; i++)
    {
        free(mhd->phi[i].c);
        free(mhd->alpha[i].c);
    }
    free(mhd->nmode);
    free(mhd->mmode);
    free(mhd->phase);
    free(mhd->omega);
    free(mhd->amplitude);
    free(mhd->phi);
    free(mhd->alpha);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void MhdDynamic_offload(MhdDynamic *mhd)
{
    SUPPRESS_UNUSED_WARNING(mhd);
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
 * @param boozer pointer to boozer data
 * @param mhd pointer to mhd data
 * @param bfield pointer to magnetic field data
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
err_t MhdDynamic_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, size_t includemode,
    Boozer *boozer, MhdDynamic *mhd, Bfield *bfield)
{

    err_t err = 0;

    real ptz[12];
    int isinside;
    if (!err)
    {
        err =
            boozer_eval_psithetazeta(ptz, &isinside, r, phi, z, bfield, boozer);
    }
    real rho[2];
    if (!err && isinside)
    {
        err = Bfield_eval_rho(rho, ptz[0], bfield);
    }

    /* Initialize values */
    for (int i = 0; i < 10; i++)
    {
        mhd_dmhd[i] = 0;
    }

    int interperr = 0;
    for (size_t i = 0; i < mhd->n; i++)
    {
        if (includemode != MHD_INCLUDE_ALL && includemode != i)
        {
            continue;
        }
        /* Get interpolated values */
        real a_da[6], phi_dphi[6];
        interperr += interp2Dcomp_eval_df(a_da, &(mhd->alpha[i]), rho[0], t);
        interperr +=
            interp2Dcomp_eval_df(phi_dphi, &(mhd->phi[i]), rho[0], t);

        /* The interpolation returns dx/drho but we require dx/dpsi. Second
           order derivatives are not needed. */
        a_da[1] *= rho[1];
        phi_dphi[1] *= rho[1];

        /* These are used frequently, so store them in separate variables */
        real mhdarg = mhd->nmode[i] * ptz[8] - mhd->mmode[i] * ptz[4] -
                      mhd->omega[i] * t + mhd->phase[i];
        real sinmhd = sin(mhdarg);
        real cosmhd = cos(mhdarg);

        /* Sum over modes to get alpha, phi */
        mhd_dmhd[0] += a_da[0] * mhd->amplitude[i] * cosmhd;
        mhd_dmhd[5] += phi_dphi[0] * mhd->amplitude[i] * cosmhd;

        /* Time derivatives */
        mhd_dmhd[1] +=
            a_da[0] * mhd->amplitude[i] * mhd->omega[i] * sinmhd +
            a_da[2] * mhd->amplitude[i] * cosmhd;
        mhd_dmhd[6] +=
            phi_dphi[0] * mhd->amplitude[i] * mhd->omega[i] * sinmhd +
            phi_dphi[2] * mhd->amplitude[i] * cosmhd;

        /* R component of gradients */
        mhd_dmhd[2] +=
            mhd->amplitude[i] * (a_da[1] * ptz[1] * cosmhd +
                                    a_da[0] * mhd->mmode[i] * ptz[5] * sinmhd -
                                    a_da[0] * mhd->nmode[i] * ptz[9] * sinmhd);
        mhd_dmhd[7] += mhd->amplitude[i] *
                       (phi_dphi[1] * ptz[1] * cosmhd +
                        phi_dphi[0] * mhd->mmode[i] * ptz[5] * sinmhd -
                        phi_dphi[0] * mhd->nmode[i] * ptz[9] * sinmhd);

        /* phi component of gradients */
        mhd_dmhd[3] += (1 / r) * mhd->amplitude[i] *
                       (a_da[1] * ptz[2] * cosmhd +
                        a_da[0] * mhd->mmode[i] * ptz[6] * sinmhd -
                        a_da[0] * mhd->nmode[i] * ptz[10] * sinmhd);
        mhd_dmhd[8] += (1 / r) * mhd->amplitude[i] *
                       (phi_dphi[1] * ptz[2] * cosmhd +
                        phi_dphi[0] * mhd->mmode[i] * ptz[6] * sinmhd -
                        phi_dphi[0] * mhd->nmode[i] * ptz[10] * sinmhd);

        /* z component of gradients */
        mhd_dmhd[4] +=
            mhd->amplitude[i] * (a_da[1] * ptz[3] * cosmhd +
                                    a_da[0] * mhd->mmode[i] * ptz[7] * sinmhd -
                                    a_da[0] * mhd->nmode[i] * ptz[11] * sinmhd);
        mhd_dmhd[9] += mhd->amplitude[i] *
                       (phi_dphi[1] * ptz[3] * cosmhd +
                        phi_dphi[0] * mhd->mmode[i] * ptz[7] * sinmhd -
                        phi_dphi[0] * mhd->nmode[i] * ptz[11] * sinmhd);
    }

    /* Omit evaluation if point outside the boozer or mhd grid. */
    if (!isinside || interperr)
    {
        for (int i = 0; i < 10; i++)
        {
            mhd_dmhd[i] = 0;
        }
    }

    return err;
}

/**
 * @brief Evaluate mhd perturbed fields Btilde, Etilde and potential
 * Phi for full orbit
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
 * @param boozer pointer to boozer data
 * @param mhd pointer to mhd data
 * @param bfield pointer to magnetic field data
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
err_t MhdDynamic_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    size_t includemode, Boozer *boozer, MhdDynamic *mhd, Bfield *bfield)
{
    err_t err = 0;
    real mhd_dmhd[10];
    if (!err)
    {
        err = MhdDynamic_eval(
            mhd_dmhd, r, phi, z, t, includemode, boozer, mhd, bfield);
    }
    /*  see example of curl evaluation in step_gc_rk4.c, ydot_gc*/
    real B_dB[15];
    if (!err)
    {
        err = Bfield_eval_b_db(B_dB, r, phi, z, t, bfield);
    }

    if (!err)
    {
        real B[3];
        B[0] = B_dB[0];
        B[1] = B_dB[4];
        B[2] = B_dB[8];

        real curlB[3];
        curlB[0] = B_dB[10] / r - B_dB[7];
        curlB[1] = B_dB[3] - B_dB[9];
        curlB[2] = (B[1] - B_dB[2]) / r + B_dB[5];

        real gradalpha[3];
        gradalpha[0] = mhd_dmhd[2];
        gradalpha[1] = mhd_dmhd[3];
        gradalpha[2] = mhd_dmhd[4];

        real gradalphacrossB[3];

        math_cross(gradalpha, B, gradalphacrossB);

        pert_field[0] = mhd_dmhd[0] * curlB[0] + gradalphacrossB[0];
        pert_field[1] = mhd_dmhd[0] * curlB[1] + gradalphacrossB[1];
        pert_field[2] = mhd_dmhd[0] * curlB[2] + gradalphacrossB[2];

        pert_field[3] = -mhd_dmhd[7] - B[0] * mhd_dmhd[1];
        pert_field[4] = -mhd_dmhd[8] - B[1] * mhd_dmhd[1];
        pert_field[5] = -mhd_dmhd[9] - B[2] * mhd_dmhd[1];
        pert_field[6] = mhd_dmhd[5];

        if (!pertonly)
        {
            pert_field[0] += B[0];
            pert_field[1] += B[1];
            pert_field[2] += B[2];
        }
    }

    return err;
}
