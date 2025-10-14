/**
 * Implements mhd_stationary.h.
 */
#include "mhd_stationary.h"
#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <stdlib.h>

int MhdStationary_init(
    MhdStationary *mhd, size_t n, size_t nrho, int nmode[n], int mmode[n],
    real rholim[2], real amplitude[n], real omega[n], real phase[n],
    real alpha[n * nrho], real phi[n * nrho])
{

    int err = 0;
    mhd->n = n;
    mhd->nmode = (int *)malloc(n * sizeof(int));
    mhd->mmode = (int *)malloc(n * sizeof(int));
    mhd->omega = (real *)malloc(n * sizeof(real));
    mhd->phase = (real *)malloc(n * sizeof(real));
    mhd->amplitude = (real *)malloc(n * sizeof(real));
    mhd->phi = (Spline1D *)malloc(n * sizeof(Spline1D));
    mhd->alpha = (Spline1D *)malloc(n * sizeof(Spline1D));
    for (size_t i = 0; i < n; i++)
    {
        mhd->nmode[i] = nmode[i];
        mhd->mmode[i] = mmode[i];
        mhd->omega[i] = omega[i];
        mhd->phase[i] = phase[i];
        mhd->amplitude[i] = amplitude[i];

        err += interp1Dcomp_setup(
            &mhd->alpha[i], &alpha[i * nrho], nrho, NATURALBC, rholim[0],
            rholim[1]);
        err += interp1Dcomp_setup(
            &mhd->phi[i], &phi[i * nrho], nrho, NATURALBC, rholim[0],
            rholim[1]);
    }
    return err;
}

void MhdStationary_free(MhdStationary *mhd)
{
    for (size_t i = 0; i < mhd->n; i++)
    {
        free(mhd->phi[i].c);
        free(mhd->alpha[i].c);
    }
    free(mhd->phi);
    free(mhd->alpha);
    free(mhd->nmode);
    free(mhd->mmode);
    free(mhd->phase);
    free(mhd->omega);
    free(mhd->amplitude);
}

void MhdStationary_offload(MhdStationary *mhd) { SUPPRESS_UNUSED_WARNING(mhd); }

err_t MhdStationary_eval_alpha_Phi(
    real alpha[5], real Phi[5], real r, real phi, real z, real t,
    size_t include_mode, MhdStationary *mhd, Bfield *bfield, Boozer *boozer)
{
    err_t err = 0;
    real ptz[12], rho[2];
    int isinside;
    err = Boozer_map_coordinates(ptz, &isinside, r, phi, z, t, boozer, bfield);
    err = err ? err : Bfield_eval_rho(rho, ptz[0], bfield);

    for (size_t i = 0; i < 5; i++)
    {
        Phi[i] = 0;
        alpha[i] = 0;
    }

    int interperr = 0;
    for (size_t i = 0; i < mhd->n; i++)
    {
        if (include_mode != MHD_INCLUDE_ALL && include_mode != i)
        {
            continue;
        }

        real a_da[3], phi_dphi[3];
        interperr += interp1Dcomp_eval_df(a_da, &(mhd->alpha[i]), rho[0]);
        interperr += interp1Dcomp_eval_df(phi_dphi, &(mhd->phi[i]), rho[0]);

        /* The interpolation returns dx/drho but we require dx/dpsi. The second
           order derivatives are not needed anywhere */
        a_da[1] *= rho[1];
        phi_dphi[1] *= rho[1];

        /* Frequently used helper variables */
        real mhdarg = mhd->nmode[i] * ptz[8] - mhd->mmode[i] * ptz[4] -
                      mhd->omega[i] * t + mhd->phase[i];
        real sinmhd = sin(mhdarg);
        real cosmhd = cos(mhdarg);

        alpha[0] += a_da[0] * mhd->amplitude[i] * cosmhd;
        Phi[0] += phi_dphi[0] * mhd->amplitude[i] * cosmhd;

        /* Time derivatives */
        alpha[1] += a_da[0] * mhd->amplitude[i] * mhd->omega[i] * sinmhd;
        Phi[1] += phi_dphi[0] * mhd->amplitude[i] * mhd->omega[i] * sinmhd;

        /* R component of gradients */
        alpha[2] +=
            mhd->amplitude[i] * (a_da[1] * ptz[1] * cosmhd +
                                 a_da[0] * mhd->mmode[i] * ptz[5] * sinmhd -
                                 a_da[0] * mhd->nmode[i] * ptz[9] * sinmhd);
        Phi[2] +=
            mhd->amplitude[i] * (phi_dphi[1] * ptz[1] * cosmhd +
                                 phi_dphi[0] * mhd->mmode[i] * ptz[5] * sinmhd -
                                 phi_dphi[0] * mhd->nmode[i] * ptz[9] * sinmhd);

        /* phi component of gradients */
        alpha[3] += (1 / r) * mhd->amplitude[i] *
                    (a_da[1] * ptz[2] * cosmhd +
                     a_da[0] * mhd->mmode[i] * ptz[6] * sinmhd -
                     a_da[0] * mhd->nmode[i] * ptz[10] * sinmhd);
        Phi[3] += (1 / r) * mhd->amplitude[i] *
                  (phi_dphi[1] * ptz[2] * cosmhd +
                   phi_dphi[0] * mhd->mmode[i] * ptz[6] * sinmhd -
                   phi_dphi[0] * mhd->nmode[i] * ptz[10] * sinmhd);

        /* z component of gradients */
        alpha[4] +=
            mhd->amplitude[i] * (a_da[1] * ptz[3] * cosmhd +
                                 a_da[0] * mhd->mmode[i] * ptz[7] * sinmhd -
                                 a_da[0] * mhd->nmode[i] * ptz[11] * sinmhd);
        Phi[4] += mhd->amplitude[i] *
                  (phi_dphi[1] * ptz[3] * cosmhd +
                   phi_dphi[0] * mhd->mmode[i] * ptz[7] * sinmhd -
                   phi_dphi[0] * mhd->nmode[i] * ptz[11] * sinmhd);
    }

    for (size_t i = 0; i < 5; i++)
    {
        Phi[i] = isinside ? Phi[i] : 0.0;
        alpha[i] = isinside ? alpha[i] : 0.0;
    }

    err = ERROR_CHECK(
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_MHD_STATIONARY_C);
    return err;
}
