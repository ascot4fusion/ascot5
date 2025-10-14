/**
 * Implements mhd_dynamic.h.
 */
#include "mhd_dynamic.h"
#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <stdlib.h>

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
            &mhd->alpha[i], &alpha[i * nrho * ntime], nrho, ntime, NATURALBC,
            NATURALBC, rholim[0], rholim[1], tlim[0], tlim[1]);
        err += interp2Dcomp_setup(
            &mhd->phi[i], &phi[i * nrho * ntime], nrho, ntime, NATURALBC,
            NATURALBC, rholim[0], rholim[1], tlim[0], tlim[1]);
    }
    return err;
}

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

void MhdDynamic_offload(MhdDynamic *mhd) { SUPPRESS_UNUSED_WARNING(mhd); }

err_t MhdDynamic_eval_alpha_Phi(
    real alpha[5], real Phi[5], real r, real phi, real z, real t,
    size_t include_mode, MhdDynamic *mhd, Bfield *bfield, Boozer *boozer)
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
            continue;

        /* Get interpolated values */
        real a_da[6], phi_dphi[6];
        interperr += interp2Dcomp_eval_df(a_da, &(mhd->alpha[i]), rho[0], t);
        interperr += interp2Dcomp_eval_df(phi_dphi, &(mhd->phi[i]), rho[0], t);

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
        alpha[0] += a_da[0] * mhd->amplitude[i] * cosmhd;
        Phi[0] += phi_dphi[0] * mhd->amplitude[i] * cosmhd;

        /* Time derivatives */
        alpha[1] += a_da[0] * mhd->amplitude[i] * mhd->omega[i] * sinmhd +
                    a_da[2] * mhd->amplitude[i] * cosmhd;
        Phi[1] += phi_dphi[0] * mhd->amplitude[i] * mhd->omega[i] * sinmhd +
                  phi_dphi[2] * mhd->amplitude[i] * cosmhd;

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
        err, interperr, ERR_INTERPOLATED_OUTSIDE_RANGE, DATA_MHD_DYNAMIC_C);

    return err;
}
