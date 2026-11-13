/**
 * Implements "atomic.h".
 */
#include "data/atomic.h"
#include "consts.h"
#include "defines.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "utils/ratecoefficients.h"
#include <math.h>
#include <stdlib.h>

void Atomic_offload(Atomic *atomic) { SUPPRESS_UNUSED_WARNING(atomic); }

err_t Atomic_eval_cx(
    real *ratecoeff, int charge, real ekin, size_t nspecies,
    const int *znum, const int *anum, real T0, real *n0, int extrapolate,
    Atomic *atomic)
{
    err_t err = 0;

    *ratecoeff = 0;
    if (atomic->CX != NULL)
    {
        real sigmav;
        int interperr = 0;
        for (size_t i = 0; i < nspecies; i++)
        {
            AtomicInterpolator *data = &atomic->CX[i];
            real x = data->logscale[0] ? log10(T0) : T0;
            real y = data->logscale[1] ? log10(ekin) : ekin;
            interperr += Linear2D_eval_f(&sigmav, &data->linear2d, x, y);
            *ratecoeff += exp(sigmav) * n0[i];
        }

        *ratecoeff = interperr && extrapolate ? 0. : *ratecoeff;
        err = ERROR_CHECK(
            err, interperr && !extrapolate, ERR_INTERPOLATED_OUTSIDE_RANGE,
            DATA_ATOMIC_C);
    }
    else
    {
        (void)charge;
        (void)znum;
        (void)anum;
        real sigmav;
        lindsaystebbings_sigmav(ekin, T0, &sigmav);
        *ratecoeff += sigmav * n0[0];
    }

    return err;
}

err_t Atomic_eval_bms(
    real ratecoeff[1], real ekin, real mass, size_t nion,
    const int *znum, const int *anum, real Te, real ne, real *ni,
    int extrapolate, Atomic *atomic)
{
    err_t err = 0;
    real mass_amu = mass / CONST_U;

    *ratecoeff = 0;
    if (atomic->BMS != NULL)
    {
        real sigmav;
        int interperr = 0;
        for (size_t i = 0; i < nion; i++)
        {
            AtomicInterpolator *data = &atomic->BMS[i];
            real x = data->logscale[0] ? log10(Te) : Te;
            real y = data->logscale[1] ? log10(ekin / mass_amu) : ekin / mass_amu;
            real z = data->logscale[2] ? log10(ni[i] * znum[i]) : ni[i] * znum[i];
            interperr += Linear3D_eval_f(&sigmav, &data->linear3d, x, y, z);
            *ratecoeff += exp(sigmav) * ni[i] * znum[i];
        }
        *ratecoeff /= ne;

        *ratecoeff = interperr && extrapolate ? 0. : *ratecoeff;
        err = ERROR_CHECK(
            err, interperr && !extrapolate, ERR_INTERPOLATED_OUTSIDE_RANGE,
            DATA_ATOMIC_C);
    }
    else
    {
        real gamma = physlib_gamma_Ekin(mass, ekin);
        real vnorm = physlib_vnorm_gamma(gamma);
        err = suzuki_sigma(ratecoeff, ekin / mass_amu, ne, Te, nion, ni, anum, znum);
        *ratecoeff *= vnorm;
    }

    return err;
}

err_t Atomic_eval_ei(
    real *ratecoeff, int z_1, int a_1, int q_1, real Te,
    real ne, int extrapolate, Atomic *atomic)
{

    err_t err = 0;

    *ratecoeff = 0;
    if (atomic->EI != NULL)
    {
        real sigmav;
        int interperr = 0;
        AtomicInterpolator *data = atomic->EI;
        real x = data->logscale[0] ? log10(Te) : Te;
        interperr += Linear1D_eval_f(&sigmav, &data->linear1d, x);
        *ratecoeff += exp(sigmav) * ne;

        *ratecoeff = interperr && extrapolate ? 0. : *ratecoeff;
        err = ERROR_CHECK(
            err, interperr && !extrapolate, ERR_INTERPOLATED_OUTSIDE_RANGE,
            DATA_ATOMIC_C);
    }
    else
    {
        (void)q_1;
        real sigmav = 0;
        err = lotzformula_sigmav(a_1, z_1, Te, &sigmav);
        *ratecoeff += sigmav * ne;
    }

    return err;
}

err_t Atomic_eval_ii(
    real *ratecoeff, int charge, real ekin, size_t nspecies,
    const int *znum, const int *anum, real Ti, real *ni, int extrapolate,
    Atomic *atomic)
{
    err_t err = 0;

    *ratecoeff = 0;
    if (atomic->II != NULL)
    {
        real sigmav;
        int interperr = 0;
        for (size_t i = 0; i < nspecies; i++)
        {
            AtomicInterpolator *data = &atomic->II[i];
            real x = data->logscale[0] ? log10(Ti) : Ti;
            real y = data->logscale[1] ? log10(ekin) : ekin;
            interperr += Linear2D_eval_f(&sigmav, &data->linear2d, x, y);
            *ratecoeff += exp(sigmav) * ni[i];
        }

        *ratecoeff = interperr && extrapolate ? 0. : *ratecoeff;
        err = ERROR_CHECK(
            err, interperr && !extrapolate, ERR_INTERPOLATED_OUTSIDE_RANGE,
            DATA_ATOMIC_C);
    }
    else
    {
        (void)charge;
        (void)znum;
        (void)anum;
        real sigmav;
        barnett_sigmav(ekin, Ti, &sigmav);
        *ratecoeff += sigmav * ni[0];
    }

    return err;
}
