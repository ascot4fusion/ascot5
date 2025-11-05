/**
 * Implements atomic.h.
 */
#include "data/atomic.h"
#include "consts.h"
#include "defines.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "utils/suzuki.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#pragma omp declare target
/** Set values outside abscissae to zero instead of raising an error. */
static int ASIGMA_EXTRAPOLATE = 0;
#pragma omp end declare target

int Atomic_init(
    Atomic *atomic, size_t nreac, int *z1, int *a1, int *z2, int *a2,
    int *reactype, int *ne, real *emin, real *emax, int *nn, real *nmin,
    real *nmax, int *nT, real *Tmin, real *Tmax, real *sigma)
{

    int err = 0;
    atomic->N_reac = nreac;
    atomic->z_1 = (int *)malloc(nreac * sizeof(int));
    atomic->a_1 = (int *)malloc(nreac * sizeof(int));
    atomic->z_2 = (int *)malloc(nreac * sizeof(int));
    atomic->a_2 = (int *)malloc(nreac * sizeof(int));
    atomic->reac_type =
        (asigma_reac_type *)malloc(nreac * sizeof(asigma_reac_type));
    atomic->sigma = (Spline1D *)malloc(nreac * sizeof(Spline1D));
    atomic->sigmav = (Spline2D *)malloc(nreac * sizeof(Spline2D));
    atomic->BMSsigmav = (Spline3D *)malloc(nreac * sizeof(Spline3D));
    for (size_t i_reac = 0; i_reac < nreac; i_reac++)
    {
        atomic->z_1[i_reac] = z1[i_reac];
        atomic->a_1[i_reac] = a1[i_reac];
        atomic->z_2[i_reac] = z2[i_reac];
        atomic->a_2[i_reac] = a2[i_reac];
        atomic->reac_type[i_reac] = reactype[i_reac];

        /* Initialize spline struct according to dimensionality of
           reaction data (and mark reaction availability) */
        int dim = (ne[i_reac] > 1) + (nn[i_reac] > 1) + (nT[i_reac] > 1);
        real *pos = sigma;
        switch (dim)
        {
        case 1:
            err = Spline1D_init(
                &atomic->sigma[i_reac], ne[i_reac], NATURALBC,
                (real[]){emin[i_reac], emax[i_reac]}, pos);
            pos += ne[i_reac];
            break;
        case 2:
            err = Spline2D_init(
                &atomic->sigmav[i_reac], ne[i_reac], nT[i_reac], NATURALBC,
                NATURALBC, (real[]){emin[i_reac], emax[i_reac]},
                (real[]){Tmin[i_reac], Tmax[i_reac]}, pos);
            pos += ne[i_reac] * nT[i_reac];
            break;
        case 3:
            err = Spline3D_init(
                &atomic->BMSsigmav[i_reac], ne[i_reac], nn[i_reac], nT[i_reac],
                NATURALBC, NATURALBC, NATURALBC,
                (real[]){emin[i_reac], emax[i_reac]},
                (real[]){nmin[i_reac], nmax[i_reac]},
                (real[]){Tmin[i_reac], Tmax[i_reac]}, pos);
            pos += ne[i_reac] * nn[i_reac] * nT[i_reac];
            break;
        default:
            err = 1;
            break;
        }
    }

    return err;
}

void Atomic_free(Atomic *atomic)
{
    for (size_t i_reac = 0; i_reac < atomic->N_reac; i_reac++)
    {
        if (atomic->reac_type[i_reac] == sigma_CX)
        {
            free(atomic->sigma[i_reac].c);
        }
        else if (atomic->reac_type[i_reac] == sigmav_CX)
        {
            free(atomic->sigmav[i_reac].c);
        }
        else if (atomic->reac_type[i_reac] == sigmav_BMS)
        {
            free(atomic->BMSsigmav[i_reac].c);
        }
    }
    free(atomic->z_1);
    free(atomic->a_1);
    free(atomic->z_2);
    free(atomic->a_2);
    free(atomic->reac_type);
    free(atomic->sigma);
    free(atomic->sigmav);
    free(atomic->BMSsigmav);
}

void Atomic_offload(Atomic *data) { SUPPRESS_UNUSED_WARNING(data); }

void Atomic_extrapolate(int extrapolate) { ASIGMA_EXTRAPOLATE = extrapolate; }

err_t Atomic_eval_cx(
    real *ratecoeff, int z_1, int a_1, real E, real mass, size_t nspec,
    const int *znum, const int *anum, real T_0, real *n_0, Atomic *atomic)
{
    err_t err = 0;
    (void)mass;

    /* Convert Joule to eV */
    E /= CONST_E;
    T_0 /= CONST_E;
    *ratecoeff = 0;
    for (size_t i_spec = 0; i_spec < nspec; i_spec++)
    {

        /* Find the matching reaction */
        size_t i_reac;
        int reac_found = -1;
        for (i_reac = 0; i_reac < atomic->N_reac; i_reac++)
        {
            if (atomic->z_1[i_reac] == z_1 && atomic->a_1[i_reac] == a_1 &&
                atomic->z_2[i_reac] == znum[i_spec] &&
                atomic->a_2[i_reac] == anum[i_spec] &&
                atomic->reac_type[i_reac] == sigmav_CX)
            {
                reac_found = i_reac;
            }
        }
        i_reac = reac_found;

        if (reac_found < 0)
        {
            /* Reaction not found. Raise error. */
            // err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
        }
        else
        {
            real sigmav;
            int interperr =
                Spline2D_eval_f(&sigmav, &atomic->sigmav[i_reac], E, T_0);

            /* Interpolation error means the data has to be extrapolated */
            if (interperr)
            {
                if (ASIGMA_EXTRAPOLATE)
                {
                    sigmav = 0.0;
                }
                else
                {
                    // err = error_raise(
                    //     ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
                }
            }
            *ratecoeff += sigmav * n_0[i_spec];
        }
    }

    return err;
}

err_t Atomic_eval_bms(
    real *ratecoeff, int z_1, int a_1, real E, real mass, size_t nion,
    const int *znum, const int *anum, real T_e, real *n_i, Atomic *atomic)
{
    err_t err = 0;

    /* Convert Joule to eV */
    real E_eV = E / CONST_E;
    T_e /= CONST_E;

    /* Find the matching reaction. Note that BMS data is same for all
     * isotopes, so we don't compare anums */
    int reac_found = -1;
    real n_e = 0;
    *ratecoeff = 0;
    for (size_t i_spec = 0; i_spec < nion; i_spec++)
    {
        n_e += znum[i_spec] * n_i[i_spec];
        for (size_t i_reac = 0; i_reac < atomic->N_reac; i_reac++)
        {
            if (atomic->z_1[i_reac] == z_1 &&
                atomic->z_2[i_reac] == znum[i_spec] &&
                atomic->reac_type[i_reac] == sigmav_BMS)
            {
                reac_found = i_reac;
                real sigmav;
                int interperr = Spline3D_eval_f(
                    &sigmav, &atomic->BMSsigmav[i_reac], E_eV / a_1,
                    znum[i_spec] * n_i[i_spec], T_e);

                /* Interpolation error means the data has to be extrapolated */
                if (interperr)
                {
                    if (ASIGMA_EXTRAPOLATE)
                    {
                        sigmav = 0.0;
                    }
                    else
                    {
                        // err = error_raise(
                        //     ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
                    }
                }
                *ratecoeff += sigmav * (znum[i_spec] * n_i[i_spec]);
            }
        }
    }
    *ratecoeff /= n_e;

    if (reac_found < 0)
    {
        /* Reaction not found. Try Suzuki before throwing error. */
        T_e *= CONST_E;
        real gamma = physlib_gamma_Ekin(mass, E);
        real vnorm = physlib_vnorm_gamma(gamma);
        if (suzuki_sigmav(
                ratecoeff, E / a_1, vnorm, n_e, T_e, nion, n_i, anum, znum))
        {
            // err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
        }
    }

    return err;
}
