/**
 * Implements neutral_arbitrary.h.
 */
#include "defines.h"
#include "neutral.h"
#include "utils/interp.h"
#include "utils/mathlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int NeutralArbitrary_init(
    NeutralArbitrary *neutral, size_t n, size_t nr, size_t nphi, size_t nz,
    real rlim[2], real philim[2], real zlim[2], int *anum, int *znum,
    real *density, real *temperature)
{
    neutral->n = n;
    neutral->anum = (int *)malloc(n * sizeof(int));
    neutral->znum = (int *)malloc(n * sizeof(int));
    neutral->n0 = (Linear3D *)malloc(n * sizeof(Linear3D));
    neutral->T0 = (Linear3D *)malloc(n * sizeof(Linear3D));
    int err = 0;
    for (size_t i = 0; i < neutral->n; i++)
    {
        neutral->anum[i] = anum[i];
        neutral->znum[i] = znum[i];

        real *c = (real *)malloc(nr * nphi * nz * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nr * nphi * nz; i++)
        {
            c[i] = density[i];
        }
        linint3D_init(
            &neutral->n0[i], c, nr, nphi, nz, NATURALBC, PERIODICBC, NATURALBC,
            rlim[0], rlim[1], philim[0], philim[1], zlim[0], zlim[1]);
        c = (real *)malloc(nr * nphi * nz * sizeof(real));
        err += c == NULL ? 1 : 0;
        for (size_t i = 0; i < nr * nphi * nz; i++)
        {
            c[i] = temperature[i];
        }
        linint3D_init(
            &neutral->T0[i], c, nr, nphi, nz, NATURALBC, PERIODICBC, NATURALBC,
            rlim[0], rlim[1], philim[0], philim[1], zlim[0], zlim[1]);
    }
    return err;
}

void NeutralArbitrary_free(NeutralArbitrary *neutral)
{
    free(neutral->anum);
    free(neutral->znum);
    for (size_t i = 0; i < neutral->n; i++)
    {
        free(neutral->n0[i].c);
        free(neutral->T0[i].c);
    }
    free(neutral->n0);
    free(neutral->T0);
}

void NeutralArbitrary_offload(NeutralArbitrary *neutral)
{
    SUPPRESS_UNUSED_WARNING(neutral);
}

err_t NeutralArbitrary_eval_n0(
    real *n0, real r, real phi, real z, NeutralArbitrary *neutral)
{
    err_t err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr += linint3D_eval_f(&n0[i], &neutral->n0[i], r, phi, z);
    }

    if (interperr)
    {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_3D);
    }

    return err;
}

err_t NeutralArbitrary_eval_T0(
    real *T0, real r, real phi, real z, NeutralArbitrary *neutral)
{
    err_t err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for (size_t i = 0; i < neutral->n; i++)
    {
        interperr += linint3D_eval_f(&T0[i], &neutral->T0[i], r, phi, z);
    }

    if (interperr)
    {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_3D);
    }

    return err;
}
